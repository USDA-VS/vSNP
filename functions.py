import os
import sys
import shutil
import subprocess
import gzip
import glob
import csv
import json
import time
import regex
import re
import numpy as np
import pandas as pd
import zipfile
import xlsxwriter
import xlrd
import pysam
import vcf
import smtplib
from multiprocessing import Pool
from prettytable import PrettyTable
from dask import delayed
from itertools import repeat as itertools_repeat
from numpy import mean
from email.utils import formatdate
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders
from distutils.dir_util import copy_tree
from datetime import datetime
from concurrent import futures
from collections import OrderedDict
from collections import Counter
from collections import defaultdict
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
import logging
import inspect

from parameters import Get_Specie_Parameters

# logging.basicConfig(format='%(levelname)s: %(message)s', filemode='w', filename='debug.log', level=logging.DEBUG)
# logging.getLogger().addHandler(logging.StreamHandler()) #print to console

# def warning_log(ex, inspect_getframeinfo, *args):
#     logging.warning(f'\nException occured, file: {inspect_getframeinfo.filename}\nfuction: {inspect.stack()[0][3]}, near line in script {inspect_getframeinfo.lineno} --> {type(ex).__name__, ex.args}\nAdditional args: {args}\n\n')

# def debug_log(ex, inspect_getframeinfo, *args):
#     logging.debug(f'\nException occured, file: {inspect_getframeinfo.filename}\nfuction: {inspect.stack()[0][3]}, near line in script {inspect_getframeinfo.lineno} --> {type(ex).__name__, ex.args}\nAdditional args: {args}\n\n')


def run_loop(arg_options):

    root_dir = arg_options['root_dir']
    limited_cpu_count = arg_options['limited_cpu_count']

    startTime = datetime.now()
    ts = time.time()
    st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
    print("Start time: %s" % st)

    list_of_files = glob.glob('*gz')

    for afile in list_of_files:
        prefix_name = re.sub('_.*', '', afile)
        prefix_name = re.sub('\..*', '', prefix_name)
        print(prefix_name)
        if not os.path.exists(prefix_name):
            os.makedirs(prefix_name)
        shutil.move(afile, prefix_name)

    # placed at root
    # get file opened and give a header
    summary_file = root_dir + '/stat_alignment_summary_' + st + '.xlsx'
    workbook = xlsxwriter.Workbook(summary_file)
    worksheet = workbook.add_worksheet()
    row = 0
    col = 0
    top_row_header = ["time_stamp", "sample_name", "species", "reference_sequence_name", "R1size", "R2size", "Q_ave_R1", "Q_ave_R2", "Q30_R1", "Q30_R2", "allbam_mapped_reads", "genome_coverage", "ave_coverage", "ave_read_length", "unmapped_reads", "unmapped_assembled_contigs", "good_snp_count", "mlst_type", "octalcode", "sbcode", "hexadecimal_code", "binarycode"]
    for header in top_row_header:
        worksheet.write(row, col, header)
        col += 1
    ###

    # Cumulative stats
    path_found = False
    if os.path.isdir("/project/mycobacteria_brucella/mycobacterium/stats"): #check bioinfo from server
        path_found = True
        copy_to = "/project/mycobacteria_brucella/mycobacterium/stats"
    elif os.path.isdir("/Volumes/bioinfo/project/mycobacteria_brucella/mycobacterium/stats"): #check bioinfo from Mac
        path_found = True
        copy_to = "/Volumes/bioinfo/project/mycobacteria_brucella/mycobacterium/stats"
    else:
        copy_to = None
        print("Bioinfo not connected")

    if path_found:
        try:
            summary_cumulative_file = copy_to + '/stat_alignment_culmulative_summary' + '.xlsx'
            summary_cumulative_file_temp = copy_to + '/stat_alignment_culmulative_summary-' + st + '-temp.xlsx'
            temp_folder = copy_to + '/temp'
        except OSError as ex:
            text = "ERROR, Bioinfo unresponsive unable to copy to stats file"
            msg = MIMEMultipart()
            msg['From'] = "tod.p.stuber@aphis.usda.gov"
            msg['To'] = "tod.p.stuber@aphis.usda.gov"
            msg['Date'] = formatdate(localtime=True)
            msg['Subject'] = "### No coverage file"
            msg.attach(MIMEText(text))
            smtp = smtplib.SMTP('10.10.8.12')
            smtp.send_message(msg)
            smtp.quit()
            # debug_log(ex, inspect.getframeinfo(inspect.currentframe()), "Bioinfo unresponsive unable to copy to stats file")
    ###

    directory_list = [item for item in os.listdir() if os.path.isdir(item)]
    total_samples = len(directory_list)
    directory_size = {}
    for folder in directory_list: #run files by size, smallest to largest
        size = sum(os.path.getsize(os.path.join(dirpath,filename)) for dirpath, dirnames, filenames in os.walk(folder) for filename in filenames)
        directory_size[folder] = size
    directory_size = {k: v for k, v in sorted(directory_size.items(), key=lambda x: x[1])}
    directory_list = [*directory_size] #ordered list
    lower_count = 0
    upper_count = 1
    row = 1
    while lower_count < total_samples:
        upper_count = lower_count + limited_cpu_count
        run_list = directory_list[lower_count:upper_count] #create a run list
        for i in run_list:
            directory_list.remove(i)
        total_samples = len(directory_list)
        print(run_list)

        print("Iterating directories")
        frames = []
        if arg_options['debug_call']: #run just one sample at a time to debug
            for sample_name in run_list:
                print("DEBUGGING, SAMPLES RAN INDIVIDUALLY")
                stat_summary = read_aligner(sample_name, arg_options)

                df_stat_summary = pd.DataFrame.from_dict(stat_summary, orient='index') #convert stat_summary to df
                frames.append(df_stat_summary) #frames to concatenate
                worksheet.write(row, 0, stat_summary.get('time_stamp', 'n/a'))
                worksheet.write(row, 1, stat_summary.get('sample_name', 'n/a'))
                worksheet.write(row, 2, stat_summary.get('species', 'n/a'))
                worksheet.write(row, 3, stat_summary.get('reference_sequence_name', 'n/a'))
                worksheet.write(row, 4, stat_summary.get('R1size', 'n/a'))
                worksheet.write(row, 5, stat_summary.get('R2size', 'n/a'))
                worksheet.write(row, 6, stat_summary.get('Q_ave_R1', 'n/a'))
                worksheet.write(row, 7, stat_summary.get('Q_ave_R2', 'n/a'))
                worksheet.write(row, 8, stat_summary.get('Q30_R1', 'n/a'))
                worksheet.write(row, 9, stat_summary.get('Q30_R2', 'n/a'))
                worksheet.write(row, 10, stat_summary.get('allbam_mapped_reads', 'n/a'))
                worksheet.write(row, 11, stat_summary.get('genome_coverage', 'n/a'))
                worksheet.write(row, 12, stat_summary.get('ave_coverage', 'n/a'))
                worksheet.write(row, 13, stat_summary.get('ave_read_length', 'n/a'))
                worksheet.write(row, 14, stat_summary.get('unmapped_reads', 'n/a'))
                worksheet.write(row, 15, stat_summary.get('unmapped_assembled_contigs', 'n/a'))
                worksheet.write(row, 16, stat_summary.get('good_snp_count', 'n/a'))
                worksheet.write(row, 17, stat_summary.get('mlst_type', 'n/a'))
                worksheet.write(row, 18, stat_summary.get('octalcode', 'n/a'))
                worksheet.write(row, 19, stat_summary.get('sbcode', 'n/a'))
                worksheet.write(row, 20, stat_summary.get('hexadecimal_code', 'n/a'))
                worksheet.write(row, 21, stat_summary.get('binarycode', 'n/a'))
                row += 1
        else:  # run all in run_list in parallel
            print("SAMPLES RAN IN PARALLEL")
            # itertools allows additional arguments to pass
            # Need to speed test which why is faster
            with futures.ProcessPoolExecutor(max_workers=limited_cpu_count) as pool:
                for stat_summary in pool.map(read_aligner, run_list, itertools_repeat(arg_options)):
                    df_stat_summary = pd.DataFrame.from_dict(stat_summary, orient='index')  #convert stat_summary to df
                    frames.append(df_stat_summary)  #frames to concatenate

                    worksheet.write(row, 0, stat_summary.get('time_stamp', 'n/a'))
                    worksheet.write(row, 1, stat_summary.get('sample_name', 'n/a'))
                    worksheet.write(row, 2, stat_summary.get('species', 'n/a'))
                    worksheet.write(row, 3, stat_summary.get('reference_sequence_name', 'n/a'))
                    worksheet.write(row, 4, stat_summary.get('R1size', 'n/a'))
                    worksheet.write(row, 5, stat_summary.get('R2size', 'n/a'))
                    worksheet.write(row, 6, stat_summary.get('Q_ave_R1', 'n/a'))
                    worksheet.write(row, 7, stat_summary.get('Q_ave_R2', 'n/a'))
                    worksheet.write(row, 8, stat_summary.get('Q30_R1', 'n/a'))
                    worksheet.write(row, 9, stat_summary.get('Q30_R2', 'n/a'))
                    worksheet.write(row, 10, stat_summary.get('allbam_mapped_reads', 'n/a'))
                    worksheet.write(row, 11, stat_summary.get('genome_coverage', 'n/a'))
                    worksheet.write(row, 12, stat_summary.get('ave_coverage', 'n/a'))
                    worksheet.write(row, 13, stat_summary.get('ave_read_length', 'n/a'))
                    worksheet.write(row, 14, stat_summary.get('unmapped_reads', 'n/a'))
                    worksheet.write(row, 15, stat_summary.get('unmapped_assembled_contigs', 'n/a'))
                    worksheet.write(row, 16, stat_summary.get('good_snp_count', 'n/a'))
                    worksheet.write(row, 17, stat_summary.get('mlst_type', 'n/a'))
                    worksheet.write(row, 18, stat_summary.get('octalcode', 'n/a'))
                    worksheet.write(row, 19, stat_summary.get('sbcode', 'n/a'))
                    worksheet.write(row, 20, stat_summary.get('hexadecimal_code', 'n/a'))
                    worksheet.write(row, 21, stat_summary.get('binarycode', 'n/a'))
                    row += 1

            if not arg_options['quiet'] and path_found:
                try:
                    open_check = open(summary_cumulative_file, 'a') #'a' is very important, 'w' will leave you with an empty file
                    open_check.close()
                    df_all = pd.read_excel(summary_cumulative_file)
                    df_all_trans = df_all.T #indexed on column headers
                    # save back the old and remake the working stats file
                    shutil.move(summary_cumulative_file, '{}' .format(temp_folder + '/stat_backup' + st + '.xlsx'))
                    sorter = list(df_all_trans.index) #list of original column order
                    frames.insert(0, df_all_trans) #put as first item in list
                    df_concat = pd.concat(frames, axis=1, sort=True) #cat frames
                    df_sorted = df_concat.loc[sorter] #sort based on sorter order
                    df_sorted.T.to_excel(summary_cumulative_file, index=False) # transpose before writing to excel, numerical index not needed
                except BlockingIOError as ex:
                    sorter = list(df_stat_summary.index) #list of original column order
                    df_concat = pd.concat(frames, axis=1, sort=True) #cat frames
                    df_sorted = df_concat.loc[sorter] #sort based on sorter order
                    df_sorted.T.to_excel(summary_cumulative_file_temp, index=False)
                    # debug_log(ex, inspect.getframeinfo(inspect.currentframe()))
                    pass
                except OSError as ex:
                    sorter = list(df_stat_summary.index) #list of original column order
                    df_concat = pd.concat(frames, axis=1, sort=True) #cat frames
                    df_sorted = df_concat.loc[sorter] #sort based on sorter order
                    # debug_log(ex, inspect.getframeinfo(inspect.currentframe()))
                    try:
                        df_sorted.T.to_excel(summary_cumulative_file_temp, index=False)
                    except OSError as ex:
                        # debug_log(ex, inspect.getframeinfo(inspect.currentframe()), "##### UNABLE TO MAKE CONNECTION TO BIOINFO")
                        pass
                    pass
            else:
                logging.debug("Path to cumulative stat summary file not found")
                print("Path to cumulative stat summary file not found")

    workbook.close()

    runtime = (datetime.now() - startTime)
    print("\n\nruntime: %s:  \n" % runtime)

    if arg_options['email_list']:
        try:
            send_email_step1(arg_options['email_list'], runtime, path_found, summary_file, st)
        except TimeoutError as ex:
            # debug_log(ex, inspect.getframeinfo(inspect.currentframe()), "Unable to send email with current smtp setting")
            pass


def reference_table():

    pretty_table = PrettyTable(['-r option', 'Species', 'NCBI identifier'])
    pretty_table.add_row(['af', 'Mycobacterium_bovis_AF2122/97', 'NC_002945.4'])
    pretty_table.add_row(['h37', 'Mycobacterium tuberculosis H37Rv', 'NC_000962.3'])
    pretty_table.add_row(['ab1', 'Brucella abortus biovar 1 str. 9-941', 'NC_006932.1, NC_006933.1'])
    pretty_table.add_row(['ab3', 'Brucella abortus strain BER', 'NZ_CP007682.1, NZ_CP007683.1'])
    pretty_table.add_row(['suis1', 'Brucella suis 1330', 'NC_017251.1, NC_017250.1'])
    pretty_table.add_row(['suis2', 'Brucella suis ATCC 23445', 'NC_010169.1, NC_010167.1'])
    pretty_table.add_row(['suis3', 'Brucella suis bv. 3 str. 686', 'NZ_CP007719.1, NZ_CP007718.1'])
    pretty_table.add_row(['mel1', 'Brucella melitensis bv. 1 str. 16M', 'NC_003317.1, NC_003318.1'])
    pretty_table.add_row(['mel1b', 'Brucella melitensis BwIM_SOM_36b', 'NZ_CP018508.1, NZ_CP018509.1'])
    pretty_table.add_row(['mel2', 'Brucella melitensis ATCC 23457', 'NC_012441.1, NC_012442.1'])
    pretty_table.add_row(['mel3', 'Brucella melitensis bv. 3 str. Ether', 'NZ_CP007760, NZ_CP007761'])
    pretty_table.add_row(['canis', 'Brucella canis ATCC 23365', 'NC_010103.1, NC_010104.1'])
    pretty_table.add_row(['ceti1', 'Bceti1Cudo', 'Bceti1Cudo inhouse'])
    pretty_table.add_row(['ceti2', 'Brucella ceti TE10759-12', 'NC_022905.1, NC_022906.1'])
    pretty_table.add_row(['ovis', 'Brucella ovis ATCC 25840', 'NC_009505.1, NC_009504.1'])
    pretty_table.add_row(['para-NC_002944', 'Mycobacterium avium subsp. paratuberculosis str. k10', 'NC_002944.2'])
    pretty_table.add_row(['para-CP033688', 'Mycobacterium avium subsp. paratuberculosis strain Telford', 'CP033688.1'])
    pretty_table.add_row(['typhimurium-14028S', 'Salmonella enterica subsp. enterica serovar Typhimurium str. 14028S', 'NC_016856.1, NC_016855.1(plasmid)'])
    pretty_table.add_row(['typhimurium-LT2', 'Salmonella enterica subsp. enterica serovar Typhimurium str. LT2', 'AE006468.2'])
    pretty_table.add_row(['heidelberg-SL476', 'Salmonella enterica subsp. enterica serovar Heidelberg str. SL476', 'NC_011083.1'])
    pretty_table.add_row(['strepequi', 'Streptococcus equi subsp. zooepidemicus ATCC 35246', 'NC_017582'])
    pretty_table.add_row(['blockley', 'Salmonella enterica subsp. enterica serovar Blockley strain 79-1229', 'GCF_002283115.1_ASM228311v1_genomic'])
    pretty_table.add_row(['infantis-FSIS1502169', 'Salmonella enterica subsp. enterica serovar Infantis strain FSIS1502169', 'CP016406.1'])
    pretty_table.add_row(['dublin-ATCC39184', 'Salmonella enterica subsp. enterica serovar Dublin str. ATCC 39184', 'CP019179.1'])
    pretty_table.add_row(['kentucky-SA20030505', 'Salmonella enterica subsp. enterica serovar Kentucky str. SA20030505 chromosome', 'CP022500.1'])
    pretty_table.add_row(['newport-USDA-ARS-USMARC-1925', 'Salmonella enterica subsp. enterica serovar Newport str. USDA-ARS-USMARC-1925 chromosome', 'CP025232.1'])
    pretty_table.add_row(['senftenberg-NCTC10080', 'Salmonella enterica subsp. enterica serovar Senftenberg strain NCTC10080 genome assembly', 'LS483465.1'])
    pretty_table.add_row(['enteritidis-Durban', 'Salmonella enterica subsp. enterica serovar Enteritidis strain Durban', 'CP007507.1'])
    pretty_table.add_row(['montevideo-507440-20', 'Salmonella enterica subsp. enterica serovar Montevideo str. 507440-20', 'CP007530.1'])
    pretty_table.add_row(['te_atcc35865', 'Taylorella equigenitalis ATCC 35865', 'NC_018108.1'])
    pretty_table.add_row(['te_09-0932', 'Taylorella equigenitalis strain 09-0932', 'NZ_CP021201.1'])
    pretty_table.add_row(['te_89-0490', 'Taylorella equigenitalis strain 89-0490', 'NZ_CP021199.1'])
    pretty_table.add_row(['te_92-0972', 'Taylorella equigenitalis strain 92-0972', 'NZ_CP021060.1'])
    pretty_table.add_row(['te_98-0554', 'Taylorella equigenitalis strain 98-0554', 'NZ_CP021246.1'])
    pretty_table.add_row(['te_mce9', 'Taylorella equigenitalis MCE9', 'NC_014914.1'])
    pretty_table.add_row(['flu', 'H7N3', 'segments 1-8'])
    pretty_table.add_row(['newcaste', '18-016505-001-fusion-HN', '18-016505-001-fusion-HN'])
    pretty_table.add_row(['belize', 'Newcastle disease virus isolate Belize (Spanish Lookout)/4224-3/2008', 'KF767466.1'])
    pretty_table.add_row(['gua', 'Newcastle disease virus isolate Guatemala chicken/Guatemala/Gua1-1407/2018', 'Gua1_1407_2018'])

    return pretty_table


def fix_vcf(each_vcf, arg_options):
    mal = []
    # Fix common VCF errors
    temp_file = each_vcf + ".temp"
    write_out = open(temp_file, 'w') #r+ used for reading and writing to the same file
    initial_file_time_stats = os.stat(each_vcf)
    with open(each_vcf, 'r') as file:
        try:
            for line in file:
                if line.rstrip(): # true if not empty line'^$'
                    line = line.rstrip() #remove right white space
                    line = re.sub('"AC=', 'AC=', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('"$', '', line)
                    line = re.sub('GQ:PL\t"', 'GQ:PL\t', line)
                    line = re.sub('[0-9]+\tGT\t.\/.$', '999\tGT:AD:DP:GQ:PL\t1/1:0,80:80:99:2352,239,0', line)
                    line = re.sub('^"', '', line)
                    if line.startswith('##') and line.endswith('"'):
                        line = re.sub('"$', '', line)
                    if line.startswith('##'):
                        line = line.split('\t')
                        line = ''.join(line[0])
                    if not line.startswith('##'):
                        line = re.sub('"', '', line)
                        line = line.split('\t')
                        line = "\t".join(line[0:10])
                        print(line, file=write_out)
                    else:
                        print(line, file=write_out)
        except IndexError as ex:
            mal.append("##### IndexError: Deleting corrupt VCF file: " + each_vcf)
            # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), each_vcf, "Deleting corrupt VCF file")
            os.remove(each_vcf)
        except UnicodeDecodeError as ex:
            mal.append("##### UnicodeDecodeError: Deleting corrupt VCF file: " + each_vcf)
            # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), each_vcf, "Deleting corrupt VCF file")
            os.remove(each_vcf)

    write_out.close()
    os.rename(temp_file, each_vcf)
    # revert timestamp to original allows elites to properly sort on file modification time
    os.utime(each_vcf, times=(initial_file_time_stats.st_mtime, initial_file_time_stats.st_mtime))
    return mal


def read_aligner(sample_name, arg_options):

    os.chdir(arg_options['root_dir'] + "/" + sample_name)
    sample_directory = str(os.getcwd())
    pair_check = len(glob.glob('*_R2*fastq.gz'))
    if pair_check > 0:
        R1 = glob.glob('*_R1*fastq.gz')
        R2 = glob.glob('*_R2*fastq.gz')
        paired = True
    else:
        R1 = glob.glob('*fastq.gz')
        R2 = None
        paired = False
    arg_options['paired'] = paired

    if len(R1) > 1:
        print("#### Check for a duplicate file in {}" .format(sample_name))
        sys.exit(0)

    os.makedirs("zips")
    shutil.move(R1[0], "zips")
    arg_options['R1'] = sample_directory + "/zips/" + R1[0]
    if paired:
        shutil.move(R2[0], "zips")
        arg_options['R2'] = sample_directory + "/zips/" + R2[0]
    else:
        arg_options['R2'] = None

    read_quality_stats = {}
    print("Getting mean for {}" .format(arg_options['R1']))
    handle = gzip.open(arg_options['R1'], "rt")
    mean_quality_list = []
    for rec in SeqIO.parse(handle, "fastq"):
        mean_q = get_read_mean(rec)
        mean_quality_list.append(mean_q)

    read_quality_stats["Q_ave_R1"] = "{:.1f}" .format(mean(mean_quality_list))
    thirty_or_greater_count = sum(i > 29 for i in mean_quality_list)
    read_quality_stats["Q30_R1"] = "{:.1%}" .format(thirty_or_greater_count / len(mean_quality_list))

    if paired:
        print("Getting mean for {}" .format(arg_options['R2']))
        handle = gzip.open(arg_options['R2'], "rt")
        mean_quality_list = []
        for rec in SeqIO.parse(handle, "fastq"):
            mean_q = get_read_mean(rec)
            mean_quality_list.append(mean_q)

        read_quality_stats["Q_ave_R2"] = "{:.1f}" .format(mean(mean_quality_list))
        thirty_or_greater_count = sum(i > 29 for i in mean_quality_list)
        read_quality_stats["Q30_R2"] = "{:.1%}" .format(thirty_or_greater_count / len(mean_quality_list))
    arg_options['read_quality_stats'] = read_quality_stats

    arg_options['sample_name'] = sample_name
    arg_options = species_selection_step1(arg_options)
    if arg_options['species'] is None:
        arg_options = {**arg_options, **read_quality_stats}
        R1size = sizeof_fmt(os.path.getsize(arg_options['R1']))
        if paired:
            R2size = sizeof_fmt(os.path.getsize(arg_options['R2']))
        arg_options['R1size'] = R1size
        arg_options['R2size'] = R2size
        return arg_options
    try:
        stat_summary = align_reads(arg_options)
        for k, v in stat_summary.items():
            print("%s: %s" % (k, v))
        return(stat_summary)
    except AttributeError as ex:
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), sample_name, "Unable to return stat_summary")
        return arg_options


def get_species(arg_options):

    #species = corresponding NCBI accession
    species_cross_reference = {}
    species_cross_reference["salmonella"] = ["016856, 016855"]
    species_cross_reference["bovis"] = ["AF2122_NC002945", "00879"]
    species_cross_reference["af"] = ["NC_002945.4"]
    species_cross_reference["h37"] = ["000962", "002755", "009525", "018143"]
    species_cross_reference["para-NC_002944"] = ["NC_002944"]
    species_cross_reference["para-CP033688"] = ["CP033688"]
    species_cross_reference["ab1"] = ["006932", "006933"]
    species_cross_reference["ab3"] = ["007682", "007683"]
    species_cross_reference["canis"] = ["010103", "010104"]
    species_cross_reference["ceti1"] = ["Bceti1Cudo"]
    species_cross_reference["ceti2"] = ["022905", "022906"]
    species_cross_reference["mel1"] = ["003317", "003318"]
    species_cross_reference["mel1b"] = ["CP018508", "CP018509"]
    species_cross_reference["mel2"] = ["012441", "012442"]
    species_cross_reference["mel3"] = ["007760", "007761"]
    species_cross_reference["ovis"] = ["009504", "009505"]
    species_cross_reference["neo"] = ["KN046827"]
    species_cross_reference["suis1"] = ["017250", "017251"]
    species_cross_reference["suis2"] = ["NC_010169", "NC_010167"]
    species_cross_reference["suis3"] = ["007719", "007718"]
    species_cross_reference["suis4"] = ["B-REF-BS4-40"]
    species_cross_reference["te_09-0932"] = ["CP021201"]
    species_cross_reference["te_89-0490"] = ["CP021199"]
    species_cross_reference["te_92-0972"] = ["CP021060"]
    species_cross_reference["te_98-0554"] = ["CP021246"]
    species_cross_reference["te_atcc35865"] = ["NC_018108"]
    species_cross_reference["te_mce9"] = ["NC_014914"]
    species_cross_reference["gua"] = ["Gua1_1407_2018"]
    vcf_list = glob.glob('*vcf')
    for each_vcf in vcf_list:
        print(each_vcf)
        vcf_reader = vcf.Reader(open(each_vcf, 'r'))
        print("single_vcf %s" % each_vcf)
        for record in vcf_reader:
            header = record.CHROM
            for key, vlist in species_cross_reference.items():
                for li in vlist:
                    if li in header:
                        return (key)


def species_selection_step1(arg_options):
    all_parameters = Get_Specie_Parameters()

    if arg_options['species']:
        species_selection = arg_options['species']
        print("Sample will be ran as:  {}" .format(species_selection))
        parameters, genotype_codes = all_parameters.choose(species_selection)
    else:
        best_ref_found = best_reference([arg_options['R1'], arg_options['R2']])
        arg_options['species'] = best_ref_found
        print("Sample will be ran as {}" .format(best_ref_found))
        parameters, genotype_codes = all_parameters.choose(best_ref_found)

    if parameters['species'] is None:
        print("\n#### ERROR #####\nNo specie parameters found for: \n\t{} \n\t{}\n\n" .format(arg_options['R1'], arg_options['R2']))
        arg_options['reference_sequence_name'] = best_ref_found
        arg_options.update(parameters)
        return arg_options
    elif parameters:
        #shutil.copy2(parameters["reference"], arg_options['root_dir'])
        arg_options.update(parameters)
        return arg_options
    else:
        print("### See species_selection_step1 function")
        arg_options['species'] = None
        arg_options['reference_sequence_name'] = best_ref_found
        arg_options.update(parameters)
        return arg_options


def best_reference(fastq_list):

    '''Use oligos to determine species.  Most often if the absents of a single oligo from a set specific for either brucella or bovis will confer species type.  Some species will the absents of more than one oligo.  Oligo findings are translated to binary patterns.'''
    print("\nFinding the best reference\n")
    write_out = open("best_reference.txt", 'w')
    '''get the species'''
    oligo_dictionary = {}
    oligo_dictionary["01_ab1"] = "AATTGTCGGATAGCCTGGCGATAACGACGC"
    oligo_dictionary["02_ab3"] = "CACACGCGGGCCGGAACTGCCGCAAATGAC"
    oligo_dictionary["03_ab5"] = "GCTGAAGCGGCAGACCGGCAGAACGAATAT"
    oligo_dictionary["04_mel"] = "TGTCGCGCGTCAAGCGGCGTGAAATCTCTG"
    oligo_dictionary["05_suis1"] = "TGCGTTGCCGTGAAGCTTAATTCGGCTGAT"
    oligo_dictionary["06_suis2"] = "GGCAATCATGCGCAGGGCTTTGCATTCGTC"
    oligo_dictionary["07_suis3"] = "CAAGGCAGATGCACATAATCCGGCGACCCG"
    oligo_dictionary["08_ceti1"] = "GTGAATATAGGGTGAATTGATCTTCAGCCG"
    oligo_dictionary["09_ceti2"] = "TTACAAGCAGGCCTATGAGCGCGGCGTGAA"
    oligo_dictionary["10_canis4"] = "CTGCTACATAAAGCACCCGGCGACCGAGTT"
    oligo_dictionary["11_canis"] = "ATCGTTTTGCGGCATATCGCTGACCACAGC"
    oligo_dictionary["12_ovis"] = "CACTCAATCTTCTCTACGGGCGTGGTATCC"
    oligo_dictionary["13_ether2"] = "CGAAATCGTGGTGAAGGACGGGACCGAACC"
    oligo_dictionary["14_63B1"] = "CCTGTTTAAAAGAATCGTCGGAACCGCTCT"
    oligo_dictionary["15_16M0"] = "TCCCGCCGCCATGCCGCCGAAAGTCGCCGT"
    oligo_dictionary["16_mel1b"] = "TCTGTCCAAACCCCGTGACCGAACAATAGA" #added 2018-01-30
    oligo_dictionary["17_tb157"] = "CTCTTCGTATACCGTTCCGTCGTCACCATGGTCCT"
    oligo_dictionary["18_tb7"] = "TCACGCAGCCAACGATATTCGTGTACCGCGACGGT"
    oligo_dictionary["19_tbbov"] = "CTGGGCGACCCGGCCGACCTGCACACCGCGCATCA"
    oligo_dictionary["20_tb5"] = "CCGTGGTGGCGTATCGGGCCCCTGGATCGCGCCCT"
    oligo_dictionary["21_tb2"] = "ATGTCTGCGTAAAGAAGTTCCATGTCCGGGAAGTA"
    oligo_dictionary["22_tb3"] = "GAAGACCTTGATGCCGATCTGGGTGTCGATCTTGA"
    oligo_dictionary["23_tb4"] = "CGGTGTTGAAGGGTCCCCCGTTCCAGAAGCCGGTG"
    oligo_dictionary["24_tb6"] = "ACGGTGATTCGGGTGGTCGACACCGATGGTTCAGA"
    oligo_dictionary["25_para"] = "CCTTTCTTGAAGGGTGTTCG"
    oligo_dictionary["26_para_sheep"] = "CGTGGTGGCGACGGCGGCGGGCCTGTCTAT"
    oligo_dictionary["27_para_cattle"] = "TCTCCTCGGTCGGTGATTCGGGGGCGCGGT"

    brucella_identifications = {}
    brucella_identifications["1111111111111111"] = "odd" #Unexpected findings
    brucella_identifications["0111111111111111"] = "ab1" #Brucella abortus bv 1, 2 or 4
    brucella_identifications["1011111111111111"] = "ab3" #Brucella abortus bv 3
    brucella_identifications["1101111111111111"] = "ab1" #Brucella abortus bv 5, 6 or 9
    brucella_identifications["1110111111111101"] = "mel1"
    brucella_identifications["0000010101101101"] = "mel1"
    brucella_identifications["1110111111111100"] = "mel1b" #added 2018-01-30
    brucella_identifications["0000010101101100"] = "mel1b" #added 2018-01-30
    brucella_identifications["1110111111111011"] = "mel2"
    brucella_identifications["0000010101101001"] = "mel2"
    brucella_identifications["0100010101101001"] = "mel2"
    brucella_identifications["1110011111101011"] = "mel2"
    brucella_identifications["1110111111110111"] = "mel3"
    brucella_identifications["1110011111100111"] = "mel3"
    brucella_identifications["1111011111111111"] = "suis1"
    brucella_identifications["1111101111111111"] = "suis2"
    brucella_identifications["1111110111111101"] = "suis3"
    brucella_identifications["1111111011111111"] = "ceti1"
    brucella_identifications["1111111001111111"] = "ceti1"
    brucella_identifications["1111111101111111"] = "ceti2"
    brucella_identifications["1111111110111101"] = "suis4"
    brucella_identifications["1111111110011101"] = "canis"
    brucella_identifications["1111111111101111"] = "ovis"

    bovis_identifications = {}
    bovis_identifications["11101111"] = "h37" #tb1
    bovis_identifications["11101101"] = "h37" #tb1
    bovis_identifications["01100111"] = "h37" #tb2
    bovis_identifications["01101011"] = "h37" #tb3
    bovis_identifications["11101011"] = "h37" #tb3
    bovis_identifications["01101111"] = "h37" #tb4a
    bovis_identifications["01101101"] = "h37" #tb4b
    bovis_identifications["11101101"] = "h37" #tb4b
    bovis_identifications["01101111"] = "h37" #tb4b
    bovis_identifications["11111111"] = "h37" #tb5
    bovis_identifications["11001111"] = "h37" #tb6
    bovis_identifications["10101110"] = "h37" #tb7
    bovis_identifications["11001110"] = "af" #bovis
    bovis_identifications["11011110"] = "af" #bovis
    bovis_identifications["11001100"] = "af"  #bovis
    
    para_identifications = {}
    para_identifications["110"] = "para-CP033688"
    para_identifications["101"] = "para-NC_002944"
    
    count_summary = {}

    with futures.ProcessPoolExecutor() as pool:
        for v, count in pool.map(finding_best_ref, oligo_dictionary.values(), itertools_repeat(fastq_list)):
            for k, value in oligo_dictionary.items():
                if v == value:
                    count_summary.update({k: count})
                    count_summary = OrderedDict(sorted(count_summary.items()))

    count_list = []
    for v in count_summary.values():
        count_list.append(v)
    brucella_sum = sum(count_list[:16])
    bovis_sum = sum(count_list[16:24])
    para_sum = sum(count_list[24:])
    print("Best reference Brucella counts:", file=write_out)
    for i in count_list[:16]:
        print(i, end=',', file=write_out)
    print("\nBest reference TB counts:", file=write_out)
    for i in count_list[16:24]:
        print(i, end=',', file=write_out)

    print("\nBest reference Para counts:", file=write_out)
    for i in count_list[24:]:
        print(i, end=',', file=write_out)

    #Binary dictionary
    binary_dictionary = {}
    for k, v in count_summary.items():
        if v > 1:
            binary_dictionary.update({k: 1})
        else:
            binary_dictionary.update({k: 0})
    binary_dictionary = OrderedDict(sorted(binary_dictionary.items()))

    binary_list = []
    for v in binary_dictionary.values():
        binary_list.append(v)
    brucella_binary = binary_list[:16]
    brucella_string = ''.join(str(e) for e in brucella_binary)
    bovis_binary = binary_list[16:24]
    bovis_string = ''.join(str(e) for e in bovis_binary)
    para_binary = binary_list[24:]
    para_string = ''.join(str(e) for e in para_binary)

    if brucella_sum > 3:
        if brucella_string in brucella_identifications:
            print("Brucella group, species %s" % brucella_identifications[brucella_string])
            print("\n\nBrucella group, species %s" % brucella_identifications[brucella_string], file=write_out)
            return(brucella_identifications[brucella_string]) # return to set parameters
        else:
            print("Brucella group, but no match")
            print("\n\nBrucella group, but no match", file=write_out)
            return("Brucella group, but no match")
    elif bovis_sum > 3:
        if bovis_string in bovis_identifications:
            print("TB group, species %s" % bovis_identifications[bovis_string])
            print("\n\nTB group, species %s" % bovis_identifications[bovis_string], file=write_out)
            return(bovis_identifications[bovis_string]) # return to set parameters
        else:
            print("TB group, but no match")
            print("\n\nTB group, but no match", file=write_out)
            return("TB group, but no match")
    elif para_sum >= 1:
        if para_string in para_identifications:
            print("Para group, species %s" % para_identifications[para_string])
            print("\n\nPara group, species %s" % para_identifications[para_string], file=write_out)
            return(para_identifications[para_string]) # return to set parameters
        else:
            print("M. paratuberculosis group, but no match")
            print("\n\nNo match", file=write_out)
            return("M. paratuberculosis group, but no match")
    else:
        print("Unable to find a best reference species or group")
        return("Unable to find a best reference species or group")

    write_out.close()


def finding_best_ref(v, fastq_list):
    count = 0
    for fastq in fastq_list:
        with gzip.open(fastq, 'rt') as in_handle:
            # all 3, title and seq and qual, were needed
            for title, seq, qual in FastqGeneralIterator(in_handle):
                count += seq.count(v)
    return(v, count)


def align_reads(arg_options):
    paired = arg_options['paired']
    working_directory = os.getcwd()
    print("Working on: {}" .format(working_directory))
    ts = time.time()
    st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
    R1 = arg_options['R1']
    if paired:
        R2 = arg_options['R2']
    if arg_options["species"] is None:
        stat_summary = {}
        R1size = sizeof_fmt(os.path.getsize(R1))
        if paired:
            R2size = sizeof_fmt(os.path.getsize(R2))
            stat_summary["R2size"] = R2size
        stat_summary["time_stamp"] = st
        stat_summary["sample_name"] = arg_options["sample_name"]
        stat_summary["species"] = "NOT_FOUND"
        stat_summary["reference_sequence_name"] = "N/A"
        stat_summary["R1size"] = R1size
        stat_summary["allbam_mapped_reads"] = "CHECK SAMPLE *****************************************"
        stat_summary.update(arg_options['read_quality_stats'])
        return(stat_summary)
    else:
        startTime = datetime.now()
        print("\n\n*** START ***\n")
        print("Start time: %s" % startTime)
        sample_name = arg_options["sample_name"]
        print("species: %s" % arg_options["species"])
        if arg_options["species"] in ["ab1", "ab3", "suis1", "suis2", "suis3", "suis4", "mel1", "mel1b", "mel2", "mel3", "canis", "ceti1", "ceti2"]:
            print("Brucella")
            mlst(arg_options)
        elif arg_options["species"] in ["h37", "af"]: #removed bovis
            print("TB")
            spoligo(arg_options)
        os.chdir(working_directory)
        shutil.copy(arg_options["reference"], working_directory)
        sample_reference = glob.glob(working_directory + '/*fasta')

        print("reference: %s" % sample_reference)
        ref = re.sub('\.fasta', '', os.path.basename(sample_reference[0]))
        if len(sample_reference) != 1:
            print("### ERROR reference not available or too many")
            sys.exit(0)
        sample_reference = sample_reference[0]

        print("--")
        print("sample name: %s" % sample_name)
        print("sample reference: %s" % sample_reference)
        print("Read 1: %s" % R1)
        if paired:
            print("Read 2: %s\n" % R2)
        print("working_directory: %s" % working_directory)
        print("--")

        loc_sam = working_directory + "/" + sample_name
        os.system("samtools faidx {}" .format(sample_reference))
        os.system("picard CreateSequenceDictionary REFERENCE={} OUTPUT={}" .format(sample_reference, working_directory + "/" + ref + ".dict"))
        os.system("bwa index {}" .format(sample_reference))
        samfile = loc_sam + ".sam"
        bamfile = loc_sam + ".bam"
        abyss_out = loc_sam + "-unmapped_contigs.fasta"
        unmapsam = loc_sam + "-unmapped.sam"
        metrics = loc_sam + "-metrics.txt"
        unmapped_read1 = loc_sam + "-unmapped_R1.fastq"
        unmapped_read1gz = loc_sam + "-unmapped_R1.fastq.gz"
        if paired:
            unmapped_read2 = loc_sam + "-unmapped_R2.fastq"
            unmapped_read2gz = loc_sam + "-unmapped_R2.fastq.gz"
        sortedbam = loc_sam + "-sorted.bam"
        nodupbam = loc_sam + "-nodup.bam"
        unfiltered_hapall = loc_sam + "-unfiltered_hapall.vcf"
        mapq_fix = loc_sam + "-mapq_fix_hapall.vcf"
        hapall = loc_sam + "-hapall.vcf"
        zero_coverage_vcf = loc_sam + "_zc.vcf"

        #########################################################
        print("\n@@@ BWA mem: {}".format(sample_name))
        if paired:
            os.system(r'bwa mem -M -R "@RG\tID:%s\tSM:%s\tPL:ILLUMINA\tPI:250" -t 16 %s %s %s > %s' % (sample_name, sample_name, sample_reference, R1, R2, samfile))
        else:
            os.system(r'bwa mem -M -R "@RG\tID:%s\tSM:%s\tPL:ILLUMINA\tPI:250" -t 16 %s %s > %s' % (sample_name, sample_name, sample_reference, R1, samfile))
        os.system("samtools view -Sb {} -o {}" .format(samfile, bamfile))
        os.system("samtools sort {} -o {}" .format(bamfile, sortedbam))
        os.system("samtools index {}" .format(sortedbam))

        print("\n@@@ Remove Duplicate Reads: {}" .format(sample_name))
        os.system("picard MarkDuplicates INPUT={} OUTPUT={} METRICS_FILE={} ASSUME_SORTED=true REMOVE_DUPLICATES=true" .format(sortedbam, nodupbam, metrics))
        os.system("samtools index {}" .format(nodupbam))

        print("\n@@@ Calling SNPs with FreeBayes Parallel: {}" .format(sample_name))
        #os.system("freebayes -f {} {} > {}" .format(sample_reference, nodupbam, unfiltered_hapall))
        chrom_ranges = open("chrom_ranges.txt", 'w')
        for record in SeqIO.parse(sample_reference, "fasta"):
            chrom = record.id
            total_len = len(record.seq)
            min_number = 0
            step = 100000
            if step < total_len:
                for chunk in range(min_number, total_len, step)[1:]:
                    print("{}:{}-{}".format(chrom, min_number, chunk), file=chrom_ranges)
                    min_number = chunk
            print("{}:{}-{}".format(chrom, min_number, total_len), file=chrom_ranges)

            # for pos in list(range(last_number, total_len, 100000))[1:]:
            #     print("{}:{}-{}" .format(chrom, last_number, pos), file=chrom_ranges)
            #     last_number = pos
            # print("{}:{}-{}" .format(chrom, pos, total_len), file=chrom_ranges)
        chrom_ranges.close()
        os.system(r'freebayes-parallel chrom_ranges.txt 8 -E -1 -e 1 -u --strict-vcf -f %s %s > %s' % (sample_reference, nodupbam, unfiltered_hapall))
        # "fix" MQ notation in VCF to match GATK output
        write_fix = open(mapq_fix, 'w+')
        with open(unfiltered_hapall, 'r') as unfiltered:
            for line in unfiltered:
                line = line.strip()
                new_line = re.sub(r';MQM=', r';MQ=', line)
                new_line = re.sub(r'ID=MQM,', r'ID=MQ,', new_line)
                print(new_line, file=write_fix)
            write_fix.close()
        # remove clearly poor positions
        os.system(r'vcffilter -f "QUAL > 20" %s > %s' % (mapq_fix, hapall))

        print("\n@@@ Assemble Unmapped Reads: {}" .format(sample_name))
        os.system("samtools view -h -f4 -T {} {} -o {}".format(sample_reference, nodupbam, unmapsam))
        if paired:
            os.system("picard SamToFastq INPUT={} FASTQ={} SECOND_END_FASTQ={}".format(unmapsam, unmapped_read1, unmapped_read2))
        else:
            os.system("picard SamToFastq INPUT={} FASTQ={}".format(unmapsam, unmapped_read1))

        abyss_contig_count = 0
        try:
            if paired:
                os.system("ABYSS --out {} --coverage 5 --kmer 64 {} {}" .format(abyss_out, unmapped_read1, unmapped_read2))
                with open(abyss_out) as f:
                    for line in f:
                        abyss_contig_count += line.count(">")
            else:
                os.system("ABYSS --out {} --coverage 5 --kmer 64 {}" .format(abyss_out, unmapped_read1))
                with open(abyss_out) as f:
                    for line in f:
                        abyss_contig_count += line.count(">")
        except FileNotFoundError as ex:
            abyss_contig_count = 0
            # debug_log(ex, inspect.getframeinfo(inspect.currentframe()), sample_name, "Zero Abyss contigs")

        # Full bam stats
        stat_out = open("stat_align.txt", 'w')
        stat_out.write(os.popen("samtools idxstats {} " .format(sortedbam)).read())
        stat_out.close()
        with open("stat_align.txt") as f:
            first_line = f.readline()
            first_line = first_line.rstrip()
            first_line = re.split(':|\t', first_line)
            reference_sequence_name = str(first_line[0])
            allbam_mapped_reads = int(first_line[2])
        # Duplicate bam stats
        duplicate_stat_file = "duplicate_stat_align.txt"
        duplicate_stat_out = open(duplicate_stat_file, 'w')
        #os.system("samtools idxstats {} > {}" .format(sortedbam, stat_out)) Doesn't work when needing to std out.
        duplicate_stat_out.write(os.popen("samtools idxstats {} " .format(nodupbam)).read())
        duplicate_stat_out.close()
        with open(duplicate_stat_file) as f:
            for line_num, line in enumerate(f):
                if line_num == 1:
                    dup_line_two = line
                    dup_line_two = dup_line_two.split()
                    unmapped_reads = int(dup_line_two[3])
        try:
            print(f"{unmapped_reads} unmapped reads")
        except:
            unmapped_reads = "none_found"
            logging.debug("Zero unmapped reads")
            print(f"no unmapped reads")
        allbam_mapped_reads = "{:,}".format(allbam_mapped_reads)

        try:
            zero_coverage_vcf, good_snp_count, ave_coverage, genome_coverage = add_zero_coverage(sample_name, sample_reference, nodupbam, hapall, zero_coverage_vcf)
        except FileNotFoundError as ex:
            # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), sample_name, "ALIGNMENT ERROR, NO COVERAGE FILE")
            text = "ALIGNMENT ERROR, NO COVERAGE FILE " + sample_name
            msg = MIMEMultipart()
            msg['From'] = "tod.p.stuber@aphis.usda.gov"
            msg['To'] = "tod.p.stuber@aphis.usda.gov"
            msg['Date'] = formatdate(localtime=True)
            msg['Subject'] = "### No coverage file"
            msg.attach(MIMEText(text))
            smtp = smtplib.SMTP('10.10.8.12')
            smtp.send_message(msg)
            smtp.quit()

        if arg_options["gbk_file"] or arg_options['no_annotation']:
            annotated_vcf = loc_sam + "-annotated.vcf"
            gbk_file = arg_options['gbk_file']

            print("Putting gbk into indexed dataframe...")
            annotation_dict = {}
            for gbk in gbk_file:
                print("gbk file: %s" % gbk)
                write_out = open('temp.csv', 'w+')
                gbk_dict = SeqIO.to_dict(SeqIO.parse(gbk, "genbank"))
                gbk_chrome = list(gbk_dict.keys())[0]
                for key, value in gbk_dict.items():
                    for feature in value.features:
                        if "CDS" in feature.type or "rRNA" in feature.type:
                            myproduct = None
                            mylocus = None
                            mygene = None
                            try:
                                myproduct = feature.qualifiers['product'][0]
                            except KeyError:
                                pass
                            try:
                                mylocus = feature.qualifiers['locus_tag'][0]
                            except KeyError:
                                pass
                            try:
                                mygene = feature.qualifiers['gene'][0]
                            except KeyError:
                                pass
                            print(key, int(feature.location.start), int(feature.location.end), mylocus, myproduct, mygene, sep='\t', file=write_out)
                write_out.close()
                df = pd.read_csv('temp.csv', sep='\t', names=["chrom", "start", "stop", "locus", "product", "gene"])
                df = df.sort_values(['start', 'gene'], ascending=[True, False])
                df = df.drop_duplicates('start')
                pro = df.reset_index(drop=True)
                pro.index = pd.IntervalIndex.from_arrays(pro['start'], pro['stop'], closed='both')
                annotation_dict[gbk_chrome] = pro

            header_out = open('v_header.csv', 'w+')
            with open(zero_coverage_vcf) as fff:
                for line in fff:
                    if re.search('^#', line):
                        print(line.strip(), file=header_out)
            header_out.close()

            vcf_df = pd.read_csv(zero_coverage_vcf, sep='\t', header=None, names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"], comment='#')
            vcf_df['ABS_VALUE'] = vcf_df['CHROM'].map(str) + '-' + vcf_df['POS'].map(str)
            vcf_df = vcf_df.set_index('ABS_VALUE')

            annotate_condense_dict = {}
            for gbk_chrome, pro in annotation_dict.items():
                matching_chrom_df = vcf_df[vcf_df['CHROM'] == gbk_chrome]
                for index, row in matching_chrom_df.iterrows():
                    pos = row.POS
                    try:
                        aaa = pro.iloc[pro.index.get_loc(int(pos))][['chrom', 'locus', 'product', 'gene']]
                        try:
                            chrom, name, locus, tag = aaa.values[0]
                            annotate_condense_dict[str(chrom) + "-" + str(pos)] = "{}, {}, {}".format(name, locus, tag)
                        except ValueError:
                            # if only one annotation entire chromosome (such with flu) then having [0] fails
                            chrom, name, locus, tag = aaa.values
                            annotate_condense_dict[str(chrom) + "-" + str(pos)] = "{}, {}, {}".format(name, locus, tag)
                    except KeyError:
                        annotate_condense_dict[str(gbk_chrome) + "-" + str(pos)] = "No annotated product"

            annotate_df = pd.DataFrame.from_dict(annotate_condense_dict, orient='index', columns=["ID"])
            annotate_df.index.name = 'ABS_VALUE'
            vcf_df.drop(['ID'], axis=1, inplace=True)

            vcf_df = vcf_df.merge(annotate_df, how='left', left_index=True, right_index=True)
            vcf_df = vcf_df[["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"]]
            vcf_df.to_csv('v_annotated_body.csv', sep='\t', header=False, index=False)

            cat_files = ['v_header.csv', 'v_annotated_body.csv']
            with open(annotated_vcf, "wb") as outfile:
                for cf in cat_files:
                    with open(cf, "rb") as infile:
                        outfile.write(infile.read())
        try:
            os.remove('temp.csv')
            os.remove('v_header.csv')
            os.remove('v_annotated_body.csv')
        except FileNotFoundError as ex:
            # debug_log(ex, inspect.getframeinfo(inspect.currentframe()), sample_name)
            pass
        os.remove(samfile)
        os.remove(bamfile)
        os.remove(unmapsam)
        os.remove(sortedbam)
        os.remove(sortedbam + ".bai")
        os.remove(unfiltered_hapall)
        os.remove(sample_reference + ".amb")
        os.remove(sample_reference + ".ann")
        os.remove(sample_reference + ".bwt")
        os.remove(sample_reference + ".pac")
        os.remove(sample_reference + ".sa")
        os.remove(ref + ".dict")
        os.remove(duplicate_stat_file)
        os.remove("stat_align.txt")

        unmapped = working_directory + "/unmapped"
        os.makedirs(unmapped)

        newZip = zipfile.ZipFile(unmapped_read1gz, 'w')
        newZip.write(unmapped_read1, compress_type=zipfile.ZIP_DEFLATED)
        newZip.close()
        os.remove(unmapped_read1)
        if paired:
            newZip = zipfile.ZipFile(unmapped_read2gz, 'w')
            newZip.write(unmapped_read2, compress_type=zipfile.ZIP_DEFLATED)
            newZip.close()
            os.remove(unmapped_read2)

        try:
            shutil.move(unmapped_read1gz, unmapped)
            if paired:
                shutil.move(unmapped_read2gz, unmapped)
            shutil.move(abyss_out, unmapped)
        except FileNotFoundError as ex:
            # debug_log(ex, inspect.getframeinfo(inspect.currentframe()), sample_name)
            pass
        alignment = working_directory + "/alignment"
        os.makedirs(alignment)
        movefiles = glob.glob('*-*')
        for i in movefiles:
            shutil.move(i, alignment)
        try:
            shutil.move(sample_reference, alignment)
            shutil.move(sample_reference + ".fai", alignment)
        except shutil.Error:
            pass
        except FileNotFoundError:
            pass
        except FileExistsError:
            pass

        runtime = (datetime.now() - startTime)
        print("\n\nruntime: %s:  \n" % runtime)
        ave_coverage = "{:0.1f}".format(float(ave_coverage))
        print("average_coverage: %s" % ave_coverage)

        R1size = sizeof_fmt(os.path.getsize(R1))
        if paired:
            R2size = sizeof_fmt(os.path.getsize(R2))

        try:
            with open("mlst/mlst.txt") as f:
                first_line = f.readline()
                mlst_type = first_line.rstrip()
                first_line = first_line.split()
                mlst_type = first_line[1:]
                mlst_type = '-'.join(mlst_type)
        except FileNotFoundError:
            mlst_type = "N/A"

        try:
            with open("spoligo.txt") as f:
                first_line = f.readline()
                first_line = first_line.rstrip()
                first_line = first_line.split()
                octalcode = first_line[0]
                sbcode = first_line[1]
                hexcode = first_line[2]
                binarycode = first_line[3]
        except FileNotFoundError:
            octalcode = "N/A"
            sbcode = "N/A"
            hexcode = "N/A"
            binarycode = "N/A"
        #Capture program versions for step 1
        try:
            verison_out = open("version_capture.txt", 'w')
            print(os.popen('conda list bwa | grep -v "^#"; \
                conda list abyss | grep -v "^#"; \
                conda list picard | grep -v "^#"; \
                conda list samtools | grep -v "^#"; \
                conda list freebayes | grep -v "^#"; \
                conda list biopython | grep -v "^#"').read(), file=verison_out)
            print("Dependent source:  {}" .format(arg_options['script_dependents']), file=verison_out)
            verison_out.close()
        except:
            logging.debug("failed version capture")
            pass

        sequence_count = 0
        total_length = 0
        with gzip.open(R1, "rt") as handle:
            for r in SeqIO.parse(handle, "fastq"):
                total_length = total_length + len(r.seq)
                sequence_count = sequence_count + 1
        ave_read_length = total_length / sequence_count
        ave_read_length = "{:0.1f}".format(float(ave_read_length))

        ts = time.time()
        st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')

        stat_summary = {}
        stat_summary["time_stamp"] = st
        stat_summary["sample_name"] = sample_name
        stat_summary["species"] = arg_options["species"]
        stat_summary["reference_sequence_name"] = reference_sequence_name
        stat_summary["R1size"] = R1size
        if paired:
            stat_summary["R2size"] = R2size
        else:
            stat_summary["R2size"] = None
        stat_summary["allbam_mapped_reads"] = allbam_mapped_reads
        stat_summary["genome_coverage"] = genome_coverage
        stat_summary["ave_coverage"] = ave_coverage
        stat_summary["ave_read_length"] = ave_read_length
        stat_summary["unmapped_reads"] = unmapped_reads
        stat_summary["unmapped_assembled_contigs"] = abyss_contig_count
        stat_summary["good_snp_count"] = good_snp_count
        stat_summary["mlst_type"] = mlst_type
        stat_summary["octalcode"] = octalcode
        stat_summary["sbcode"] = sbcode
        stat_summary["hexadecimal_code"] = hexcode
        stat_summary["binarycode"] = binarycode
        ###
        # Create a sample stats file in the sample's script1 directory
        summary_file = loc_sam + "_" + st + '.xlsx'
        workbook = xlsxwriter.Workbook(summary_file)
        worksheet = workbook.add_worksheet()
        row = 0
        col = 0

        top_row_header = ["time_stamp", "sample_name", "species", "reference_sequence_name", "R1size", "R2size", "Q_ave_R1", "Q_ave_R2", "Q30_R1", "Q30_R2", "allbam_mapped_reads", "genome_coverage", "ave_coverage", "ave_read_length", "unmapped_reads", "unmapped_assembled_contigs", "good_snp_count", "mlst_type", "octalcode", "sbcode", "hexadecimal_code", "binarycode"]
        for header in top_row_header:
            worksheet.write(row, col, header)
            col += 1
            # worksheet.write(row, col, v)
        try:
            stat_summary.update(arg_options['read_quality_stats'])
        except KeyError as ex:
            # debug_log(ex, inspect.getframeinfo(inspect.currentframe()), sample_name)
            pass
        worksheet.write(1, 0, stat_summary.get('time_stamp', 'n/a'))
        worksheet.write(1, 1, stat_summary.get('sample_name', 'n/a'))
        worksheet.write(1, 2, stat_summary.get('species', 'n/a'))
        worksheet.write(1, 3, stat_summary.get('reference_sequence_name', 'n/a'))
        worksheet.write(1, 4, stat_summary.get('R1size', 'n/a'))
        worksheet.write(1, 5, stat_summary.get('R2size', 'n/a'))
        worksheet.write(1, 6, stat_summary.get('Q_ave_R1', 'n/a'))
        worksheet.write(1, 7, stat_summary.get('Q_ave_R2', 'n/a'))
        worksheet.write(1, 8, stat_summary.get('Q30_R1', 'n/a'))
        worksheet.write(1, 9, stat_summary.get('Q30_R2', 'n/a'))
        worksheet.write(1, 10, stat_summary.get('allbam_mapped_reads', 'n/a'))
        worksheet.write(1, 11, stat_summary.get('genome_coverage', 'n/a'))
        worksheet.write(1, 12, stat_summary.get('ave_coverage', 'n/a'))
        worksheet.write(1, 13, stat_summary.get('ave_read_length', 'n/a'))
        worksheet.write(1, 14, stat_summary.get('unmapped_reads', 'n/a'))
        worksheet.write(1, 15, stat_summary.get('unmapped_assembled_contigs', 'n/a'))
        worksheet.write(1, 16, stat_summary.get('good_snp_count', 'n/a'))
        worksheet.write(1, 17, stat_summary.get('mlst_type', 'n/a'))
        worksheet.write(1, 18, stat_summary.get('octalcode', 'n/a'))
        worksheet.write(1, 19, stat_summary.get('sbcode', 'n/a'))
        worksheet.write(1, 20, stat_summary.get('hexadecimal_code', 'n/a'))
        worksheet.write(1, 21, stat_summary.get('binarycode', 'n/a'))
        workbook.close()
        return(stat_summary)


def sizeof_fmt(num, suffix='B'):
    for unit in ['', 'K', 'M', 'G', 'T', 'P', 'E', 'Z']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)


def get_annotations(line, in_annotation_as_dict):
    #pos_found = False
    line = line.rstrip()
    if line.startswith("#"): # save headers to file
        return(line)
    elif not line.startswith("#"): # position rows
        #pos_found = False
        split_line = line.split('\t')
        chrom = split_line[0]
        position = split_line[1] # get position
    #print("Getting annotations")
        for each_key, each_value in in_annotation_as_dict.items():
            pos_found = False
            if chrom == each_key:
                for feature in each_value.features:
                    position = int(position)
                    #print(position)
                    if position in feature and ("CDS" in feature.type or "rRNA" in feature.type):
                        myproduct = "none list"
                        mylocus = "none list"
                        mygene = "none list"
                        for p in feature.qualifiers['product']:
                            myproduct = p
                        for l in feature.qualifiers['locus_tag']:
                            mylocus = l
                        if "gene" in feature.qualifiers:
                            gene = feature.qualifiers['gene']
                            for g in gene:
                                mygene = g
                        annotation_found = myproduct + ", gene: " + mygene + ", locus_tag: " + mylocus
                        pos_found = True
            if not pos_found:
                annotation_found = "No annotated product"
            #print(annotation_found)
            split_line[2] = annotation_found
            annotated_line = "\t".join(split_line)
            return(annotated_line)


def mlst(arg_options):

    if arg_options['debug_call']:
        with open("mlst-arg_options.json", 'w') as outfile:
                json.dump(arg_options, outfile)

    # with open("mlst-arg_options.json") as infile:
    #     arg_options = json.load(infile)

    sample_directory = str(os.getcwd())
    R1 = arg_options['R1']
    R2 = arg_options['R2']
    sample_name = arg_options['sample_name']

    #https://bmcmicrobiol.biomedcentral.com/articles/10.1186/1471-2180-7-34
    write_ref = open("ST1-MLST.fasta", 'w')
    print(">ST1-MLST", file=write_ref)
    print("CGTTTCCCCAAGGAAGTGGAGGTTGCAGGCGATACGATCGATGTTGGCTACGGCCCGATCAAGGTTCATGCCGTCCGCAACCCGGCCGAACTGCCGTGGAAGGAAGAAAACGTCGATATCGCCCTTGAATGCACCGGCATTTTCACCTCGCGCGACAAGGCAGCACTTCATCTTGAAGCTGGCGCCAAGCGCGTCATCGTCTCCGCTCCCGCAGACGGTGCCGATCTCACCGTCGTCTATGGTGTCAACAACGACAAGCTGACGAAGGACCATCTGGTCATCTCCAACGCTTCGTGTACCACCAACTGCCTTGCGCCGGTGGCTCAGGTTCTCAACGATACTATCGGTATCGAAAAGGGCTTCATGACCACGATCCACTCCTATACGGGCGACCAGCCGACGCTGGACACCATGCACAAGGATCTCTACCGCGCCCGCGCCGCTGCCCTTTCCATGATCCCGACCTCGACGGGTGCGGCCAAGGCCGTCGGTCTCGTTCTGCCGGAACTGAAAGGCAAGCTCGACGGCGTTGCCATTCGCGTCCCGACCCCAAATGTCTCGGTCGTTGATCTCACCTTCATCGCCAAGCGTGAAACCACCGTTGAAGAAGTCAACAATGCGATCCGCGAAGCCGCCAATGGCCGCCTCAAGGGCATTCTCGGCTATACCGATGAGAAGCTCGTCTCGCACGACTTCAACCACGATTCCCATTCCTCGGTCTTCCACACCGACCAGACCAAGGTTATGGACGGCACCATGGTGCGTATCCTGTCGTGGTACGACAATGAATGGGGCTTCTCCAGCCGCATGAGCGACACCGCCGTCGCTTTGGGCAAGCTGATCTGATAACGGCAACGCCTCTCCTTCACTGGCGAGGCGTTTTCATTTCTTGATAAGGACCGAGAGAAGAAACATGATGTTCCGCACCCTTGACGATGCCAATGTCCAATCCAAGCGCGTGCTGGTCCGTGTTGACCTCAACGTGCCGAAAATCCGCCGTTCTGCTCGCTGGTCTTAACACCCCGGGCGTCACCACCGTGATCGAGCCGGTCATGACGCGCGATCATACGGAAAAGATGCTGCAAGACTTTGGCGCAGACCTGACGGTTGAAACCGATAAGGATGGTGTGCGCCATATCCGTATTGTCGGCCAGGGCAAGCTTACCGGCCAGACCATCGACGTGCCGGGTGATCCCTCGTCAACGGCTTTTCCGCTGGTGGCCGCCCTTCTGGTCGAAGGTTCGGAGGTCACCATCCGCAATGTGCTGATGAACCCGACCCGCACCGGCCTGATCCTGACGTTGCAGGAAATGGGGGCGGATATCGAGATCATCGATCCACGCCTTGCCGGCGGCGAGGATGTCGCCGATCTGCGCGTCAAGGCCTCGAAGCTGAAAGGCGTTGTCGTTCCGCCGGAACGTGCGCCTTCGATGATCGATGAATATCCGGTTCTGGCCATTGCCGCGTCTTTTGCGGAAGGCGAAACCGTGATGGACGGTCTCGATGAACTGCGCGTCAAGGAATCGGATCGTCTGGCGGCCGTTGCGCGCGGCCTTGAAGCCAATGGTGTCGATTGTACCGAAGGCGAGATGTCGCTGACGGTTCGTGGCCGCCCCGGCGGCAAGGGGCTGGGCGGTGGCACGGTTGCAACCCACCTCGACCACCGCATCGCGATGAGTTTCCTCGTCATGGGCCTTGCATCGGAAAAGCCGGTTACGGTGGATGACAGCACCATGATCGCCACCTCTTTCCCGGAATTCATGGGCATGATGGCGGGGCTGGGGGCGAAGATTGCCGAAAGCGGTGCAGAATGAAATCGTTCGTCGTCGCCCCGTTCATTGTCGCCATTGACGGACCGGCCGCCTCGGGCAAGGGAACCCTTGCCCGGCGGATCGCGACACATTACGGGATGCCGCATCTCGATACGGGCCTGACCTATCGCGCGGTCGCCAAGAGCCGCGCTCTGTCATTCTGGCCGTGGCAGGCCCGGTGGACGGCGACGAGATCGACCTCACCAATTGCGACTGGGTCGTGCGTCCTAAAAAGATGATCGCTGATCTGGGCTTTGAAGACGTGACCGTCCTCAATGATTTCGAGGCGCAGGCCCTTGCCGTGGTTTCGCTGGAAGGCCACCATATGGAACAGATCGGCGGCAAACCGGAGGAGGCTGTTGCCACCCGCGTCGTGCTCGGCCCCGGCACGGGCCTTGGCGTGGCAGGTCTGTTTCGCACACGTCATGCATGGGTTCCGGTTCCCGGTGAAGGCGGTCATATCGATATCGGTCCACGCACCGAACGCGACTACCAGATTTTCCCGCATATCGAACGCATCGAAGGGCGTGTCACCGGCGAGCAAATTCTTAGCGGGCGGGGCCTGCGCAACCTCTATCTGGGCATCTGCGCCGCCGACAAGATCACGCCCACCCTTGAGACGCCAGTAGACATTACATCCGCCGGACTGGACGGCAGCAATCCACAAGCCGCAGAAACGCTTGACCTCTTCGCCACCTATCTGGGGCGGCTTGCGGGCGACCTTGCGCTCATTTTCATGGCGCATGGCGGCGTTTATCTTTCGGGTGGCATCCCGGTGCGCATCCTTTCCGCCCTCAAGGCCGGTTCGTTCCGCGCAACCTTCGAGGACAAGGCCCCGCACAAGGCCATCATGCGCGACATACCGGTCCGCGTTATCACATATCAACTGGCGGCCTTAACCGGGCTTTCCGCTTTCGCCCGCACCCCCTCGCGCTTTGAAGTTTCGACCGAGGGCCGCCGCTGGCGCATGCGCCGCTAGAGCATTTCCGAGCCAAAAGTGCGAAGCGGTTCCGTTTCCCAACGAGCCGACCGCGGCTGCGCTTGCCTATGGTCTCGACAAGAGCGAAGGCAAGACCATCGCTGTCTATGACCTTGGCGGCGGTACTTTCGACGTGTCGGTTCTGGAAATCGGCGACGGCGTTTTTGAAGTGAAGTCCACCAATGGCGACACGTTCCTTGGCGGTGAAGACTTCGATATTCGTCTGGTCGAATATCTGGTTGCCGAGTTCAAGAAGGAAAGTGGCATCGACCTGAAGAACGACAAGCTTGCCCTGCAGCGCCTCAAGGAAGCTGCCGAAAAGGCCAAGATCGAACTGTCGTCCTCGCAGCAGACCGAAATCAACCTGCCGTTCATCACGGCTGACCAGACTGGCCCGAAGCATCTGGCGATCAAGCTGTCGCGCGCCAAGTTTGAAAGCCTGGTCGATGATCTCGTGCAGCGCACGGTCGAGCCGTGCAAGGCGGCGCTCAAGGATGCCGGCCTCAAGGCTGGCGAAATTGACGAAGTGGTTCTGGTCGGCGGCATGACCCGCATGCCCAAGATTCAGGAAGTCGTGAAGGCCTTCTTCGGCAAGGAACCGCACAAGGGCGTGAACCCGGATGAAGTCGTGGCCATGGGCGCGGCGATCCAGGGCGGCGTTTTGCAGGGCGACGTCAAGGACGTGCTGCTGCTCGACGTGACCCCGCTTTCGCTCGGCATTGAAACGCTGGGCGGCGTGTTCACCCGCCTGATCGAACGCAACACCACTATCCCGACCAAGAAGTCGCAGACCTTCTCCACGGCTGAGGACAACCAGTCGGCCGTGACGATCCGCGTCTTCCAGGGCGAGCGTGAAATGGCAGCCGATAACAAGCTGCTTGGACAGTTCGACCTCGTTGGCATTCCGCCACGTCCCTGCCCGGAAAGCTTGCCGATTGCCAGGAGCGCGATCCGGCCAAGTCCGAAATCTTCATCGTCGAGGGCGATTCGGCAGGCGGTTCCGCCAAGAGCGGGCGCTCGCGCCAGAATCAGGCCATTCTGCCGCTGCGCGGCAAAATCCTGAACGTGGAACGCGTGCGTTTCGACCGGATGATTTCATCCGATCAGGTGGGCACCCTCATCACGGCGCTTGGCACCTCCATCGGCAAGGATGAAACGCACGGCTTCAACGCCGACAAGCTGCGTTATCACAAGATCATCATCATGACCGACGCCGACGTCGATGGCGCCCATATTCGTACGCTTCTGCTCACCTTCTTCTTCCGGCAGATGCCGGAACTGATCGAACGCGGGCATATCTATATCGCGCAGCCGCCGCTCTATAAGGTGACACGCGGCAAGTCTTCGCAATATATCAAGAACGAAGCCGCCTTTGAGGATTTCCTCATCGAAACCGGCCTTGAAGAAACGACACTGGAACTGGTGACTGGCGAAATGCGCGCCGGGCCGGATTTGCGCTCGGTGGTGGAGGATGCGCGCACGCTGCGTCAGCTTCTGCACGGCCTGCACACCCGCTATGACCGCAGCGTGGTGGAACAGGCGGCAATTGCCGGCCTGCTCAACCCCGATGCCTCAAGGGACAATGCAACGGCACAGCATTCCGCCGATACGGTTGCCAAGCGTCTCGACATGATTTCGGAAGAGACCGAGCGCGGCTGGAGCGGCCATGTGATGGAAGACGGCGGCTATCGCTTCGAGCGTATGGTGCGCGGTGTAAAGGATATCGCCATTCTCGACATGGCCCTGCTCGGCTCGGCCGATGCCCGCCAGGTCGACCGAGATCGAGATGTATTCCCGCCTGATCCATACGGTCGATCATATCGAAGGCCGCCTGCGTGACGGCATGGATGCGTTTGACGGCTTCCTCAGCCATGCATGGGCTGTGACGGTGACAGGCGCGCCGAAGCTGTGGGCAATGCGCTTTCTTGAGGAAAACGAACGCAGCCCGCGCGCATGGTATGGCGGCGCGATCGGCATGATGCATTTCAATGGCGATATGAATACAGGGCTGACGCTGCGCACCATCCGCATCAAGGATGGTGTGGCGGAAATCCGTGCAGGGGCGACGCTTCTGTTCGATTCCAACCCTGACGAGGAAGAAGCCGAGACCGAATTGAAGGCATCGGCCATGATTGCGGCTGTGCGGGACGCACAGAAGAGCAATCAGATCGCGGAAGAAAGTGTGGCGGCAAAGGTGGGTGAGGGGGTTTCGATCCTGCTGGTCGATCACGAGGATTCCTTCGTCCATACGCTTGCCAATTATTTCCGCCAGACGGGCGCCAAGGTTTCCACCGTGCGTTCACCGGTGGCAGAGGAGATATTCGACCGCGTCAATCCCGATCTGGTGGTGTTATCGCCGGGACCGGGCTCGCCGCAGGATTTCGATTGCAAGGCGACCATCGATAAGGCGCGCAAGCGCCAGCTTCCGATTTTTGGCGTCTGCCTCGGCCTTCAGGCACTGGCGGAAGCCTATGGCGGGGCGTTGCGCCAGCTTCGCGTTCCGGTGCATGGCAAGCCTTCACGCATCCGCGTATCAAAGCCGGAGCGCATTTTCTCCGGCTTGCCGGAGGAAGTGACGGTGGGGCGTTATCATTCGATCTTCGCCGATCCTGAACGCCTGCCGGATGATTTTCTCGTCACAGCCGAAACGGAAGACGGGATCATAGCCTGCGGTGGAGGTGGTGATGGTGCCGCCGGGCTCCAGCCTGCCTGCGGATGCGGGGCTTGTCGTGTTGCCCGGCACCAAATCCACGATTGCCGATCTGCTGGCGCTGCGTGAAAACGGCTGGGACCGCGAATTGGTCGCCCATGTGAAGCGGGGCGGGCATGTGCTTGGTATTTGCGGCGGGTTTCAAATGCTTGGACGGCGGATCAGTGACCCGGCGGGTATTGAAGGCAATGTGCGCGATATCGAGGGGCTGGGCCTTCTCGATATCGAGACGATGACGGAGCCGGAAAAAGTGGTTCGCAATGTTGAGGCGGTGTCGCTGCTGCATGATGAGCCGCTGGAGGGCTATGAAATCCACATCGGGCGCACCAGCGGGCCGGATATGGCGCGGCCATTTGCGCGTATCGGCGATCATGATGATGGGGCCGTCTCGCCCGATGGTCGTATCATGGGAACCTATCTCCACGGTATTTTCAGTGCGGATCGTTTCCGCCACCACTTTTTGCGCGCGCTGGGTGTGGAAGGCGGCCAGATGAATTATCGCGAGAGCGTCGAAGAGGCTCTGGGCGAACTGGCTGAAGGGCTGGAAGCCTCGCTGGATATTGATGGCCTGTTTGCGCTGGCATGATTGACGCCGCGAAGCCGAAAGCCTAGTGTCAAACCATGTGACAGGTTTTGCCGGAACGAATCCCCGGCAATACCAAAAGGGAATGCGACGGACGGACCCACGCCGGGCGTCTTTATCGCAGCCGACCCCGCGACTGTAGAGCGGAGAGGGAAGAGGCAAGCCGGGCAACCGGCAGCCACTGGAAATCAGATGCGATAATGCAACATCGCATTTTTGCCATCTTCTCGACAGATTATCTCCACACAATGGGGCATTTCGTGCCGCAATTACCCTCGATATGTCACCCCTGTCAGCGCGGCATGGGCGGTTTACTCCCGATGCTGCCCGCCCGATAAGGGACCGCGCAAAACGTAATTTGTGTAAGGAGAATGCCATGCGCACTCTTAAGTCTCTCGTAATCGTCTCGGCTGCGCTGCTGCCGTTCTCTGCGACCGCTTTTGCTGCCGACGCCATCCAGGAACAGCCTCCGGTTCCGGCTCCGGTTGAAGTAGCTCCCCAGTATAGCTGGGCTGGTGGCTATACCGGTCTTTACCTTGGCTATGGCTGGAACAAGGCCAAGACCAGCACCGTTGGCAGCATCAAGCCTGACGATTGGAAGGCTGGCGCCTTTGCTGGCTGGAACTTCCAGCAGGACCAGATCGTATACGGTGTTGAAGGTGATGCAGGTTATTCCTGGGCCAAGAAGTCCAAGGACGGCCTGGAAGTCAAGCAGGGCTTTGAAGGCTCGCTGCGTGCCCGCGTCGGCTACGACCTGAACCCGGTTATGCCGTACCTCACGGCTGGTATTGCCGGTTCGCAGATCAAGCTTAACAACGGCTTGGACGACGAAAGCAAGTTCCGCGTGGGTTGGACGGCTGGTGCCGGTCTCGAAGCCAAGCTGACGGACAACATCCTCGGCCGCGTTGAGTACCGTTACACCCAGTACGGCAACAAGAACTATGATCTGGCCGGTACGACTGTTCGCAACAAGCTGGACACGCAGGATATCCGCGTCGGCATCGGCTACAAGTTCTAATTATAGCATAATTGGACACGGAAAACCGGACAGGCAACTGTCCGGTTTTTTGTTGTCTGCAAAGGTCGAGAAAGCGCGGCAGAGCAACGGCGGCAGCCTGATTTTCAGGGGAAATGAAGTGGAGGCTTCTGTTGCCAGGTGCCTCCGAACCCCGCCTTAAGGGGCTAACCCTAAGGACTTTAGAGTGGGTTTCCCGCACCGCCATTAGGCAGCGAGAGCATAACCCTGAGCATTGTTGTCATTTGCAACTACTCTGTTGACCCGATAACGGTGGTATCATGCCGAGTAAAAGAGCGATCTTTACACCCTTGTCGATCCTGTTTCGCCCCCGCCACAACACAGCCTGATCGGCAAGCTGTGCTGTGGTGGAGGCGCCGGGTACCGCCCCCGGGTCCAATGGGTTTATTACACCGTCCGTTTATCACCATAGTCGGCTTGCGCCGACAGGACGTATATAGGCGTGGTTTTTACCGATTGGAAGGGGGCTTGTGCGTTTTCGCGCAAGACCGACAGAGGTGGTGCGGCCCTTCCGTTCATTTTCCATTGACAGCTTCCGCGTGCTGGTCAATCCTCACAATATATCGGGATCGGCCTTGAAGAGGCTTGGCGCAGCCGGGGCGGAAACCATGGCTGAAACGGGGACGATATGCCCCAATCGAAGGAGAGTGGATATATGAGTGAATATCTCGCGGATGTCCGTCGCTATGATGCTGCCGCCGATGAGGCCGTTGTCGAGAAAATCGTCAAGCATCTTGGCATTGCGCTTCGCAATCGCGATTCCTCGCTCGTTTCGGCAAGC", file=write_ref)
    write_ref.close()

    directory = str(os.getcwd())
    print(directory)
    sample_reference_mlst_location = directory + "/ST1-MLST.fasta"
    sample_name_mlst = sample_name + "-mlst"
    print("mlst reference: %s" % sample_reference_mlst_location)
    ref = re.sub('\.fasta', '', os.path.basename(sample_reference_mlst_location))
    print(ref)

    loc_sam_mlst = directory + "/" + sample_name_mlst
    print("\n--")
    print("sample name mlst: %s" % sample_name_mlst)
    print("sample reference mlst: %s" % sample_reference_mlst_location)
    print("ref, no ext: %s " % ref)
    print("Read 1: %s" % R1)
    print("Read 2: %s\n" % R2)
    print("directory: %s" % directory)
    print("loc_sam_mlst: %s" % loc_sam_mlst)
    print("--\n")

    os.system("samtools faidx {}" .format(sample_reference_mlst_location))
    os.system("picard CreateSequenceDictionary REFERENCE={} OUTPUT={}" .format(sample_reference_mlst_location, directory + "/" + ref + ".dict"))
    os.system("bwa index {}" .format(sample_reference_mlst_location))
    print("\n@@@ BWA mem")
    samfile_mlst = loc_sam_mlst + ".sam"
    os.system("bwa mem -M -t 16 {} {} {} > {}" .format(sample_reference_mlst_location, R1, R2, samfile_mlst))
    print("\nAdd read group and out all BAM")
    all_bam_mlst = loc_sam_mlst + "-all.bam"
    os.system("picard AddOrReplaceReadGroups INPUT={} OUTPUT={} RGLB=lib1 RGPU=unit1 RGSM={} RGPL=illumina" .format(samfile_mlst, all_bam_mlst, sample_name_mlst))

    print("\n@@@ Samtools mapped")
    mapbam = loc_sam_mlst + "-unmapped.bam"
    os.system("samtools view -h -F4 -b -T {} {} -o {}" .format(sample_reference_mlst_location, all_bam_mlst, mapbam))

    print("\n@@@ Sort BAM")
    sortedbam = loc_sam_mlst + "-sorted.bam"
    os.system("samtools sort {} -o {}" .format(mapbam, sortedbam))

    print("\n@@@ Index BAM")
    os.system("samtools index {}" .format(sortedbam))

    print("\n@@@ Calling SNPs with freebayes")
    unfiltered_vcf_mlst = directory + "/" + sample_name + "-unfiltered_mlst" + ".vcf"
    mapq_fix = loc_sam_mlst + "-mapq_fix_mlst.vcf"
    vcf_mlst = directory + "/" + sample_name + "_mlst" + ".vcf"

    os.system("freebayes -E -1 -f {} {} > {}" .format(sample_reference_mlst_location, sortedbam, unfiltered_vcf_mlst))
    # "fix" MQ notation in VCF to match GATK output
    write_fix = open(mapq_fix, 'w+')
    with open(unfiltered_vcf_mlst, 'r') as unfiltered:
        for line in unfiltered:
            line = line.strip()
            new_line = re.sub(r';MQM=', r';MQ=', line)
            print(new_line, file=write_fix)
        write_fix.close()
    # remove clearly poor positions
    os.system(r'vcffilter -f "QUAL > 20" %s > %s' % (mapq_fix, vcf_mlst))

    pos_call_dict = {}
    vcf_reader = vcf.Reader(open(vcf_mlst, 'r'))
    for record in vcf_reader:
        if record.ALT[0]:
            pos_call_dict.update({record.POS: str(record.ALT[0])})

    # Position 1629 was too close to the end of glk sequence.  Reads would not assemble properly to call possilbe SNP, therefore 100 bases of the gene were added.  Because of this all positions beyond this point are 100 more.  Same with position 1645 and 2693.

    target_pos_ref = {231: 'C', 297: 'T', 363: 'C', 398: 'C', 429: 'C', 523: 'G', 631: 'G', 730: 'G', 1247: 'G', 1296: 'C', 1342: 'G', 1381: 'A', 1648: 'C', 1685: 'C', 1741: 'C', 1754: 'G', 2165: 'A', 2224: 'T', 2227: 'C', 2297: 'G', 2300: 'A', 2344: 'A', 2352: 'G', 2403: 'C', 2530: 'G', 2557: 'G', 2578: 'G', 2629: 'A', 3045: 'A', 3054: 'G', 3118: 'G', 3295: 'C', 3328: 'C', 3388: 'A', 3966: 'C', 3969: 'G', 4167: 'G', 4271: 'C', 4296: 'G', 4893: 'C', 4996: 'G', 4998: 'T', 5058: 'G', 5248: 'A', 5672: 'G', 5737: 'C', 5928: 'A', 5963: 'G', 5984: 'C', 5987: 'C', 6025: 'G', 6045: 'G', 6498: 'G', 6499: 'C', 6572: 'A', 6627: 'T', 6715: 'C', 6735: 'T', 6745: 'G', 6785: 'T', 6810: 'C', 6828: 'C', 6845: 'C', 6864: 'G', 6875: 'C', 7382: 'G', 7432: 'G', 7464: 'G', 7594: 'G', 7660: 'T', 7756: 'A'}

    #pos_call_dict will replace target_pos_ref
    for key, value in pos_call_dict.items():
        if key in target_pos_ref.keys():
            target_pos_ref[key] = value
    ordered_combined_dict = OrderedDict(sorted(target_pos_ref.items()))
    combined_value_list = list(ordered_combined_dict.values())
    mlst_join = ''.join(combined_value_list)

    mlst_dictionary = {}
    mlst_dictionary["CTCCCGGGGCGACCCGATCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 01"
    mlst_dictionary["CTCCCGGGGCGACCCGAGCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 02"
    mlst_dictionary["CTCCCGTGGCGACCCGAGCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 03"
    mlst_dictionary["CTCCCGGGGCGACCCGAGCGAAGCGGGAAGGCCAAGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 04"
    mlst_dictionary["CTCCCGGGGCGACCCGATCGAAGCGGGAAGGCCACGGCGAGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 05"
    mlst_dictionary["TTCCTGGGGCAACCCGAGCGAGGCAGGGAGGCCGCGGCTCGTGAGCGGTCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 06"
    mlst_dictionary["CTTCCTGGCCGAGCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT"] = "MLST type 07"
    mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCCGTCTCGCGGTGCT"] = "MLST type 08"
    mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCTGTCTCGTGGTGCT"] = "MLST type 09"
    mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGGGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT"] = "MLST type 10"
    mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCTGTCTCGCGGTGCT"] = "MLST type 11"
    mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT"] = "MLST type 12"
    mlst_dictionary["CCCCCGGGCCGACTCGAGCGAAGCGAAGAGGCCACGGCGCGTGAGTGACCAGGCACCTATCCCACGGGGTA"] = "MLST type 13"
    mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAGTGCGTGAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 14"
    mlst_dictionary["CCCCCGGGCCGACCCGGGCGAAGCGGGGAGGCTACGGTGCGTGAGTGGCCAGGCACCTGTCCCGCGAGGTA"] = "MLST type 15"
    mlst_dictionary["CCCCCGGCCCGACCCGGGCGAAGCGGGGAGGCTACGGTGCGTGAGTGGCCAGGCACCTGTCCCGCGAGGTA"] = "MLST type 16"
    mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACCTGTCCCGCAGGGTA"] = "MLST type 17"
    mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACCTGTCCCGCAGGCTA"] = "MLST type 18"
    mlst_dictionary["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGACACGGCGCGTGAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 19"
    mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACATGTCCCGCAGGGTA"] = "MLST type 20"
    mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACATGCCCCGCAGGGTA"] = "MLST type 21"
    mlst_dictionary["CCCCCGGGCCGACCCGAGCGAGGCGGGGAGGCCACGGCGCGGGAGTGGCCAGACACCTGTCCTGCGGGGTA"] = "MLST type 22"
    mlst_dictionary["CCCCCGGGCTGACCCGAGCGAAACGGGGAAGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 23"
    mlst_dictionary["CCCCCGGGCTGACCCGAGCGGAACGGGGAAGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 23x"
    mlst_dictionary["CCCCCGGGCCGACCCGAGCAAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 24"
    mlst_dictionary["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 25"
    mlst_dictionary["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAAGCACCTGTTCCGCGGGGTA"] = "MLST type 26"
    mlst_dictionary["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGACCACGGCGCATAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 27"
    mlst_dictionary["CCCTCGGGCCGACCTGAGCGAAGCGGGGAGACCACGGCGCATAAGTGGCCAGGCTCCTGTCCCGCGGGGTA"] = "MLST type 28"

    remove_files = glob.glob('ST1-MLST*')
    for i in remove_files:
        os.remove(i)
    remove_files = glob.glob('*-mlst*')
    for i in remove_files:
        os.remove(i)
    remove_files = glob.glob('*_mlst.vcf.idx')
    for i in remove_files:
        os.remove(i)
    os.remove(unfiltered_vcf_mlst)

    write_out = open("mlst.txt", 'w')
    if mlst_join in mlst_dictionary:
        mlst_type = mlst_dictionary[mlst_join]
        print(mlst_type)
        print(mlst_type, file=write_out)
    else:
        print("NO MLST MATCH FOUND\n")
        print("NO MLST MATCH FOUND", file=write_out)
    write_out.close()

    os.makedirs("mlst")
    shutil.move(vcf_mlst, "mlst")
    shutil.move("mlst.txt", "mlst")
    os.chdir(sample_directory)


def finding_sp(v):
    total = 0
    total_finds = 0
    #if total < 6: # doesn't make a big different.  Might as well get full counts
    #total += sum(seq.count(x) for x in (v)) #v=list of for and rev spacer
    total_finds = [len(regex.findall("(" + spacer + "){s<=1}", seq_string)) for spacer in v]
    for number in total_finds:
        total += number
    return(total)


def binary_to_octal(binary):
    #binary_len = len(binary)
    i = 0
    ie = 1
    octal = ""
    while ie < 43:
        ie = i + 3
        print(binary[i:ie])
        region = binary[i:ie]
        region_len = len(region)
        i += 3
        if int(region[0]) == 1:
            if region_len < 2: # for the lone spacer 43.  When present needs to be 1 not 4.
                oct = 1
            else:
                oct = 4
        else:
            oct = 0
        try:
            if int(region[1]) == 1:
                oct += 2
            if int(region[2]) == 1:
                oct += 1
        except IndexError:
            pass
        octal = octal + str(oct)
    return(octal)


def binary_to_hex(binary):
    section1 = binary[0:7]
    section2 = binary[7:14]
    section3 = binary[14:21]
    section4 = binary[21:28]
    section5 = binary[28:36]
    section6 = binary[36:43]

    hex_section1 = hex(int(section1, 2))
    hex_section2 = hex(int(section2, 2))
    hex_section3 = hex(int(section3, 2))
    hex_section4 = hex(int(section4, 2))
    hex_section5 = hex(int(section5, 2))
    hex_section6 = hex(int(section6, 2))

    return(hex_section1.replace('0x', '').upper() + "-" + hex_section2.replace('0x', '').upper() + "-" + hex_section3.replace('0x', '').upper() + "-" + hex_section4.replace('0x', '').upper() + "-" + hex_section5.replace('0x', '').upper() + "-" + hex_section6.replace('0x', '').upper())


def spoligo(arg_options):
    sample_directory = str(os.getcwd())
    R1 = arg_options['R1']
    R2 = arg_options['R2']
    print("\nFinding spoligotype pattern...\n")
    '''spoligo spacers'''
    spoligo_dictionary = {}
    spoligo_dictionary["spacer01"] = ["TGATCCAGAGCCGGCGACCCTCTAT", "ATAGAGGGTCGCCGGCTCTGGATCA"]
    spoligo_dictionary["spacer02"] = ["CAAAAGCTGTCGCCCAAGCATGAGG", "CCTCATGCTTGGGCGACAGCTTTTG"]
    spoligo_dictionary["spacer03"] = ["CCGTGCTTCCAGTGATCGCCTTCTA", "TAGAAGGCGATCACTGGAAGCACGG"]
    spoligo_dictionary["spacer04"] = ["ACGTCATACGCCGACCAATCATCAG", "CTGATGATTGGTCGGCGTATGACGT"]
    spoligo_dictionary["spacer05"] = ["TTTTCTGACCACTTGTGCGGGATTA", "TAATCCCGCACAAGTGGTCAGAAAA"]
    spoligo_dictionary["spacer06"] = ["CGTCGTCATTTCCGGCTTCAATTTC", "GAAATTGAAGCCGGAAATGACGACG"]
    spoligo_dictionary["spacer07"] = ["GAGGAGAGCGAGTACTCGGGGCTGC", "GCAGCCCCGAGTACTCGCTCTCCTC"]
    spoligo_dictionary["spacer08"] = ["CGTGAAACCGCCCCCAGCCTCGCCG", "CGGCGAGGCTGGGGGCGGTTTCACG"]
    spoligo_dictionary["spacer09"] = ["ACTCGGAATCCCATGTGCTGACAGC", "GCTGTCAGCACATGGGATTCCGAGT"]
    spoligo_dictionary["spacer10"] = ["TCGACACCCGCTCTAGTTGACTTCC", "GGAAGTCAACTAGAGCGGGTGTCGA"]
    spoligo_dictionary["spacer11"] = ["GTGAGCAACGGCGGCGGCAACCTGG", "CCAGGTTGCCGCCGCCGTTGCTCAC"]
    spoligo_dictionary["spacer12"] = ["ATATCTGCTGCCCGCCCGGGGAGAT", "ATCTCCCCGGGCGGGCAGCAGATAT"]
    spoligo_dictionary["spacer13"] = ["GACCATCATTGCCATTCCCTCTCCC", "GGGAGAGGGAATGGCAATGATGGTC"]
    spoligo_dictionary["spacer14"] = ["GGTGTGATGCGGATGGTCGGCTCGG", "CCGAGCCGACCATCCGCATCACACC"]
    spoligo_dictionary["spacer15"] = ["CTTGAATAACGCGCAGTGAATTTCG", "CGAAATTCACTGCGCGTTATTCAAG"]
    spoligo_dictionary["spacer16"] = ["CGAGTTCCCGTCAGCGTCGTAAATC", "GATTTACGACGCTGACGGGAACTCG"]
    spoligo_dictionary["spacer17"] = ["GCGCCGGCCCGCGCGGATGACTCCG", "CGGAGTCATCCGCGCGGGCCGGCGC"]
    spoligo_dictionary["spacer18"] = ["CATGGACCCGGGCGAGCTGCAGATG", "CATCTGCAGCTCGCCCGGGTCCATG"]
    spoligo_dictionary["spacer19"] = ["TAACTGGCTTGGCGCTGATCCTGGT", "ACCAGGATCAGCGCCAAGCCAGTTA"]
    spoligo_dictionary["spacer20"] = ["TTGACCTCGCCAGGAGAGAAGATCA", "TGATCTTCTCTCCTGGCGAGGTCAA"]
    spoligo_dictionary["spacer21"] = ["TCGATGTCGATGTCCCAATCGTCGA", "TCGACGATTGGGACATCGACATCGA"]
    spoligo_dictionary["spacer22"] = ["ACCGCAGACGGCACGATTGAGACAA", "TTGTCTCAATCGTGCCGTCTGCGGT"]
    spoligo_dictionary["spacer23"] = ["AGCATCGCTGATGCGGTCCAGCTCG", "CGAGCTGGACCGCATCAGCGATGCT"]
    spoligo_dictionary["spacer24"] = ["CCGCCTGCTGGGTGAGACGTGCTCG", "CGAGCACGTCTCACCCAGCAGGCGG"]
    spoligo_dictionary["spacer25"] = ["GATCAGCGACCACCGCACCCTGTCA", "TGACAGGGTGCGGTGGTCGCTGATC"]
    spoligo_dictionary["spacer26"] = ["CTTCAGCACCACCATCATCCGGCGC", "GCGCCGGATGATGGTGGTGCTGAAG"]
    spoligo_dictionary["spacer27"] = ["GGATTCGTGATCTCTTCCCGCGGAT", "ATCCGCGGGAAGAGATCACGAATCC"]
    spoligo_dictionary["spacer28"] = ["TGCCCCGGCGTTTAGCGATCACAAC", "GTTGTGATCGCTAAACGCCGGGGCA"]
    spoligo_dictionary["spacer29"] = ["AAATACAGGCTCCACGACACGACCA", "TGGTCGTGTCGTGGAGCCTGTATTT"]
    spoligo_dictionary["spacer30"] = ["GGTTGCCCCGCGCCCTTTTCCAGCC", "GGCTGGAAAAGGGCGCGGGGCAACC"]
    spoligo_dictionary["spacer31"] = ["TCAGACAGGTTCGCGTCGATCAAGT", "ACTTGATCGACGCGAACCTGTCTGA"]
    spoligo_dictionary["spacer32"] = ["GACCAAATAGGTATCGGCGTGTTCA", "TGAACACGCCGATACCTATTTGGTC"]
    spoligo_dictionary["spacer33"] = ["GACATGACGGCGGTGCCGCACTTGA", "TCAAGTGCGGCACCGCCGTCATGTC"]
    spoligo_dictionary["spacer34"] = ["AAGTCACCTCGCCCACACCGTCGAA", "TTCGACGGTGTGGGCGAGGTGACTT"]
    spoligo_dictionary["spacer35"] = ["TCCGTACGCTCGAAACGCTTCCAAC", "GTTGGAAGCGTTTCGAGCGTACGGA"]
    spoligo_dictionary["spacer36"] = ["CGAAATCCAGCACCACATCCGCAGC", "GCTGCGGATGTGGTGCTGGATTTCG"]
    spoligo_dictionary["spacer37"] = ["CGCGAACTCGTCCACAGTCCCCCTT", "AAGGGGGACTGTGGACGAGTTCGCG"]
    spoligo_dictionary["spacer38"] = ["CGTGGATGGCGGATGCGTTGTGCGC", "GCGCACAACGCATCCGCCATCCACG"]
    spoligo_dictionary["spacer39"] = ["GACGATGGCCAGTAAATCGGCGTGG", "CCACGCCGATTTACTGGCCATCGTC"]
    spoligo_dictionary["spacer40"] = ["CGCCATCTGTGCCTCATACAGGTCC", "GGACCTGTATGAGGCACAGATGGCG"]
    spoligo_dictionary["spacer41"] = ["GGAGCTTTCCGGCTTCTATCAGGTA", "TACCTGATAGAAGCCGGAAAGCTCC"]
    spoligo_dictionary["spacer42"] = ["ATGGTGGGACATGGACGAGCGCGAC", "GTCGCGCTCGTCCATGTCCCACCAT"]
    spoligo_dictionary["spacer43"] = ["CGCAGAATCGCACCGGGTGCGGGAG", "CTCCCGCACCCGGTGCGATTCTGCG"]

    count_summary = {}

    global seq_string
    sequence_list = []
    try:
        for fastq in R1, R2:
            with gzip.open(fastq, "rt") as in_handle:
                # all 3, title and seq and qual, were needed
                for title, seq, qual in FastqGeneralIterator(in_handle):
                    sequence_list.append(seq)
    except TypeError:
        # TypeError if not paired
        pass

    if len(seq) > 99:
        #Three 10bp sequences dispersed across repeat region, forward and reverse
        capture_spacer_sequence = re.compile(".*TTTCCGTCCC.*|.*GGGACGGAAA.*|.*TCTCGGGGTT.*|.*AACCCCGAGA.*|.*TGGGTCTGAC.*|.*GTCAGACCCA.*")
        sequence_list = list(filter(capture_spacer_sequence.match, sequence_list))
        seq_string = "".join(sequence_list)
    else:
        #if < 100 then search all reads, not just those with repeat regions.
        seq_string = "".join(sequence_list)

    # for spacer_id, spacer_sequence in spoligo_dictionary.items():
    #     count = finding_sp(spacer_sequence)
    #     count_summary.update({spacer_id: count})
    # count_summary = OrderedDict(sorted(count_summary.items()))
    # print("count_summary {}" .format(count_summary))

    for spacer_id, spacer_sequence in spoligo_dictionary.items():
        count = delayed(finding_sp)(spacer_sequence)
        count_summary.update({spacer_id: count})
    pull = delayed(count_summary)
    count_summary = pull.compute()
    count_summary = OrderedDict(sorted(count_summary.items()))
    print("count_summary {}".format(count_summary))

    seq_string = ""

    spoligo_binary_dictionary = {}
    for k, v in count_summary.items():
        if v > 4:
            spoligo_binary_dictionary.update({k: 1})
        else:
            spoligo_binary_dictionary.update({k: 0})
    spoligo_binary_dictionary = OrderedDict(sorted(spoligo_binary_dictionary.items()))
    spoligo_binary_list = []
    for v in spoligo_binary_dictionary.values():
        spoligo_binary_list.append(v)
    bovis_string = ''.join(str(e) for e in spoligo_binary_list) #bovis_string correct
    hexadecimal = binary_to_hex(bovis_string)
    write_out = open("spoligo.txt", 'w')
    found = False
    with open(arg_options["spoligo_db"]) as f: # put into dictionary or list
        for line in f:
            line = line.rstrip()
            octalcode = line.split()[0] #no arg splits on whitespace
            sbcode = line.split()[1]
            binarycode = line.split()[2]
            if bovis_string == '0000000000000000000000000000000000000000000':
                found = True
                octalcode = "spoligo not found"
                sbcode = "spoligo not found"
                hexadecimal = "SB2277 ???"
                binarycode = "0000000000000000000000000000000000000000000"
                print("CHECK SAMPLE!  NO SPACERS FOUND.  LIKELY NOT TB COMPLEX.  ALTHOUGH SB2277 IS A ZERO STRING BINARY\n")
                print("CHECK SAMPLE!  NO SPACERS FOUND.  LIKELY NOT TB COMPLEX.  ALTHOUGH SB2277 IS A ZERO STRING BINARY", file=write_out)
                print("\nOne mismatch allowed spacers search against both R1 and R2 reads.\n", file=write_out)
                for k, v in count_summary.items():
                    print(k, v, file=write_out)
            elif bovis_string == binarycode:
                found = True
                print("Pattern found:")
                print("%s %s %s %s" % (octalcode, sbcode, hexadecimal, binarycode))
                print("%s %s %s %s" % (octalcode, sbcode, hexadecimal, binarycode), file=write_out)
                print("\One mismatch allowed spacer search against both R1 and R2 reads.\n", file=write_out)
                for k, v in count_summary.items():
                    print(k, v, file=write_out)

                print("bovis_string: %s" % bovis_string, file=write_out)
                print("binarycode  : %s" % binarycode, file=write_out)

        if not found:
            octal = binary_to_octal(bovis_string)
            sbcode = "N/A"
            print("%s %s %s %s" % (octal, sbcode, hexadecimal, bovis_string))
            print("%s %s %s %s" % (octal, sbcode, hexadecimal, bovis_string), file=write_out)
            print("SPOLIGO SB NUMBER NOT FOUND\n")
            print("\nSPOLIGO SB NUMBER NOT FOUND\n", file=write_out)
            print("\nOne mismatch allowed spacer search against both R1 and R2 reads.\n", file=write_out)
            for k, v in count_summary.items():
                print(k, v, file=write_out)

        write_out.close()
    os.chdir(sample_directory)


def add_zero_coverage(sample_name, sample_reference, nodupbam, hapall, zero_coverage_vcf):
    print("\n@@@ Depth of coverage using pysam: {}"  .format(sample_name))
    coverage_dict = {}
    coverage_list = pysam.depth(nodupbam, split_lines=True)
    for line in coverage_list:
        chrom, position, depth = line.split('\t')
        coverage_dict[chrom + "-" + position] = depth
    coverage_df = pd.DataFrame.from_dict(coverage_dict, orient='index', columns=["depth"])
    zero_dict = {}
    for record in SeqIO.parse(sample_reference, "fasta"):
        chrom = record.id
        total_len = len(record.seq)
        for pos in list(range(1, total_len + 1)):
            zero_dict[str(chrom) + "-" + str(pos)] = 0
    zero_df = pd.DataFrame.from_dict(zero_dict, orient='index', columns=["depth"])
    #df with depth_x and depth_y columns, depth_y index is NaN
    coverage_df = zero_df.merge(coverage_df, left_index=True, right_index=True, how='outer')
    #depth_x "0" column no longer needed
    coverage_df = coverage_df.drop(columns=['depth_x'])
    coverage_df = coverage_df.rename(columns={'depth_y': 'depth'})
    #covert the NaN to 0 coverage
    coverage_df = coverage_df.fillna(0)
    coverage_df['depth'] = coverage_df['depth'].apply(int)

    print("...coverage found")
    total_length = len(coverage_df)
    ave_coverage = coverage_df['depth'].mean()

    zero_df = coverage_df[coverage_df['depth'] == 0]
    total_zero_coverage = len(zero_df)
    print("Total zero coverage positions: {:,}" .format(total_zero_coverage))
    total_coverage = total_length - total_zero_coverage
    genome_coverage = "{:.2%}".format(total_coverage / total_length)

    vcf_df = pd.read_csv(hapall, sep='\t', header=None, names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"], comment='#')
    good_snp_count = len(vcf_df[(vcf_df['ALT'].str.len() == 1) & (vcf_df['REF'].str.len() == 1) & (vcf_df['QUAL'] > 150)])

    if total_zero_coverage > 0:
        header_out = open('v_header.csv', 'w+')
        with open(hapall) as fff:
            for line in fff:
                if re.search('^#', line):
                    print(line.strip(), file=header_out)
        header_out.close()
        vcf_df_snp = vcf_df[vcf_df['REF'].str.len() == 1]
        vcf_df_snp = vcf_df_snp[vcf_df_snp['ALT'].str.len() == 1]
        vcf_df_snp['ABS_VALUE'] = vcf_df_snp['CHROM'].map(str) + '-' + vcf_df_snp['POS'].map(str)
        vcf_df_snp = vcf_df_snp.set_index('ABS_VALUE')
        cat_df = pd.concat([vcf_df_snp, zero_df], axis=1, sort=False)
        cat_df = cat_df.drop(columns=['CHROM', 'POS', 'depth'])
        cat_df[['ID', 'ALT', 'QUAL', 'FILTER', 'INFO']] = cat_df[['ID', 'ALT', 'QUAL', 'FILTER', 'INFO']].fillna('.')
        cat_df['REF'] = cat_df['REF'].fillna('N')
        cat_df['FORMAT'] = cat_df['FORMAT'].fillna('GT')
        cat_df['Sample'] = cat_df['Sample'].fillna('./.')
        cat_df['temp'] = cat_df.index.str.split('-')
        cat_df[['CHROM', 'POS']] = pd.DataFrame(cat_df.temp.values.tolist(), index=cat_df.index)
        cat_df = cat_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Sample']]
        cat_df['POS'] = cat_df['POS'].astype(int)
        cat_df = cat_df.sort_values(['CHROM', 'POS'])
        cat_df.to_csv('v_annotated_body.csv', sep='\t', header=False, index=False)
        cat_files = ['v_header.csv', 'v_annotated_body.csv']
        with open(zero_coverage_vcf, "wb") as outfile:
            for cf in cat_files:
                with open(cf, "rb") as infile:
                    outfile.write(infile.read())
    else:
        shutil.copyfile(hapall, zero_coverage_vcf)
    return (zero_coverage_vcf, good_snp_count, ave_coverage, genome_coverage)


def send_email_step1(email_list, runtime, path_found, summary_file, st):
    text = "See attached:  "
    send_from = "tod.p.stuber@aphis.usda.gov"
    send_to = email_list
    msg = MIMEMultipart()
    msg['From'] = send_from
    msg['To'] = send_to
    msg['Date'] = formatdate(localtime=True)
    if not path_found:
        msg['Subject'] = "###CUMULATIVE STATS NOT UPDATED - Script1 stats summary"
    else:
        msg['Subject'] = "Script1 stats summary, runtime: {}" .format(runtime)
    msg.attach(MIMEText(text))

    part = MIMEBase('application', "octet-stream")
    part.set_payload(open(summary_file, "rb").read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition', 'attachment; filename="stat_summary_{}.xlsx"' .format(st))
    msg.attach(part)

    #context = ssl.SSLContext(ssl.PROTOCOL_SSLv3)
    #SSL connection only working on Python 3+
    smtp = smtplib.SMTP('10.10.8.12')

    smtp.send_message(msg)
    #smtp.sendmail(send_from, send_to, msg.as_string())
    smtp.quit()


def group_files(each_vcf, arg_options):

    mal = ""
    list_pass = []
    list_amb = []
    dict_amb = {}
    group_calls = []

    try:
        vcf_reader = vcf.Reader(open(each_vcf, 'r'))
        # PUT VCF NAME INTO LIST, capturing for htmlfile
        group_calls.append(each_vcf)
        # for each single vcf getting passing position
        for record in vcf_reader:
            try:
                # Freebayes VCFs place MQ values are placed into a list.  GATK as a float
                record.INFO['MQ'] = record.INFO['MQ'][0]
            except TypeError:
                pass
            except KeyError:
                pass
            chrom = record.CHROM
            position = record.POS
            absolute_positon = str(chrom) + "-" + str(position)
            # find quality SNPs and put absolute positions into list
            try:
                record_alt_length = len(record.ALT[0])
            except TypeError:
                record_alt_length = 0
            try:
                record_ref_length = len(record.REF)
            except TypeError:
                record_alt_length = 0
            try:
                if str(record.ALT[0]) != "None" and record_ref_length == 1 and record_alt_length == 1 and record.INFO['AC'][0] == 2 and record.QUAL > arg_options['qual_threshold'] and record.INFO['MQ'] > 45:
                    list_pass.append(absolute_positon)
                # capture ambigous defining SNPs in htmlfile
                elif str(record.ALT[0]) != "None" and record.INFO['AC'][0] == 1:
                    list_amb.append(absolute_positon)
            except ZeroDivisionError as ex:
                # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), each_vcf, absolute_positon)
                print("bad line in %s at %s" % (each_vcf, absolute_positon))

        for key in arg_options['inverted_position'].keys():
            if key not in list_pass:
                print("key %s not in list_pass" % key)
                directory = arg_options['inverted_position'][key]
                print("*** INVERTED POSITION FOUND *** PASSING POSITION FOUND: \t%s\t\t%s" % (each_vcf, directory))
                if not os.path.exists(directory):
                    try:
                        os.makedirs(directory)
                    except FileExistsError:
                        pass
                shutil.copy(each_vcf, directory)
                # ADD GROUP TO LIST
                group_calls.append(directory)

        #if passing:
        # if a passing position is in the defining SNPs
        defining_snps = arg_options['defining_snps']
        for passing_position in list_pass:
            # normal grouping
            if passing_position in defining_snps:
                directory = defining_snps[passing_position]
                print("PASSING POSITION FOUND: \t%s\t\t%s" % (each_vcf, directory))
                if not os.path.exists(directory):
                    try:
                        os.makedirs(directory)
                    except FileExistsError:
                        pass
                shutil.copy(each_vcf, directory)
                # ADD GROUP TO LIST
                group_calls.append(directory)
        # find mixed isolates if defining snp is ambigous
        for amb_position in list_amb:
            if amb_position in defining_snps:
                directory = defining_snps[amb_position]
                dict_amb.update({each_vcf + "\t" + directory: amb_position})
                # ADD AMBIGIOUS CALL TO LIST
                group_calls.append("*" + directory + "-mix")
                if not os.path.exists(directory):
                    try:
                        os.makedirs(directory)
                    except FileExistsError:
                        pass
                shutil.copy(each_vcf, directory)
        # if -a or -e (non elites already deleted from the analysis) copy all vcfs to All_VCFs
        if arg_options['all_vcf'] or arg_options['elite']:
            if not os.path.exists("All_VCFs"):
                os.makedirs("All_VCFs")
            shutil.move(each_vcf, "All_VCFs")
        else:
            try:
                os.remove(each_vcf)
            except FileNotFoundError:
                pass
        #print(dict_amb, group_calls, malformed)

    except ZeroDivisionError as ex:
        os.remove(each_vcf)
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), each_vcf)
        mal = "ZeroDivisionError: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")
    except ValueError as ex:
        os.remove(each_vcf)
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), each_vcf)
        mal = "ValueError: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")
    except UnboundLocalError as ex:
        os.remove(each_vcf)
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), each_vcf)
        mal = "UnboundLocalError: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")
    except TypeError as ex:
        os.remove(each_vcf)
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), each_vcf)
        mal = "TypeError: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")
    except SyntaxError as ex:
        os.remove(each_vcf)
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), each_vcf)
        mal = "SyntaxError: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")
    except KeyError as ex:
        os.remove(each_vcf)
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), each_vcf)
        mal = "KeyError: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), each_vcf)
    except StopIteration as ex:
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), each_vcf)
        mal = "StopIteration: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")
    except IndexError as ex:
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), each_vcf)
        mal = "IndexError: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")

    the_sample_name = group_calls[0:1]
    list_of_groups = sorted(group_calls[1:]) # order the groups
    for i in list_of_groups:
        the_sample_name.append(i) # a is group_calls
        group_calls = the_sample_name
    return dict_amb, group_calls, mal


def run_script2(arg_options):

    # IF AVX2 IS AVAILABE (CHECK WITH `cat /proc/cpuinfo | grep -i "avx"`). CREATE A LINK TO: `ln -s path_to_raxmlHPC-PTHREADS-AVX2 raxml.  Place "raxml" in your path.  This will allow "raxml" to be found first which will call AVX2 version of RAxML
    try:
        subprocess.call("raxml", stdout=open(os.devnull, 'wb'))
        sys_raxml = "raxml"
        #print("%s found" % sys_raxml)
    except OSError:
        print("looking for RAxML")
        try:
            subprocess.call("raxmlHPC-PTHREADS")
            sys_raxml = "raxmlHPC-PTHREADS"
            print("%s found" % sys_raxml)
        except OSError:
            try:
                subprocess.call("raxmlHPC-SSE3")
                sys_raxml = "raxmlHPC-SSE3"
                print("%s found" % sys_raxml)
            except OSError:
                print("looking for RAxML")
                try:
                    subprocess.call("raxmlHPC")
                    sys_raxml = "raxmlHPC"
                    print("RAxML found")
                except OSError as ex:
                    # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), "#####RAxML is not in you PATH")
                    sys.exit(0)
    arg_options['sys_raxml'] = sys_raxml
    print("\n\n----> RAxML found in $PATH as: %s <-----" % arg_options['sys_raxml'])
    if arg_options['cpu_count'] < 20:
        raxml_cpu = 2
    else:
        raxml_cpu = int(arg_options['cpu_count'] / 10)
    arg_options['raxml_cpu'] = raxml_cpu

    all_parameters = Get_Specie_Parameters() # Class of possible parameters
    print("Sample will be ran as {}" .format(arg_options['species']))
    parameters, genotype_codes = all_parameters.choose(arg_options['species'])
    if parameters['qual_threshold'] is None:
        print("### See species_selection_step2 function")
        sys.exit(0)
    arg_options.update(parameters)

    htmlfile_name = arg_options['root_dir'] + "/summary_log.html"
    arg_options['htmlfile_name'] = htmlfile_name
    htmlfile = open(htmlfile_name, 'at')

    startTime = datetime.now()
    print("Start time: %s" % startTime)

    # DIRECTORY TEST AND BACKUP
    if getattr(sys, 'frozen', False):
        script_used = os.path.realpath(sys.executable)
    elif __file__:
        script_used = os.path.realpath(__file__)

    # make backup
    os.makedirs('starting_files')
    all_starting_files = glob.glob('*vcf')
    for i in all_starting_files:
        shutil.copy(i, 'starting_files')

    test_duplicate()

    print("\ndefiningSNPs: %s " % arg_options['definingSNPs'])
    print("filter_file: %s " % arg_options['filter_file'])
    print("remove_from_analysis: %s " % arg_options['remove_from_analysis'])
    print("step2_upload: %s \n" % arg_options['step2_upload'])

    if genotype_codes:
        print("\nUpdating VCF file names")
        arg_options = change_names(arg_options, genotype_codes)
        malformed = arg_options['malformed']
        names_not_changed = arg_options['names_not_changed']
    else:
        print("Genotypingcode file unavailable.  VCF file names not updated")
        names_not_changed = glob.glob("*.vcf")
        arg_options['malformed'] = []
        arg_options['names_not_changed'] = []

    malformed = arg_options['malformed']
    names_not_changed = arg_options['names_not_changed']

    files = glob.glob('*vcf')
    print("REMOVING FROM ANALYSIS...")
    wb = xlrd.open_workbook(arg_options['remove_from_analysis'])
    ws = wb.sheet_by_index(0)
    for each_sample in ws.col_values(0):
        each_sample = str(each_sample)
        each_sample = re.sub(r'(.*?)[._].*', r'\1', each_sample)
        #print("each sample %s" % each_sample)
        myregex = re.compile(each_sample + '.*') # create regular expression to search for in VCF list
        #print("myregex %s" % myregex)
        for i in files:
            if myregex.search(i):
                print("### --> %s removed from the analysis" % i)
                #print(files)
                #print("\n<h4>### --> %s removed from the analysis</h4>" % i, file=htmlfile)
                try:
                    os.remove(i)
                except FileNotFoundError:
                    print("FileNotFoundError:")
    vcf_starting_list = glob.glob("*.vcf")

    print("CHECKING FOR EMPTY FILES...")
    for filename in vcf_starting_list:
        if os.stat(filename).st_size == 0:
            print("### %s is an empty file and has been deleted" % filename)
            malformed.append("File was empty %s" % filename)
            os.remove(filename)

    all_starting_files = glob.glob('*vcf')
    file_number = len(all_starting_files)

    print("SORTING FILES...")
    defining_snps = {}
    inverted_position = {}
    wb = xlrd.open_workbook(arg_options['definingSNPs'])
    ws = wb.sheet_by_index(0)

    if arg_options['only_all_vcf']:
        print("Only running an All_VCF tree")
    else:
        print("Grouping files...")
        for rownum in range(ws.nrows):
            position = ws.row_values(rownum)[1:][0]
            grouping = ws.row_values(rownum)[:1][0]
            # inverted positions will NOT be found in the passing positions
            # inverted positions are indicated in Defining SNPs by ending with "!"
            if position.endswith('!'):
                position = re.sub('!', '', position)
                inverted_position.update({position: grouping})
            else:
                defining_snps.update({position: grouping})
    files = glob.glob('*vcf')

    arg_options['inverted_position'] = inverted_position
    arg_options['defining_snps'] = defining_snps

    all_list_amb = {}
    group_calls_list = []
  
    if arg_options['debug_call'] and not arg_options['get']:
        for i in files:
            dict_amb, group_calls, mal = group_files(i, arg_options)
            all_list_amb.update(dict_amb)
            group_calls_list.append(group_calls)
            malformed.append(mal)
    else:
        with futures.ProcessPoolExecutor() as pool:
            for dict_amb, group_calls, mal in pool.map(group_files, files, itertools_repeat(arg_options)):
                all_list_amb.update(dict_amb)
                group_calls_list.append(group_calls) # make list of list
                malformed.append(mal)
    malformed = [x for x in malformed if x] # remove empty sets from listn

    print("Getting directory list\n")
    directory_list = next(os.walk('.'))[1] # get list of subdirectories
    directory_list.remove('starting_files')

    print("Placing positions to filter into dictionary...")
    filter_dictionary = get_filters(arg_options)
    arg_options['filter_dictionary'] = filter_dictionary

    if arg_options['gbk_file'] and not arg_options['no_annotation']:
        print("Putting gbk into indexed dataframe...")
        annotation_dict = {}
        for gbk in arg_options['gbk_file']:
            gbk_dict = SeqIO.to_dict(SeqIO.parse(gbk, "genbank"))
            gbk_chrome = list(gbk_dict.keys())[0]
            write_out = open('temp.csv', 'w+')
            for key, value in gbk_dict.items():
                for feature in value.features:
                    if "CDS" in feature.type or "rRNA" in feature.type:
                        myproduct = None
                        mylocus = None
                        mygene = None
                        try:
                            myproduct = feature.qualifiers['product'][0]
                        except KeyError:
                            pass
                        try:
                            mylocus = feature.qualifiers['locus_tag'][0]
                        except KeyError:
                            pass
                        try:
                            mygene = feature.qualifiers['gene'][0]
                        except KeyError:
                            pass
                        print(key, int(feature.location.start), int(feature.location.end), mylocus, myproduct, mygene, sep='\t', file=write_out)
            write_out.close()

            df = pd.read_csv('temp.csv', sep='\t', names=["chrom", "start", "stop", "locus", "product", "gene"])
            os.remove('temp.csv')
            df = df.sort_values(['start', 'gene'], ascending=[True, False])
            df = df.drop_duplicates('start')
            pro = df.reset_index(drop=True)
            pro.index = pd.IntervalIndex.from_arrays(pro['start'], pro['stop'], closed='both')
            annotation_dict[gbk_chrome] = pro
        arg_options['annotation_dict'] = annotation_dict

    samples_in_output = []
    print("Getting SNPs in each directory")
    if arg_options['debug_call']:
        for i in directory_list:
            samples_in_fasta = get_snps(i, arg_options)
            samples_in_output.append(samples_in_fasta)
    else:
        cpu_restriction = int(arg_options['cpu_count'] / 2)
        if cpu_restriction < 1:
            cpu_restriction = 2
        with futures.ProcessPoolExecutor(max_workers=cpu_restriction) as pool:
            for samples_in_fasta in pool.map(get_snps, directory_list, itertools_repeat(arg_options), chunksize=5):
                samples_in_output.append(samples_in_fasta)

    arg_options.pop('filter_dictionary', None) # filters no longer need, get rid of them to make arg_option more managable.

    flattened_list = [item for sublist in samples_in_output for item in sublist]
    flattened_list = set(flattened_list)

    count_flattened_list = len(flattened_list)
    count_vcf_starting_list = len(vcf_starting_list)
    start_end_file_diff_count = count_vcf_starting_list - count_flattened_list

    pretext_flattened_list = get_pretext_list(flattened_list)
    pretext_vcf_starting_list = get_pretext_list(vcf_starting_list)
    pretext_vcf_starting_list = set(pretext_vcf_starting_list)
    try:
        pretext_flattened_list.remove('root')
    except ValueError as ex:
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), "Defining SNPs needed.  If there are no defining SNP then rerun using -a option")
        exit(0)
    difference_start_end_file = pretext_vcf_starting_list.symmetric_difference(pretext_flattened_list)
    difference_start_end_file = list(difference_start_end_file)
    difference_start_end_file.sort()

    # Zip dependency files
    dependents_dir = arg_options['root_dir'] + "/dependents"
    os.makedirs(dependents_dir)
    shutil.copy(arg_options['definingSNPs'], dependents_dir)
    shutil.copy(arg_options['filter_file'], dependents_dir)
    zip(dependents_dir, dependents_dir)
    shutil.rmtree(dependents_dir)

    # remove empty list elements
    arg_options['malformed'] = [x for x in arg_options['malformed'] if x]
    arg_options['names_not_changed'] = [x for x in arg_options['names_not_changed'] if x]
    #############################################
    #MAKE HTML FILE:
    print("<html>\n<head><style> table { font-family: arial, sans-serif; border-collapse: collapse; width: 40%; } td, th { border: 1px solid #dddddd; padding: 4px; text-align: left; font-size: 11px; } </style></head>\n<body style=\"font-size:12px;\">", file=htmlfile)
    print("<h2>Script ran using <u>%s</u> variables</h2>" % arg_options['species'].upper(), file=htmlfile)
    print("<h4>There are %s VCFs in this run</h4>" % file_number, file=htmlfile)

    #OPTIONS
    print(f"Additional options ran: email: {arg_options['email_list']}, filter: {arg_options['filter_finder']}, all_vcf: {arg_options['all_vcf']}, elite: {arg_options['elite']}, no annotation: {arg_options['no_annotation']}, debug: {arg_options['debug_call']}, get: {arg_options['get']}, uploaded: {arg_options['upload']}, ignore filters: {arg_options['ignore_filters']}", file=htmlfile)
    if arg_options['all_vcf']:
        print("\n<h4>All_VCFs is available</h4>", file=htmlfile)
    elif arg_options['elite']:
        print("\n<h4>Elite VCF comparison available</h4>", file=htmlfile)

    #TIME
    print("\n<h4>Start time: %s <br>" % startTime, file=htmlfile)
    print("End time: %s <br>" % datetime.now(), file=htmlfile)
    runtime = (datetime.now() - startTime)
    print("Total run time: %s: </h4>" % runtime, file=htmlfile)

    # ERROR LIST
    if len(arg_options['malformed']) < 1:
        print("<h2>No corrupt VCF removed</h2>", file=htmlfile)
    else:
        print("\n<h2>Corrupt VCF removed</h2>", file=htmlfile)
        for i in arg_options['malformed']:
            print("%s <br>" % i, file=htmlfile)
        print("<br>", file=htmlfile)

    # AMBIGIOUS DEFINING SNPS
    if len(all_list_amb) < 1:
        print("\n<h2>No ambiguous defining SNPs</h2>", file=htmlfile)
    else:
        print("\n<h2>Defining SNPs are ambiguous.  They may be mixed isolates.</h2>", file=htmlfile)
        print("<table>", file=htmlfile)
        print("<tr align=\"left\"><th>Sample Name</th><th>Division</th><th>Absolute Position</th><tr>", file=htmlfile)
        ordered_all_list_amb = OrderedDict(sorted(all_list_amb.items()))
        for k, v in ordered_all_list_amb.items():
            k_split = k.split('\t')
            print("<tr><td>%s</td><td>%s</td><td>%s</td></tr>" % (k_split[0], k_split[1], v), file=htmlfile)
        print("</table>", file=htmlfile)
        print("<br>", file=htmlfile)

    #GROUPING TABLE
    print("<h2>Groupings</h2>", file=htmlfile)
    print("<table>", file=htmlfile)
    print("<tr align=\"left\"><th>Sample Name</th><tr>", file=htmlfile)

    group_calls_list = list(filter(None, group_calls_list))
    try:
        group_calls_list.sort(key=lambda x: x[0]) # sort list of list by first element
    except IndexError as ex:
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), "Unable to sort grouping list")
        pass

    for i in group_calls_list:
        print("<tr>", file=htmlfile)
        for x in i:
            print("<td>%s</td>" % x, end='\t', file=htmlfile)
        print("</tr>", file=htmlfile)
    print("</table>", file=htmlfile)

    # REPORT DIFFERENCES BETWEEN STARTING FILES AND ENDING FILES REPRESENTED IN ALIGNMENTS AND TABLES
    if start_end_file_diff_count < 1:
        print("\n<h2>No files dropped from the analysis.  Input files are equal to those represented in output.</h2>", file=htmlfile)
    else:
        print("\n<h2>{} files have been dropped.  They either need a group, mixed and not finding a group or an error occured.</h2>" .format(start_end_file_diff_count), file=htmlfile)
        print("<table>", file=htmlfile)
        print("<tr align=\"left\"><th>Sample Name</th><tr>", file=htmlfile)
        for i in difference_start_end_file:
            print("<tr><td>{}</td></tr>" .format(i), file=htmlfile)
        print("</table>", file=htmlfile)
        print("<br>", file=htmlfile)
    #Capture program versions for step 2
    try:
        print("\n<h2>Program versions:</h2>", file=htmlfile)
        versions = os.popen('conda list biopython | grep -v "^#"; \
        conda list numpy | egrep -v "^#|numpydoc"; \
        conda list pandas | grep -v "^#"; \
        conda list pysam | grep -v "^#"; \
        conda list pyvcf | grep -v "^#"; \
        conda list xlrd | grep -v "^#"; \
        conda list xlsxwriter | grep -v "^#"; \
        conda list raxml | grep -v "^#"').read()
        versions = versions.split('\n')
        for i in versions:
            print("%s<br>" % i, file=htmlfile)
    except:
        logging.debug("Unable to capture versions")
        pass
    print("Dependent source:  {}<br>" .format(arg_options['script_dependents']), file=htmlfile)

    #FILES NOT RENAMED
    if names_not_changed is None:
        print("\n<h2>File names did not get changed:</h2>", file=htmlfile)
        for i in sorted(names_not_changed):
            print("%s<br>" % i, file=htmlfile)

    print("</body>\n</html>", file=htmlfile)
    #############################################
    os.chdir(arg_options['root_dir'])
    print("Zipping files...")
    zip("starting_files", "starting_files") # zip starting files directory
    shutil.rmtree("starting_files")

    htmlfile.close()

    print("\n\nruntime: %s:  \n" % runtime)

    if arg_options['email_list'] is None:
        print("\n\tEmail not sent")
    elif arg_options['email_list']:
        send_email_step2(arg_options)
        print("\n\tEmail sent to: {}" .format(arg_options['email_list']))
    else:
        print("\n\tEmail not sent")

    if arg_options['upload']:
        print("Uploading Samples...")

        def copytree(src, dst, symlinks=False, ignore=None): #required to ignore permissions
            try:
                for item in os.listdir(src):
                    s = os.path.join(src, item)
                    d = os.path.join(dst, item)
                    try:
                        if os.path.isdir(s):
                            shutil.copytree(s, d, symlinks, ignore)
                        else:
                            shutil.copy2(s, d)
                    except shutil.Error:
                        pass
            except FileNotFoundError:
                print("except FileNotFoundError: file not found")

        #upload to bioinfoVCF
        src = arg_options['root_dir']
        try:
            dst = arg_options['step2_upload'] + "/" + os.path.basename(os.path.normpath(arg_options['root_dir']))
            print("\n\t%s is copying to %s" % (src, dst))
            os.makedirs(dst, exist_ok=True)
            copy_tree(src, dst, preserve_mode=0, preserve_times=0)
            print("Samples were uploaded to {}" .format(dst))
        except TypeError as ex:
            pass
            # debug_log(ex, inspect.getframeinfo(inspect.currentframe()), "No place to upload, check parameters")
    else:
        # logging.debug("Samples were not copied or uploaded to additional location")
        print("\tSamples were not copied or uploaded to additional location")

    print("\n\tComparisons have been made with no obvious error.\n")


def send_email_step2(arg_options):
    htmlfile_name = arg_options['htmlfile_name']
    email_list = arg_options['email_list']

    msg = MIMEMultipart()
    msg['From'] = "tod.p.stuber@aphis.usda.gov"
    msg['To'] = email_list
    msg['Subject'] = "Script 2 " + arg_options['species']
    with open(htmlfile_name) as fp:
        msg.attach(MIMEText(fp.read(), 'html'))

    part = MIMEBase('application', "octet-stream")
    part.set_payload(open("summary_log.html", "r").read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition', 'attachment; filename="summary_log.html"')
    msg.attach(part)

    smtp = smtplib.SMTP('10.10.8.12')
    smtp.send_message(msg)
    smtp.quit()


def get_pretext_list(in_list):
    outlist = []
    for i in in_list:
        pretext = re.sub('[_.].*', '', i)
        outlist.append(pretext)
    return outlist


# def flatten(l):
#     for el in l:
#         if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
#             yield from flatten(el)
#         else:
#             yield el


def zip(src, dst):
    zf = zipfile.ZipFile("%s.zip" % (dst), "w", zipfile.ZIP_DEFLATED)
    abs_src = os.path.abspath(src)
    for dirname, subdirs, files in os.walk(src):
        for filename in files:
            absname = os.path.abspath(os.path.join(dirname, filename))
            arcname = absname[len(abs_src) + 1:]
            zf.write(absname, arcname)
    zf.close()


def test_duplicate():
    dup_list = []
    list_of_files = glob.glob('*vcf')
    for line in list_of_files:
        line = re.sub(r'(.*)[_.].*', r'\1', line)
        dup_list.append(line)
    # find duplicates in list
    duplicates = [k for k, v in Counter(dup_list).items() if v > 1]
    if len(duplicates) > 0:
        print("Duplicates Found: %s " % duplicates)
        print("\n***Error:  Duplicate VCFs")
        sys.exit(0)
    else:
        pass


def change_names(arg_options, genotype_codes):
    
    names_not_changed = []
    list_of_files = glob.glob('*vcf')
    name_found = False
    for filename in list_of_files:
        each_vcf = filename.replace("‐", "-")
        vcf_pretext = re.sub(r'(.*?)[._].*', r'\1', each_vcf) # ? was needed to make greedy, in my view the regex was searching right to left without it.
        vcf_pretext = vcf_pretext.rstrip()
        #Added '^' because h37 18-2397 was finding bovis 18-011018-2397, 2018-06-19
        myregex = re.compile('^' + vcf_pretext + '_.*')  #underscore required to make myregex.search below greedy.  so it finds exact match and not all matches. ex: 10-01 must match 10-01 not 10-010 also
        name_found = False
        try:
            prename = filename.replace(".vcf", "")
            foundname = genotype_codes[prename]
            name_found = True
        except KeyError:
            for key, value in genotype_codes.items():
                try:
                    if myregex.search(key):
                        name_found = True
                        foundname = key.strip('_')
                except TypeError:
                    pass
        if name_found:
            os.rename(filename, foundname + ".vcf")
            print("Name Changed {} --> {}" .format(filename, foundname + ".vcf"))
        else:
            os.rename(filename, each_vcf)
            names_not_changed.append(each_vcf)
            print("File NOT Changed: {} --> {}" .format(filename, each_vcf))
    names_not_changed = set(names_not_changed) # remove duplicates
    arg_options['names_not_changed'] = names_not_changed

    if arg_options['elite']:
        list_of_files = []
        list_of_files = glob.glob('*vcf')
        if not os.path.exists("temp_hold"):
            print("making temp_hold directory")
            os.makedirs("temp_hold") # make all_vcfs if none exists
        for each_vcf in list_of_files:
            # Default 1 * 24 * 60 *60
            time_test = time.time() - os.path.getmtime(each_vcf) < (1 * 24 * 60 * 60) # 1day * (24*60*60)sec in day
            print("%s each_vcf" % each_vcf)
            vcf_pretext = re.sub(r'(.*?)[._].*', r'\1', each_vcf) # ? was needed to make greedy, in my view the regex was searching right to left without it.
            vcf_pretext = vcf_pretext.rstrip()
            myregex = re.compile(vcf_pretext + '.*')
            if time_test:
                print("time_test true %s" % each_vcf)
                shutil.copy(each_vcf, "temp_hold")
            else:
                for k, v in genotype_codes.items():
                    if myregex.search(k):
                        try:
                            print("##### %s" % time_test)
                            if v == "Yes": # if marked yes in column 2 of genotyping codes
                                print("marked yes %s" % each_vcf)
                                shutil.copy(each_vcf, "temp_hold") # if "Yes" then moved to temp_hold
                            else:
                                print("file will be discarded %s" % each_vcf)
                        except FileNotFoundError:
                            print("except FileNotFoundError %s" % each_vcf)
                os.remove(each_vcf)
        shutil.rmtree('starting_files')
        os.makedirs('starting_files')
        os.renames('temp_hold', 'starting_files')
        list_of_files = glob.glob('starting_files/*vcf')
        file_number = len(list_of_files) # update the file_number to present on summary
        for each_vcf in list_of_files:
            shutil.copy(each_vcf, arg_options['root_dir'])
        print(file_number)

    return arg_options


def get_filters(arg_options):
    #get first header to apply all filters to vcf
    worksheet = pd.read_excel(arg_options['filter_file'])
    arg_options["first_column_header"] = worksheet.dtypes.index[0]
    filter_dictionary = defaultdict(list) # key: group_name, values: expanded list
    wb = xlrd.open_workbook(arg_options['filter_file'])
    sheets = wb.sheet_names()
    for sheet in sheets:
        ws = wb.sheet_by_name(sheet)
        for colnum in range(ws.ncols): # for each column in worksheet
            group_name = ws.col_values(colnum)[0] # column header naming file
            mylist = ws.col_values(colnum)[1:] # list of each field in column, minus the header
            mylist = [x for x in mylist if x] # remove blank cells
            for value in mylist:
                value = str(value)
                value = value.replace(sheet + "-", '')
                if "-" not in value:
                    value = int(float(value)) # change str to float to int
                    filter_dictionary[group_name].append(str(sheet) + "-" + str(value))
                elif "-" in value:
                    value = value.split("-")
                    for position in range(int(value[0]), int(value[1]) + 1):
                        filter_dictionary[group_name].append(str(sheet) + "-" + str(position))
    return(filter_dictionary)


def get_read_mean(rec):
    mean_q = int(mean(rec.letter_annotations['phred_quality']))
    return mean_q


def find_filter_dict(each_vcf):
    dict_qual = {}
    dict_map = {}
    vcf_reader = vcf.Reader(open(each_vcf, 'r'))
    for record in vcf_reader:
        try:
            # Freebayes VCFs place MQ values are placed into a list.  GATK as a float
            record.INFO['MQ'] = record.INFO['MQ'][0]
        except TypeError:
            pass
        except KeyError:
            pass
        absolute_positon = str(record.CHROM) + "-" + str(record.POS)
        try:
            returned_qual = []
            returned_map = []
            if int(record.QUAL) > 0:
                returned_qual.append(record.QUAL)
                returned_map.append(record.INFO['MQ'])
                dict_qual[absolute_positon] = returned_qual
                dict_map[absolute_positon] = returned_map
        except Exception:
            pass
    return dict_qual, dict_map


def find_positions(filename, arg_options):
    found_positions = {}
    vcf_reader = vcf.Reader(open(filename, 'r'))
    try:
        for record in vcf_reader:
            try:
                # Freebayes VCFs place MQ values are placed into a list.  GATK as a float
                record.INFO['MQ'] = record.INFO['MQ'][0]
            except TypeError:
                pass
            except KeyError:
                pass
            chrom = record.CHROM
            position = record.POS
            absolute_positon = str(chrom) + "-" + str(position)
            # Usable positins are those that:
            # ADD PARAMETERS HERE TO CHANGE WHAT'S SNP WILL BE USED
            # IF NOT FOUND HERE THE SNP WILL BE IGNORED.  WILL NOT BE REPRESENTED.  HARD REMOVAL
            # parameters
            # str(record.ALT[0]) != "None" --> filter deletions
            # len(record.REF) == 1 --> filter bad ref call with 2 nt present
            # len(record.ALT[0]) == 1 --> filter bad alt call with 2 nt present
            # record.heterozygosity == 0.0 --> filter AC=1, heterozygosity.
            # record.QUAL > 150 --> filter poor quality
            # record.INFO['MQ'] --> filter low map quality
            try:
                if arg_options['species'] == 'flu':
                    # use both AC=1 and AC=2 as valid position
                    if str(record.ALT[0]) != "None" and len(record.REF) == 1 and record.QUAL > arg_options['qual_threshold']:
                        found_positions.update({absolute_positon: record.REF})
                else:
                    if str(record.ALT[0]) != "None" and record.INFO['AC'][0] == 2 and len(record.REF) == 1 and record.QUAL > arg_options['qual_threshold'] and record.INFO['MQ'] > 56:
                        found_positions.update({absolute_positon: record.REF})
            except KeyError as ex:
                # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), filename, absolute_positon)
                pass
    except ZeroDivisionError as ex:
        pass
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), filename, absolute_positon)
    except ValueError as ex:
        pass
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), filename, absolute_positon)
    except UnboundLocalError as ex:
        pass
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), filename, absolute_positon)
    except TypeError as ex:
        pass
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), filename, absolute_positon)
    return found_positions


def check_label(chromesome_position, label, annotation_dict):
    chromesome, position = chromesome_position.split('-')
    pattern_check_list = []
    for gbk_chrome, pro in annotation_dict.items():
        if gbk_chrome == chromesome:
            try:
                aaa=pro.iloc[pro.index.get_loc(int(position))][['chrom', 'locus', 'product', 'gene']]
                chrom, name, locus, tag = aaa.values[0]
                pattern = re.compile(label, re.IGNORECASE)
                pattern_check = pattern.search(str(chrom) + str(position) + str(name) + str(locus) + str(tag))
                ind_annotation = (f"{chrom} {position} {name} {locus} {tag}")
                return pattern_check, ind_annotation
            except KeyError:
                return None


def get_snps(directory, arg_options):

    time_mark = datetime.fromtimestamp(time.time()).strftime('D%Y%m%d_%H%M')

    os.chdir(arg_options['root_dir'] + "/" + directory)
    print("\n----------------------------")
    print("\nworking on: %s " % directory)
    outdir = str(os.getcwd()) + "/"

    filter_dictionary = arg_options['filter_dictionary']
    first_column_header = arg_options["first_column_header"]

    files = glob.glob('*vcf')
    all_positions = {}
    if arg_options['debug_call'] and not arg_options['get']:
        for i in files:
            found_positions = find_positions(i, arg_options)
            all_positions.update(found_positions)
    else:
        with futures.ProcessPoolExecutor() as pool:
            for found_positions in pool.map(find_positions, files, itertools_repeat(arg_options)):
                all_positions.update(found_positions)

    print("Directory %s found positions %s" % (directory, len(all_positions)))
    presize = len(all_positions)

    # Filter applied to all positions
    if not arg_options['ignore_filters']: #and not arg_options['label']:
        try:
            for pos in filter_dictionary[first_column_header]: #filter_list
                all_positions.pop(pos, None)
        except KeyError as ex:
            # Allow keyerror if group is not represented in filter worksheet
            # debug_log(ex, inspect.getframeinfo(inspect.currentframe()), f'Not in filter worksheet: {first_column_header}')
            pass

    # Filter applied to group
    if not arg_options['ignore_filters']: #and not arg_options['label']:
        try:
            for pos in filter_dictionary[directory]: #filter_list
                all_positions.pop(pos, None)
        except KeyError as ex:
            # debug_log(ex, inspect.getframeinfo(inspect.currentframe()), f'Not in filter worksheet: {directory}')
            pass

    if arg_options['label']:
        filter_label_list = []
        label = arg_options['label']
        label = label.replace('/', '|')
        label = label.replace('-', '|')
        label = label.replace(' ', '|')
        label = label.replace(',', '|')
        label_filter = open("label_filter.txt", 'w')
        print(f"Annotation labels filtered: {label}", file=label_filter)
        for chromesome_position in all_positions:
            try:
                pattern_check, ind_annotation = check_label(chromesome_position, label, arg_options['annotation_dict'])
            except TypeError:
                pass
            if pattern_check:
                print(f"{chromesome_position} {ind_annotation}", file=label_filter)
                filter_label_list.append(chromesome_position)
        for chromesome_position in filter_label_list:
            all_positions.pop(chromesome_position, None)
        label_filter.close()


    print("\nDirectory: {}" .format(directory))
    print("Total positions found: {}" .format(presize))
    print("Possible positions filtered {}" .format(len(filter_dictionary)))
    print("Positions after filtering {}" .format(len(all_positions)))

    if arg_options['filter_finder']:
        #write to files
        positions_to_filter = "positions_to_filter.txt"
        positions_to_filter_details = "positions_to_filter_details.txt"
        good_snps = "good_snps_details.txt"
        write_out_positions = open(positions_to_filter, 'w')
        write_out_details = open(positions_to_filter_details, 'w')
        write_out_good_snps = open(good_snps, 'w')

        files = glob.glob('*vcf')

        #calculate mean/max qual and map at all possible positions
        dd_qual = {}
        dd_map = {}
        if arg_options['debug_call']:
            for each_vcf in files:
                print("working on: %s" % each_vcf)
                dict_qual, dict_map = find_filter_dict(each_vcf)
                keys = set(dd_qual).union(dict_qual)
                no = []
                #make position (key) and qual/maps list (value)
                dd_qual = dict((k, dd_qual.get(k, no) + dict_qual.get(k, no)) for k in keys)
                keys = set(dd_map).union(dict_map)
                no = []
                dd_map = dict((k, dd_map.get(k, no) + dict_map.get(k, no)) for k in keys)
        else:
            with Pool(maxtasksperchild=4) as pool:
                for dict_qual, dict_map in pool.map(find_filter_dict, files, chunksize=8):
                    keys = set(dd_qual).union(dict_qual)
                    no = []
                    dd_qual = dict((k, dd_qual.get(k, no) + dict_qual.get(k, no)) for k in keys)
                    keys = set(dd_map).union(dict_map)
                    no = []
                    dd_map = dict((k, dd_map.get(k, no) + dict_map.get(k, no)) for k in keys)

        #dict_qual=dict((k, v) for k, v in dict_qual.items() if v)
        #dict_map=dict((k, v) for k, v in dict_map.items() if v)

        ave_qual = {}
        max_qual = {}
        for k, v in dd_qual.items():
            #only use if > 3 positions have been called
            if len(v) > 3:
                ave_qual[k] = np.mean(v)
                max_qual[k] = np.max(v)

        #provides dictionary as key -> absolute poisiton, value -> average qual/map
        ave_map = {}
        max_map = {}
        for k, v in dd_map.items():
            if len(v) > 3:
                ave_map[k] = np.mean(v)
                max_map[k] = np.max(v)

        # get all possible used positions
        all_maybe_filter = []
        for k in ave_qual.keys():
            all_maybe_filter.append(k)
        for k in max_qual.keys():
            all_maybe_filter.append(k)
        for k in ave_map.keys():
            all_maybe_filter.append(k)
        for k in max_map.keys():
            all_maybe_filter.append(k)
            # remove duplicates
            all_maybe_filter = list(set(all_maybe_filter))

        # Removing those already from all positions to filter
        if arg_options['filter_file']:
            for pos in filter_dictionary[first_column_header]: #filter_list
                try:
                    all_maybe_filter.remove(pos)
                except ValueError:
                    pass
            # Removing those already being filtered for specific group
            try:
                for pos in filter_dictionary[directory]: #filter_list
                    try:
                        all_maybe_filter.remove(pos)
                    except ValueError:
                        pass
            except KeyError as ex:
                pass
                # debug_log(ex, inspect.getframeinfo(inspect.currentframe()), f'Not in filter worksheet: {directory}')
        # for each possible posible position check if to filter.
        for absolute_positon in all_maybe_filter:
            ave_qual_value = ave_qual[absolute_positon]
            max_qual_value = max_qual[absolute_positon]
            ave_map_value = ave_map[absolute_positon]
            max_map_value = max_map[absolute_positon]
            print("%s, max_qual_value: %s, ave_qual_value: %s, max_map_value: %s, ave_map_value: %s" % (absolute_positon, max_qual_value, ave_qual_value, max_map_value, ave_map_value))
            if max_qual_value < 1300 and ave_qual_value < 700 or ave_map_value < 56:
                print("%s, max_qual_value: %s, ave_qual_value: %s, max_map_value: %s, ave_map_value: %s" % (absolute_positon, max_qual_value, ave_qual_value, max_map_value, ave_map_value), file=write_out_details)
                print(absolute_positon, file=write_out_positions)
            else:
                print("%s, max_qual_value: %s, ave_qual_value: %s, max_map_value: %s, ave_map_value: %s" % (absolute_positon, max_qual_value, ave_qual_value, max_map_value, ave_map_value), file=write_out_good_snps)
        write_out_positions.close()
        write_out_details.close()
        write_out_good_snps.close()

    table_location = outdir + directory + "-table.txt"
    table = open(table_location, 'wt')

    # write absolute positions to table
    # order before adding to file to match with ordering of individual samples below
    # all_positions is abs_pos:REF
    all_positions = OrderedDict(sorted(all_positions.items()))
    # Add the positions to the table
    print("reference_pos", end="\t", file=table)
    for k, v in all_positions.items():
        print(k, end="\t", file=table)
    print("", file=table)

    list_of_files = glob.glob('*vcf')

    # for each vcf
    all_map_qualities = {}
    for file_name in list_of_files:
        sample_map_qualities = {}
        just_name = file_name.replace('.vcf', '')
        just_name = re.sub('\..*', '*', just_name) # if after the .vcf is removed there is stilll a "." in the name it is assumed the name did not get changed
        print(just_name, end="\t", file=table)
        # for each line in vcf
        vcf_reader = vcf.Reader(open(file_name, 'r'))
        sample_dict = {}
        for record in vcf_reader:
            try:
                # Freebayes VCFs place MQ values are placed into a list.  GATK as a float
                record.INFO['MQ'] = record.INFO['MQ'][0]
            except TypeError:
                pass
            except KeyError:
                pass
            record_position = str(record.CHROM) + "-" + str(record.POS)
            if record_position in all_positions:
                #print("############, %s, %s" % (file_name, record_position))
                # NOT SURE THIS IS THE BEST PLACE TO CAPTURE MQ AVERAGE
                # MAY BE FASTER AFTER PARSIMONY SNPS ARE DECIDED, BUT THEN IT WILL REQUIRE OPENING THE FILES AGAIN.
                if str(record.ALT[0]) != "None" and str(record.INFO['MQ']) != "nan": #on rare occassions MQ gets called "NaN" thus passing a string when a number is expected when calculating average.
                    #print("getting map quality:    %s          %s      %s" % (record.INFO['MQ'], file_name, str(record.POS)))
                    sample_map_qualities.update({record_position: record.INFO['MQ']})
                # ADD PARAMETERS HERE TO CHANGE WHAT'S EACH VCF REPRESENTS.
                # SNP IS REPRESENTED IN TABLE, NOW HOW WILL THE VCF REPRESENT THE CALLED POSITION
                # str(record.ALT[0]) != "None", which means a deletion as ALT
                # not record.FILTER, or rather PASSED.
                # check record.QUAL
                # In GATK VCFs "!= None" not used.
                if str(record.ALT[0]) != "None" and len(record.ALT[0]) == 1 and record.INFO['AC'][0] == 2 and record.QUAL > arg_options['N_threshold']:
                    sample_dict.update({record_position: record.ALT[0]})
                elif str(record.ALT[0]) != "None" and len(record.ALT[0]) == 1 and record.INFO['AC'][0] == 1 and int(record.QUAL) > arg_options['N_threshold']:
                    ref_alt = str(record.ALT[0]) + str(record.REF[0])
                    if ref_alt == "AG":
                        sample_dict.update({record_position: "R"})
                    elif ref_alt == "CT":
                        sample_dict.update({record_position: "Y"})
                    elif ref_alt == "GC":
                        sample_dict.update({record_position: "S"})
                    elif ref_alt == "AT":
                        sample_dict.update({record_position: "W"})
                    elif ref_alt == "GT":
                        sample_dict.update({record_position: "K"})
                    elif ref_alt == "AC":
                        sample_dict.update({record_position: "M"})
                    elif ref_alt == "GA":
                        sample_dict.update({record_position: "R"})
                    elif ref_alt == "TC":
                        sample_dict.update({record_position: "Y"})
                    elif ref_alt == "CG":
                        sample_dict.update({record_position: "S"})
                    elif ref_alt == "TA":
                        sample_dict.update({record_position: "W"})
                    elif ref_alt == "TG":
                        sample_dict.update({record_position: "K"})
                    elif ref_alt == "CA":
                        sample_dict.update({record_position: "M"})
                    else:
                        sample_dict.update({record_position: "N"})
                    # Poor calls
                elif str(record.ALT[0]) != "None" and int(record.QUAL) <= 50:
                    sample_dict.update({record_position: record.REF[0]})
                elif str(record.ALT[0]) != "None" and int(record.QUAL) <= arg_options['N_threshold']:
                    sample_dict.update({record_position: "N"})
                elif str(record.ALT[0]) != "None": #Insurance -- Will still report on a possible SNP even if missed with above statement
                    sample_dict.update({record_position: str(record.REF[0])})
                elif str(record.ALT[0]) == "None":
                    sample_dict.update({record_position: "-"})

        # After iterating through VCF combine dict to nested dict
        all_map_qualities.update({just_name: sample_map_qualities})

        # merge dictionaries and order
        merge_dict = {}
        merge_dict.update(all_positions) #abs_pos:REF
        merge_dict.update(sample_dict) # abs_pos:ALT replacing all_positions, because keys must be unique
        merge_dict = OrderedDict(sorted(merge_dict.items())) #OrderedDict of ('abs_pos', ALT_else_REF), looks like a list of lists
        for k, v in merge_dict.items():
            #print("k %s, v %s" % (k, v))
            print(str(v) + "\t", file=table, end="")
        print("", file=table) # sample printed to file
    table.close() #end of loop.  All files done

    # Select parsimony informative SNPs
    mytable = pd.read_csv(table_location, sep='\t')
    # drop NaN rows and columns
    mytable = mytable.dropna(axis=1)

    # SELECT PARISOMONY INFORMATIVE SNPSs
    # removes columns where all fields are the same
    parsimony = mytable.loc[:, (mytable != mytable.iloc[0]).any()]
    parsimony_positions = list(parsimony)
    #write over table (table_location) containing all snps
    parsimony.to_csv(table_location, sep="\t", index=False)
    table = open(table_location, 'a')
    # The reference calls are added after the parsimony positions are selected.
    # added corresponding reference to parsimony table
    print("reference_call", end="\t", file=table)
    #all_positions_list=list(all_positions)
    try: #if there is only one file in the group exception is needed to return a value
        parsimony_positions.remove('reference_pos')
    except ValueError:
        samples_in_fasta = []
        return(samples_in_fasta)

    list_of_ref = []
    for abs_pos in parsimony_positions:
        list_of_ref.append(all_positions.get(abs_pos))
    string_of_ref = "\t".join(list_of_ref)
    print(string_of_ref, file=table)
    table.close()

    samples_in_fasta = []
    #Print out fasta alignment file from table
    alignment_file = outdir + directory + "_" + time_mark + ".fasta"
    write_out = open(alignment_file, 'wt')
    with open(table_location, 'rt') as f:
        count = 0
        for line in f:
            if count > 0:
                line = re.sub('^', '>', line)
                line = line.replace('reference_call', 'root')
                line = line.replace('\t', '\n', 1)
                samples_in_fasta.append(line.split('\n')[0].replace('>', ''))
                line = line.replace('\t', '')
                print(line, end="", file=write_out)
            count = count + 1
    write_out.close()

    try: #if there are no SNP is the table
        mytable = pd.read_csv(table_location, sep='\t')
    except:
        samples_in_fasta = []
        return(samples_in_fasta)

    # move reference to top row
    myref = mytable[-1:]
    myother = mytable[:-1]
    frames = [myref, myother]
    mytable = pd.concat(frames)
    mytable.to_csv(table_location, sep="\t", index=False)

    print("\n%s table dimensions: %s" % (directory, str(mytable.shape)))

    print("%s RAxML running..." % directory)
    try:
        if arg_options['only_all_vcf']:
            os.system("{} -s {} -n raxml -m GTRCATI -o root -p 12345 -T {} > /dev/null 2>&1" .format(arg_options['sys_raxml'], alignment_file, arg_options['cpu_count']))
        else:
            os.system("{} -s {} -n raxml -m GTRCATI -o root -p 12345 -T {} > /dev/null 2>&1" .format(arg_options['sys_raxml'], alignment_file, arg_options['raxml_cpu']))
    except:
        logging.warning('RAxML failed')
        write_out = open('RAXML_FAILED', 'w+')
        write_out.close()
        pass
    try:
        ordered_list_from_tree = outdir + directory + "-cleanedAlignment.txt"
        write_out = open(ordered_list_from_tree, 'w+')
        print("reference_pos", file=write_out)
        print("reference_call", file=write_out)
        if os.path.isfile("RAxML_bestTree.raxml"):
            with open("RAxML_bestTree.raxml", 'rt') as f:
                for line in f:
                    line = re.sub('[:,]', '\n', line)
                    line = re.sub('[)(]', '', line)
                    line = re.sub('[0-9].*\.[0-9].*\n', '', line)
                    line = re.sub('root\n', '', line)
                    write_out.write(line)
            best_raxml_tre = directory + "_" + time_mark + "-RAxML-bestTree.tre"
            os.rename("RAxML_bestTree.raxml", best_raxml_tre)
            write_out.close()
        best_raxml_svg = directory + "_" + time_mark + "-RAxML-bestTree.svg"
        try:
            os.system("cat {} | nw_display -s -S -w 1300 -t -v 30 -i 'opacity:0' -b 'opacity:0' -l 'font-size:14;font-family:serif;font-style:italic' -d 'stroke-width:1;stroke:blue' - > {}" .format(best_raxml_tre, best_raxml_svg)) #-s produces svg, -S suppress scale bar, -w to set the number of columns available for display, -t tab format, -v vertical spacing, -i inner node label, -b branch style
        except:
            logging.debug(f'{directory} nw_display failed, likely RAxML related, not enough samples in group to build tree')
            pass
        out_org = str(os.getcwd()) + "/" + directory + "_" + time_mark + "-organized-table.txt"
        out_sort = str(os.getcwd()) + "/" + directory + "_" + time_mark + "-sorted-table.txt"

        sort_table(table_location, ordered_list_from_tree, out_org) #function

        print("%s Getting map quality..." % directory)
        average = lambda x: x.mean()
        all_map_qualities = pd.DataFrame(all_map_qualities)
        #ave_mq = Type: Series
        ave_mq = all_map_qualities.apply(average, axis=1)
        ave_mq = ave_mq.astype(int)
        ave_mq.to_csv('outfile.txt', sep='\t', header='True') # write to csv

        write_out = open('map_quality.txt', 'w+')
        print('reference_pos\tmap-quality', file=write_out)
        with open('outfile.txt', 'rt') as f:
            for line in f:
                write_out.write(line)
        write_out.close()
        #seemed pooling did not like a function with no parameters given
        quality = pd.read_csv('map_quality.txt', sep='\t')

        mytable = pd.read_csv(table_location, sep='\t')
        mytable = mytable.set_index('reference_pos')

        # order list is from tree file
        # gives order for samples to be listed in table to be phylogenetically correct
        ordered_list = []
        with open(ordered_list_from_tree) as infile:
            for i in infile:
                i = i.rstrip()
                ordered_list.append(i)
        # sinces this is set as the mytable index do not include in ordering
        ordered_list.remove('reference_pos')

        # reorder table based on order of list
        mytable = mytable.reindex(ordered_list)
        mytable.to_csv(table_location, sep='\t')

        mytable_sort = pd.read_csv(table_location, sep='\t') #sorted
        mytable_sort = mytable_sort.set_index('reference_pos') #sorted
        mytable_sort = mytable_sort.transpose() #sort
        mytable_sort.to_csv(out_sort, sep='\t', index_label='reference_pos') #sort

        mytable = pd.read_csv(out_org, sep='\t') #org
        mytable = mytable.set_index('reference_pos') #org
        mytable = mytable.transpose() #org
        mytable.to_csv(out_org, sep='\t', index_label='reference_pos') #org

        if arg_options['gbk_file'] and not arg_options['no_annotation']:

            print("{} annotating from annotation dictionary... {}" .format(directory, time_mark))
            mytable_sort = pd.read_csv(out_sort, sep='\t') #sort
            mytable_sort = mytable_sort.merge(quality, on='reference_pos', how='inner')  #sort
            mytable_sort.to_json('mytable_sort.json')

            annotation_dict = arg_options['annotation_dict']
            for gbk_chrome, pro in annotation_dict.items():
                ref_pos = mytable_sort[['reference_pos']]
                ref_pos = ref_pos.rename(columns={'index': 'reference_pos'})
                ref_pos = pd.DataFrame(ref_pos.reference_pos.str.split('-', expand=True).values, columns=['reference', 'position'])
                ref_pos = ref_pos[ref_pos['reference'] == gbk_chrome]

                write_out = open('annotations.csv', 'a')
                positions = ref_pos.position.to_frame()
                for index, row in positions.iterrows():
                    pos = row.position
                    try:
                        aaa = pro.iloc[pro.index.get_loc(int(pos))][['chrom', 'locus', 'product', 'gene']]
                        try:
                            chrom, name, locus, tag = aaa.values[0]
                            print("{}-{}\t{}, {}, {}".format(chrom, pos, locus, tag, name), file=write_out)
                        except ValueError:
                            # if only one annotation entire chromosome (such with flu) then having [0] fails
                            chrom, name, locus, tag = aaa.values
                            print("{}-{}\t{}, {}, {}".format(chrom, pos, locus, tag, name), file=write_out)
                    except KeyError:
                        print("{}-{}\tNo annotated product" .format(gbk_chrome, pos), file=write_out)
                write_out.close()

                annotations_df = pd.read_csv('annotations.csv', sep='\t', header=None, names=['index', 'annotations'], index_col='index')

            annotations_df.index.names = ['reference_pos']
            mytable_sort = mytable_sort.set_index('reference_pos')
            annotations_df.index.names = ['reference_pos']
            mytable_sort = mytable_sort.merge(annotations_df, left_index=True, right_index=True)
            mytable_sort = mytable_sort.transpose() #sort
            mytable_sort.to_csv(out_sort, sep='\t', index_label='reference_pos') #sort

            mytable_org = pd.read_csv(out_org, sep='\t') #org
            mytable_org = mytable_org.merge(quality, on='reference_pos', how='inner') #org
            mytable_org = mytable_org.set_index('reference_pos')
            mytable_org = mytable_org.merge(annotations_df, left_index=True, right_index=True)
            mytable_org = mytable_org.transpose() #org
            mytable_org.to_csv(out_org, sep='\t', index_label='reference_pos') #org

        else:
            print("No gbk file or no table to annotate")
            mytable_sort = pd.read_csv(out_sort, sep='\t') #sort
            mytable_sort = mytable_sort.merge(quality, on='reference_pos', how='inner') #sort
            mytable_sort = mytable_sort.set_index('reference_pos') #sort
            mytable_sort = mytable_sort.transpose() #sort
            mytable_sort.to_csv(out_sort, sep='\t', index_label='reference_pos') #sort
            # add when no annotation
            with open(out_sort, 'rt') as f:
                line = f.readline()
            f.close()
            column_count = line.count('\t') #sort
            column_count = column_count - 1 #sort
            #print("column_count: %s" % column_count)
            with open(out_sort, 'at') as f:
                print("no_annotation", end='', file=f)
                print('\t' * column_count, file=f)
            f.close()

            print("No gbk file or no table to annotate")
            mytable = pd.read_csv(out_org, sep='\t') #org
            mytable = mytable.merge(quality, on='reference_pos', how='inner') #org
            mytable = mytable.set_index('reference_pos') #org
            mytable = mytable.transpose() #org
            mytable.to_csv(out_org, sep='\t', index_label='reference_pos') #org
            # add when no annotation
            with open(out_org, 'rt') as f:
                line = f.readline()
            f.close()
            column_count = line.count('\t')
            column_count = column_count - 1
            #print("column_count: %s" % column_count)
            with open(out_org, 'at') as f:
                print("no_annotation", end='', file=f)
                print('\t' * column_count, file=f)
            f.close()

        excelwriter(out_sort) #***FUNCTION CALL #sort
        excelwriter(out_org) #***FUNCTION CALL #org

        for r in glob.glob('*vcf'):
            os.remove(r)

    except ValueError as ex:
        # warning_log(ex, inspect.getframeinfo(inspect.currentframe()), file_name, 'Possible table creation failure')
        return

    try:
        os.remove(ordered_list_from_tree)
        if arg_options['gbk_file']:
            os.remove("annotations.csv")
        os.remove("outfile.txt")
        os.remove("map_quality.txt")
        os.remove(out_sort)
        os.remove(out_org) # organized.txt table
        os.remove(table_location) # unorganized table
        os.remove('RAxML_info.raxml')
        os.remove('RAxML_log.raxml')
        os.remove('RAxML_parsimonyTree.raxml')
        os.remove('RAxML_result.raxml')
        os.remove(directory + "_" + time_mark + '.fasta.reduced')

    except FileNotFoundError:
        pass

    # PANDA NOTES
    # get the index: mytable.index
    # get columns: mytable.columns
    # get a column: mytable.AF2122_NC002945_105651, shows index (sample names)
    # get a row: mytable.ix['reference'], shows columns (positions and SNPs)
    # values: mytable.values, SNPs - series
    # strip off the bottom row: mytable[:-1]
    # get the bottom row: mytable[-1:]

    with open(directory + "_" + time_mark + "-samples_in_fasta.json", 'w') as outfile:
        json.dump(samples_in_fasta, outfile)

    return(samples_in_fasta)


def sort_table(table_location, ordered, out_org):
    mytable = pd.read_csv(table_location, sep='\t')
    #mytable=mytable.set_index('reference_pos')

    # order list is from tree file
    # gives order for samples to be listed in table to be phylogenetically correct
    ordered_list = []
    with open(ordered) as infile:
        for i in infile:
            i = i.rstrip()
            ordered_list.append(i)

    # Convert reference_pos-column to category and in set the ordered_list as categories hierarchy
    mytable.reference_pos = mytable.reference_pos.astype("category")
    mytable.reference_pos.cat.set_categories(ordered_list, inplace=True)
    mytable = mytable.sort_values(["reference_pos"]) # 'sort' changed to 'sort_values'

    # count number of SNPs in each column
    snp_per_column = []
    for column_header in mytable:
        count = 0
        column = mytable[column_header]
        # for each element in the column
        for element in column:
            if element != column[0]:
                count = count + 1
        snp_per_column.append(count)
        #print("the count is: %s" % count)
    row1 = pd.Series(snp_per_column, mytable.columns, name="snp_per_column")
    #row1 = row1.drop('reference_pos')

    # get the snp count per column
    # for each column in the table
    snp_from_top = []
    for column_header in mytable:
        count = 0
        column = mytable[column_header]
        # for each element in the column
        # skip the first element
        for element in column[1:]:
            if element == column[0]:
                count = count + 1
            else:
                break
        snp_from_top.append(count)
    row2 = pd.Series(snp_from_top, mytable.columns, name="snp_from_top")
    #row2 = row2.drop('reference_pos')

    mytable = mytable.append([row1])
    mytable = mytable.append([row2])
    #In pandas=0.18.1 even this does not work:
    #    abc = row1.to_frame()
    #    abc = abc.T --> mytable.shape (5, 18), abc.shape (1, 18)
    #    mytable.append(abc)
    #Continue to get error: "*** ValueError: all the input arrays must have same number of dimensions"

    mytable = mytable.T
    mytable = mytable.sort_values(['snp_from_top', 'snp_per_column'], ascending=[True, False])
    mytable = mytable.T

    # remove snp_per_column and snp_from_top rows
    mytable = mytable[:-2]
    mytable.to_csv(out_org, sep='\t', index=False)


def excelwriter(filename):
    orginal_name = filename
    filename = filename.replace(".txt", ".xlsx")
    wb = xlsxwriter.Workbook(filename)
    ws = wb.add_worksheet("Sheet1")
    with open(orginal_name, 'r') as csvfile:
        table = csv.reader(csvfile, delimiter='\t')
        i = 0
        for row in table:
            ws.write_row(i, 0, row)
            i += 1

    col = len(row)
    col = col + 1
    #print(i, "x", col)

    formatA = wb.add_format({'bg_color': '#58FA82'})
    formatG = wb.add_format({'bg_color': '#F7FE2E'})
    formatC = wb.add_format({'bg_color': '#0000FF'})
    formatT = wb.add_format({'bg_color': '#FF0000'})
    formatnormal = wb.add_format({'bg_color': '#FDFEFE'})
    formatlowqual = wb.add_format({'font_color': '#C70039', 'bg_color': '#E2CFDD'})
    formathighqual = wb.add_format({'font_color': '#000000', 'bg_color': '#FDFEFE'})
    formatambigous = wb.add_format({'font_color': '#C70039', 'bg_color': '#E2CFDD'})
    formatN = wb.add_format({'bg_color': '#E2CFDD'})

    ws.conditional_format(i - 2, 1, i - 2, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 60, 'format': formathighqual})
    ws.conditional_format(i - 2, 1, i - 2, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 59, 'format': formathighqual})
    ws.conditional_format(i - 2, 1, i - 2, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 58, 'format': formathighqual})
    ws.conditional_format(i - 2, 1, i - 2, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 57, 'format': formathighqual})
    ws.conditional_format(i - 2, 1, i - 2, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 56, 'format': formathighqual})
    ws.conditional_format(i - 2, 1, i - 2, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 55, 'format': formathighqual})
    ws.conditional_format(i - 2, 1, i - 2, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 54, 'format': formathighqual})
    ws.conditional_format(i - 2, 1, i - 2, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 53, 'format': formathighqual})
    ws.conditional_format(i - 2, 1, i - 2, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 52, 'format': formathighqual})
    ws.conditional_format(i - 2, 1, i - 2, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 51, 'format': formathighqual})
    ws.conditional_format(i - 2, 1, i - 2, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 50, 'format': formathighqual})
    ws.conditional_format(i - 2, 1, i - 2, col - 2, {'type': 'text', 'criteria': 'not containing', 'value': 100, 'format': formatlowqual})
    ws.conditional_format(2, 1, i - 3, col - 2, {'type': 'cell', 'criteria': '==', 'value': 'B$2', 'format': formatnormal})
    ws.conditional_format(2, 1, i - 3, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 'A', 'format': formatA})
    ws.conditional_format(2, 1, i - 3, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 'G', 'format': formatG})
    ws.conditional_format(2, 1, i - 3, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 'C', 'format': formatC})
    ws.conditional_format(2, 1, i - 3, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 'T', 'format': formatT})
    ws.conditional_format(2, 1, i - 3, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 'S', 'format': formatambigous})
    ws.conditional_format(2, 1, i - 3, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 'Y', 'format': formatambigous})
    ws.conditional_format(2, 1, i - 3, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 'R', 'format': formatambigous})
    ws.conditional_format(2, 1, i - 3, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 'W', 'format': formatambigous})
    ws.conditional_format(2, 1, i - 3, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 'K', 'format': formatambigous})
    ws.conditional_format(2, 1, i - 3, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 'M', 'format': formatambigous})
    ws.conditional_format(2, 1, i - 3, col - 2, {'type': 'text', 'criteria': 'containing', 'value': 'N', 'format': formatN})
    ws.conditional_format(2, 1, i - 3, col - 2, {'type': 'text', 'criteria': 'containing', 'value': '-', 'format': formatN})

    ws.set_column(0, 0, 30)
    ws.set_column(1, col - 2, 2)
    ws.freeze_panes(2, 1)
    format_rotation = wb.add_format({'rotation': '90'})
    ws.set_row(0, 140, format_rotation)
    formatannotation = wb.add_format({'font_color': '#0A028C', 'rotation': '-90', 'align': 'top'})
    #set last row
    ws.set_row(i - 1, 400, formatannotation)

    wb.close()
