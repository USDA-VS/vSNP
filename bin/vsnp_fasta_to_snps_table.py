#!/usr/bin/env python

__version__ = "2.03"

import os
import shutil
import subprocess
import sys
import re
import glob
import time
import json
import random
import zipfile
from datetime import datetime
from pathlib import Path
from collections import OrderedDict
import numpy as np
import pandas as pd
import multiprocessing
multiprocessing.set_start_method('spawn', True)
from concurrent import futures
import argparse
import textwrap
import vcf
from Bio import SeqIO

# from vsnp_reference_options import Ref_Options
from vsnp_file_management import File_Management
from vsnp_chromosome_reference import Reference_Chromosome
from vsnp_remove_from_analysis import Remove_From_Analysis
from vsnp_html_step2_summary import Step2_Summary


class Get_Snps:
    ''' 
    '''
    def find_raxml(self):
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
                    except OSError:
                        print("#####RAxML is not in you PATH")
                        print("conda install raxml")
                        sys.exit(0)
        return sys_raxml

    def __init__(self, reference=None, working_directory='.', gbk=None, debug=False, all_vcf=False, table=False, qual_threshold=150, MQ=56, AC=2, N_threshold=50):
        if working_directory == '.':
            #default if no path given, assume current working directory
            self.cwd = str(os.getcwd())
        elif working_directory.startswith('~'):
            #path relative to user
            home = str(Path.home())
            self.cwd = working_directory.replace('~', home)
        elif working_directory.startswith('/'):
            #assume full path given
            self.cwd = working_directory
        else:
            print(f'Need a legitimate directory path')
            print(f'path given: {working_directory}')
            sys.exit(0)
        self.vcf_list = glob.glob(f'{self.cwd}/*.vcf')
        self.all_vcf = all_vcf
        self.debug = debug
        self.gbk = gbk
        self.qual_threshold = qual_threshold
        self.MQ = MQ
        self.AC = AC
        self.N_threshold = N_threshold
        self.directory_name = self.cwd.split('/')[-1]
        self.sample_path_name = f'{self.cwd}'
        self.time_mark = datetime.fromtimestamp(time.time()).strftime('D%Y%m%d_%H%M')
        self.cpu_count_half = int(multiprocessing.cpu_count() / 2)  #best speed with 2
        ts = time.time()
        self.st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
        self.htmlsummary = {}
        self.htmlsummary["st"] = self.st
        self.all_mal = []

        script_path = os.path.dirname(os.path.realpath(__file__))

        if working_directory == '.':
            working_directory = os.getcwd()
        else:
            working_directory = working_directory


    def checksum_fastas_for_tree(self, infile):

        unique_number = ''.join(random.choice('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ') for i in range(16))

        uniq_id_fasta = re.sub('\..*', '-uniqid.fasta', infile)

        checksum_dict = {}
        record_iterator = SeqIO.parse(infile, "fasta")
        outfasta = open(uniq_id_fasta , 'at')
        for fasta_file in record_iterator:
            if fasta_file.description == "reference_seq":
                print(">%s\n%s" % (fasta_file.description, fasta_file.seq), file=outfasta)
                checksum_dict.update({fasta_file.description:fasta_file.description})
            else:
                unique_number = ''.join(random.choice('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ') for i in range(16))
                print ("description: %s " % fasta_file.description)
                print("checksum %s " % unique_number)
                if fasta_file.description in checksum_dict.values():
                    dup_header = ''.join(random.choice('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ') for i in range(2))
                    checksum_dict.update({unique_number:fasta_file.description + "-DUPLICATE_HEADER_NAME-" + dup_header})
                else:
                    checksum_dict.update({unique_number:fasta_file.description})
                print (">%s\n%s" % (unique_number, fasta_file.seq), file=outfasta)
        outfasta.close()

        idtable = open("idtable.txt" , 'wt')
        for k, v in checksum_dict.items():
            print("%s\t%s" % (k,v), file=idtable)
        idtable.close()

        return (uniq_id_fasta)

    def checksum_match_to_text(self, tree):
        # read entire tree into variable as string obj
        with open(tree, 'rt') as open_tree:
            entire_file = open_tree.read()
            print (entire_file)
        with open("idtable.txt" , 'rt') as f:
            for line in f:
                line = line.strip('\n')
                line = line.split("\t")
                org = str(line[1])
                org = org.rstrip()
                check = str(line[0])
                check = check.rstrip()
                print (org)
                print (check)
                entire_file = re.sub(check, "'" + org + "'", entire_file)
            f.close()
        outfile = "NAMES-UPDATED-" + os.path.basename(tree)
        write_out = open(outfile , 'wt')
        write_out.write(entire_file)
        write_out.close()
        return outfile

    def df_to_fasta(self, df, alignment_file):
        test_duplicates=[] # if duplicate name in alignment fasta raxml with error and exit
        with open(alignment_file, 'w') as write_out:
            for index, row in df.iterrows():
                test_duplicates.append(row.name)
                if test_duplicates.count(row.name) < 2:
                    print(f'>{row.name}', file=write_out)
                    for pos in row:
                        print(pos, end='', file=write_out)
                    print("", file=write_out)

    def get_parsimonious_pos(self, in_df):
        try:
            ref_series = in_df.loc['reference_seq']
            in_df = in_df.drop(['reference_seq']) #in all_vcf reference_seq needs to be removed
        except KeyError:
            print('Check that there is a "reference_seq" nameed')
            sys.exit(0)
        # print(f'in_df size: {in_df.shape}')
        parsimony = in_df.loc[:, (in_df != in_df.iloc[0]).any()]
        parsimony_positions = list(parsimony)
        parse_df = in_df[parsimony_positions]
        ref_df = ref_series.to_frame()
        ref_df = ref_df.T
        out_df = pd.concat([parse_df, ref_df], join='inner')
        return out_df

    def gather_and_filter(self, alignment, excel_path=None, directory=None):
        try:
            if alignment.endswith('.json'):
                alignment = pd.read_json(alignment)
        except AttributeError:
            if isinstance(alignment, pd.DataFrame):
                pass
            else:
                print(f'### Provide a proper alignment file as dataframe or json.  Must include absolute positions, therefore fasta unexceptable.')
        if directory is None:
            sample_path_name = f'{self.sample_path_name}'
            group = "fasta"
        else:
            sample_path_name = f'{self.sample_path_name}/{directory}'
            group = directory
        raxml = self.find_raxml()
        st = self.st
        if excel_path is None:
            filtered_all_df = alignment
            # ref_series = filtered_all_df.loc['reference_seq']
            sheet_names = None
        elif excel_path:
            #filter positions to be removed from all
            xl = pd.ExcelFile(excel_path)
            sheet_names = xl.sheet_names
            exclusion_list_all = self.get_position_list(sheet_names, excel_path, 0) #Use the first column to filter "all" postions
            exclusion_list_group = self.get_position_list(sheet_names, excel_path, directory) #Use the first column to filter "all" postions
            exclusion_list = exclusion_list_all + exclusion_list_group
            filtered_all_df = alignment.drop(columns=exclusion_list, errors='ignore')  #filters for all applied
            # ref_series = filtered_all_df.loc['reference_seq']
            # print(f'{filtered_all_df.shape} Table size after filtering')
        parsimonious_df = self.get_parsimonious_pos(filtered_all_df)
        parsimonious_json =f'{sample_path_name}/{group}_parsimonious_alignment-{st}.json'
        parsimonious_df.to_json(parsimonious_json, orient='split')
        alignment_file = f'{sample_path_name}/{group}_parsimonious_alignment-{st}.fasta'
        tree_file = f'{sample_path_name}/{group}_parsimonious_alignment-{st}.tre'
        self.df_to_fasta(parsimonious_df, alignment_file)
        samples_number, columns = parsimonious_df.shape
        if samples_number < 4:
            with open(f'{sample_path_name}/TOO_FEW_SAMPLES_TO_BUILD_TREE', 'w') as message_out:
                print(f'check sample numbers', file=message_out)
        else:
            os.system(f'{raxml} -s {alignment_file} -n raxml -m GTRCATI -o reference_seq -w {sample_path_name} -p 456123 -T 4 > /dev/null 2>&1')
            try:
                os.rename(f'{sample_path_name}/RAxML_bestTree.raxml', tree_file)
                raxml_to_remove = glob.glob(f'{sample_path_name}/RAxML*')
                for each in raxml_to_remove:
                    os.remove(each)
                try:
                    os.remove(f'{alignment_file}.reduced')
                except FileNotFoundError:
                    pass
            except FileNotFoundError:
                with open(f'{sample_path_name}/SEE_RAXML_INFO', 'w') as message_out:
                    print(f'check sample numbers', file=message_out)
        return parsimonious_json, tree_file

    def build_tables(self, parsimonious_json, tree, mq=None):
        sample_path_name = self.sample_path_name
        st = self.st
        alignment = pd.read_json(parsimonious_json, orient='split')
        # os.remove(alignment_json)
        with open(tree, 'rt') as tree_file: #must be the single line newick format.  Not Nexus which will be mutliline often with formating
            for line in tree_file:
                line = re.sub('[:,]', '\n', line)
                line = re.sub('[)(]', '', line)
                line = re.sub('[0-9].*\.[0-9].*\n', '', line)
                line = re.sub('reference_seq\n', '', line)
        sample_order = line.split('\n')
        sample_order = list(filter(None, sample_order))
        sample_order.insert(0, 'reference_seq')
        tree_order = alignment.loc[sample_order]
            # count number of SNPs in each column
        snp_per_column = []
        for column_header in tree_order:
            count = 0
            column = tree_order[column_header]
            # for each element in the column
            for element in column:
                if element != column[0]:
                    count = count + 1
            snp_per_column.append(count)
            #print("the count is: %s" % count)
        row1 = pd.Series(snp_per_column, tree_order.columns, name="snp_per_column")
        #row1 = row1.drop('reference_seq')

        # get the snp count per column
        # for each column in the table
        snp_from_top = []
        for column_header in tree_order:
            count = 0
            column = tree_order[column_header]
            # for each element in the column
            # skip the first element
            for element in column[1:]:
                if element == column[0]:
                    count = count + 1
                else:
                    break
            snp_from_top.append(count)
        row2 = pd.Series(snp_from_top, tree_order.columns, name="snp_from_top")
        #row2 = row2.drop('reference_seq')
        tree_order = tree_order.append([row1])
        tree_order = tree_order.append([row2])
        tree_order = tree_order.T
        tree_order = tree_order.sort_values(['snp_from_top', 'snp_per_column'], ascending=[True, False])
        tree_order = tree_order.T

        # remove snp_per_column and snp_from_top rows
        cascade_order = tree_order[:-2]

        ###Break
        max_size=10000 #max columns allowed in tables
        count=0
        chunk_start=0
        chunck_end=0
        column_count = cascade_order.shape[1]
        if column_count > max_size:
            while column_count > max_size:
                count += 1
                # print(f'{column_count} columns > {max_size}, cascade table break {count}')
                chunck_end += max_size
                df = cascade_order.iloc[:, chunk_start:chunck_end]
                df.to_json(f'{sample_path_name}/cascade_order{count}.json', orient='split')
                self.excel_formatter(f'{sample_path_name}/cascade_order{count}.json', f'{sample_path_name}/cascade_table{count}-{st}.xlsx')
                os.remove(f'{sample_path_name}/cascade_order{count}.json')
                chunk_start += max_size
                column_count -= max_size
            count += 1
            # print(f'Last break {column_count} columns, cascade table break {count}')
            df = cascade_order.iloc[:, chunk_start:]
            df.to_json(f'{sample_path_name}/cascade_order{count}.json', orient='split')
            self.excel_formatter(f'{sample_path_name}/cascade_order{count}.json', f'{sample_path_name}/cascade_table{count}-{st}.xlsx')
            os.remove(f'{sample_path_name}/cascade_order{count}.json')
        else: # no break needed
            cascade_order.to_json(f'{sample_path_name}/cascade_order.json', orient='split')
            outfile = snp_alignment.checksum_match_to_text(f'{sample_path_name}/cascade_order.json')
            self.excel_formatter(outfile, f'{sample_path_name}/cascade_table-{st}.xlsx')
            os.remove(f'{sample_path_name}/cascade_order.json')

    def excel_formatter(self, df_json, write_to, group=None):
        import pandas.io.formats.excel
        pandas.io.formats.excel.header_style = None
        sample_path_name = self.sample_path_name
        st = self.st
        table_df = pd.read_json(df_json, orient='split')
        if self.gbk:
            table_df = self.annotate_table(table_df, sample_path_name, group, self.gbk)
        else:
            table_df = table_df.append(pd.Series(name='no annotations'))
        writer = pd.ExcelWriter(write_to, engine='xlsxwriter')
        table_df.to_excel(writer, sheet_name='Sheet1')
        wb = writer.book
        ws = writer.sheets['Sheet1']
        formatA = wb.add_format({'bg_color': '#58FA82'})
        formatG = wb.add_format({'bg_color': '#F7FE2E'})
        formatC = wb.add_format({'bg_color': '#0000FF'})
        formatT = wb.add_format({'bg_color': '#FF0000'})
        formatnormal = wb.add_format({'bg_color': '#FDFEFE'})
        formatlowqual = wb.add_format({'font_color': '#C70039', 'bg_color': '#E2CFDD'})
        formathighqual = wb.add_format({'font_color': '#000000', 'bg_color': '#FDFEFE'})
        formatambigous = wb.add_format({'font_color': '#C70039', 'bg_color': '#E2CFDD'})
        formatN = wb.add_format({'bg_color': '#E2CFDD'})
        rows, cols = table_df.shape

        ws.set_column(0, 0, 30)
        ws.set_column(1, cols, 2.1)
        ws.freeze_panes(2, 1)
        formatannotation = wb.add_format({'font_color': '#0A028C', 'rotation': '-90', 'align': 'top'})
        #set last row
        ws.set_row(rows + 1, cols + 1, formatannotation)

        #'first_row', 'first_col', 'last_row', and 'last_col'
        # Careful that row/column locations don't overlap
        ws.conditional_format(rows - 2, 1, rows - 1, cols, {'type': 'cell', 'criteria': '<', 'value': 55, 'format': formatlowqual})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'cell', 'criteria': '==', 'value': 'B$2', 'format': formatnormal})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'A', 'format': formatA})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'G', 'format': formatG})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'C', 'format': formatC})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'T', 'format': formatT})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'S', 'format': formatambigous})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'Y', 'format': formatambigous})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'R', 'format': formatambigous})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'W', 'format': formatambigous})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'K', 'format': formatambigous})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'M', 'format': formatambigous})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'N', 'format': formatN})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': '-', 'format': formatN})

        format_rotation = wb.add_format({})
        format_rotation.set_rotation(90)
        # ws.set_row(0, None, format_rotation)
        for columnnum, columnname in enumerate(list(table_df.columns)):
            ws.write(0, columnnum + 1, columnname, format_rotation)
        formatannotation = wb.add_format({'font_color': '#0A028C', 'rotation': '-90', 'align': 'top'})
        #set last row
        ws.set_row(rows, 400, formatannotation)
        writer.save()
	
    def annotate_table(self, table_df, sample_path_name, group, gbk_in):
        gbk_dict = SeqIO.to_dict(SeqIO.parse(gbk_in, "genbank"))
        annotation_dict = {}
        for chromosome in list(gbk_dict.keys()):
            write_out = open(f'{sample_path_name}/{group}/{group}_temp.csv', 'w+')
            for feature in gbk_dict[chromosome].features:
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
                    print(chromosome, int(feature.location.start), int(feature.location.end), mylocus, myproduct, mygene, sep='\t', file=write_out)
                    
            write_out.close()

            df = pd.read_csv(f'{sample_path_name}/{group}/{group}_temp.csv', sep='\t', names=["chrom", "start", "stop", "locus", "product", "gene"])
            os.remove(f'{sample_path_name}/{group}/{group}_temp.csv')
            df = df.sort_values(['start', 'gene'], ascending=[True, False])
            df = df.drop_duplicates('start')
            pro = df.reset_index(drop=True)
            pro.index = pd.IntervalIndex.from_arrays(pro['start'], pro['stop'], closed='both')
            annotation_dict[chromosome] = pro
        for gbk_chrome, pro in annotation_dict.items():
            ref_pos = list(table_df)
            ref_series = pd.Series(ref_pos)
            ref_df = pd.DataFrame(ref_series.str.split(':', expand=True).values, columns=['reference', 'position'])
            all_ref = ref_df[ref_df['reference'] == gbk_chrome]
            write_out = open(f'{sample_path_name}/{group}/{group}annotations.csv', 'a')
            positions = all_ref.position.to_frame()
            for index, row in positions.iterrows():
                pos = row.position
                try:
                    aaa = pro.iloc[pro.index.get_loc(int(pos))][['chrom', 'locus', 'product', 'gene']]
                    try:
                        chrom, name, locus, tag = aaa.values[0]
                        print("{}:{}\t{}, {}, {}".format(chrom, pos, locus, tag, name), file=write_out)
                    except ValueError:
                        # if only one annotation entire chromosome (such with flu) then having [0] fails
                        chrom, name, locus, tag = aaa.values
                        print("{}:{}\t{}, {}, {}".format(chrom, pos, locus, tag, name), file=write_out)
                except KeyError:
                    print("{}:{}\tNo annotated product" .format(gbk_chrome, pos), file=write_out)
            write_out.close()

        annotations_df = pd.read_csv(f'{sample_path_name}/{group}/{group}annotations.csv', sep='\t', header=None, names=['index', 'annotations'], index_col='index')
        os.remove(f'{sample_path_name}/{group}/{group}annotations.csv')
        table_df_transposed = table_df.T
        table_df_transposed.index = table_df_transposed.index.rename('index')
        table_df_transposed = table_df_transposed.merge(annotations_df, left_index=True, right_index=True)
        table_df = table_df_transposed.T
        return table_df

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
        ---------------------------------------------------------
        2019-12-12 Beta:
        Provide an aligned FASTA, build a corresponding parsimonious table and tree.
        Must have one FASTA header labeled as "reference_seq".  This will be used as the reference/root/outgroup. Should be the most likely index isolate/most close precursor to outbreak. 
        vsnp_fasta_to_snps_table.py -f *fasta
        ---------------------------------------------------------


    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--fasta', action='store', dest='fasta', default=None, required=True, help='provide an aligned fasta')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    startTime = datetime.now()
    snp_alignment = Get_Snps()

    with open('reformated.fasta', 'w') as reformat:
        sequence = SeqIO.parse(args.fasta, "fasta")
        for each in sequence:
            print(f'>{each.description}\n{each.seq}', file=reformat)

    uniq_id_fasta = snp_alignment.checksum_fastas_for_tree('reformated.fasta')

    df = pd.read_csv(uniq_id_fasta, header=None, sep='^')
    seq = df.iloc[1::2].reset_index(drop=True)
    header = df.iloc[0::2].reset_index(drop=True).replace(to_replace=r'>', value='', regex=True)
    seq = seq.rename({0:"seq"}, axis='columns')
    header = header.rename({0:"header"}, axis='columns')
    df = pd.concat([header, seq], axis='columns', ignore_index=True)
    df = df.rename({0: 'header', 1: 'seq'}, axis='columns')
    seq = df['seq'].apply(lambda x: pd.Series(list(x)))
    df = pd.concat([header, seq], axis='columns', ignore_index=True)
    df = df.set_index(0)
    prepared_from_fasta_json = f'{snp_alignment.sample_path_name}/prepared_from_infasta.json'
    df.to_json(prepared_from_fasta_json)
    parsimonious_json, tree_file = snp_alignment.gather_and_filter(prepared_from_fasta_json)
    print("Building tables...")
    snp_alignment.build_tables(parsimonious_json, tree_file)

    outfile = snp_alignment.checksum_match_to_text(tree_file) ### CALL FUNCTION

    runtime = datetime.now() - startTime
    print(f'\n\n Total Runtime: {runtime}:  \n')