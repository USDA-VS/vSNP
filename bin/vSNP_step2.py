#!/usr/bin/env python

__version__ = "0.2.01"

import os
import shutil
import subprocess
import sys
import re
import glob
import time
import json
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
from Bio import Phylo
from cpuinfo import get_cpu_info
import pylab

from vsnp_reference_options import Ref_Options
from vsnp_file_management import File_Management
from vsnp_chromosome_reference import Reference_Chromosome
from vsnp_filter_finder import Filter_Finder
from vsnp_remove_from_analysis import Remove_From_Analysis
from vsnp_html_step2_summary import Step2_Summary


class Get_Snps:
    ''' 
    Files must have .vcf extension
    Called or given a working directory collect quality parsimonious SNPs from VCF files and output alignment file in FASTA format.
    '''

    '''
    Malformed VCF rejects tracker summary:
    init()
        self.all_mal = []
    fix_vcfs()
        fix_each_vcf() --> pool
            mal=[] --c return
    group_vcfs()
        self.all_mal.append()
    for each_directory: --> pool
        get_snps() --> return
            find_initial_position() --> return
    '''

    def __init__(self, reference=None, working_directory='.', gbk=None, filter_finder=False, debug=False, no_filters=False, all_vcf=False, table=False, qual_threshold=150, MQ=56, AC=2, N_threshold=50):
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
        self.filter_finder = filter_finder
        self.no_filters = no_filters
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

        if reference:
            reference_options = Ref_Options(reference)

        if table:
            #print and exit
            reference_options.print_options()
            sys.exit(0)
        try:
            if reference_options.gbk:
                self.gbk = reference_options.gbk
        except UnboundLocalError:
            print(f'reference_options not set - message 1')

        script_path = os.path.dirname(os.path.realpath(__file__))

        if working_directory == '.':
            working_directory = os.getcwd()
        else:
            working_directory = working_directory
        
        try:
            self.excel_path = reference_options.excel
        except AttributeError:
            print(f'{reference} reference type not available, see available options:')
            reference_options.print_options()
        except UnboundLocalError:
            print(f'reference_options not set - message 2')

        if not reference:
            ref_options = glob.glob(f'{script_path}/defining_snps_filters/*.xlsx')
            print(f'\nUse -t option for valid reference options.')
            for option in ref_options:
                print(f'\t{os.path.basename(option).replace(".xlsx", "")}')
            # print('\n### Whoa!  Are you sure you want to run without specifying a reference, "-r"\n')
            # input("--> Press Enter to continue with an all_vcf table and tree\n--> ctrl-c to exit\n")
            self.all_vcf = True
            self.excel_path = None
            print("Continuing...")
        
        '''
        Find an optimal compiled version of RAxML in conda
        '''
        cpu_info_dict = get_cpu_info()
        flags_list = cpu_info_dict["flags"]
        try:
            subprocess.call("raxml", stdout=open(os.devnull, 'wb'))
            raxml = 'raxml'
        except OSError:
            if 'avx2' in flags_list:
                print(f'AVX2 available')
                raxml = 'raxmlHPC-PTHREADS-AVX2'
            elif 'sse3' in flags_list:
                print(f'SSE3 available')
                raxml = 'raxmlHPC-PTHREADS-SSE3'
            else:
                print(f'Neither SSE3 or AVX2 are available')
                raxml = 'raxmlHPC'
            print(f'set RAxML to {raxml}')
        self.raxml = raxml

    def zipit(self, src, dst):
        zf = zipfile.ZipFile("%s.zip" % (dst), "w", zipfile.ZIP_DEFLATED)
        abs_src = os.path.abspath(src)
        for dirname, subdirs, files in os.walk(src):
            for filename in files:
                absname = os.path.abspath(os.path.join(dirname, filename))
                arcname = absname[len(abs_src) + 1:]
                zf.write(absname, arcname)
        zf.close()
        shutil.rmtree(src)

    def find_initial_positions(self, filename):
        found_positions = {}
        found_positions_mix = {}
        AC = self.AC
        qual_threshold = self.qual_threshold
        MQ = self.MQ
        mal = []
        try:
            vcf_reader = vcf.Reader(open(filename, 'r'))
            try:
                for record in vcf_reader:
                    try:
                        record_qual = int(record.QUAL)
                    except TypeError:
                        record_qual = 0 
                    try:
                        # Freebayes VCFs place MQ values are placed into a list.  GATK as a float
                        record.INFO['MQ'] = record.INFO['MQ'][0]
                    except TypeError:
                        pass
                    except KeyError:
                        pass
                    chrom = record.CHROM
                    position = record.POS
                    absolute_positon = str(chrom) + ":" + str(position)
                    try:
                        if str(record.ALT[0]) != "None" and record.INFO['AC'][0] == AC and len(record.REF) == 1 and record_qual > qual_threshold and record.INFO['MQ'] > MQ:
                            found_positions.update({absolute_positon: record.REF})
                        if str(record.ALT[0]) != "None" and record.INFO['AC'][0] == 1 and len(record.REF) == 1 and record_qual > qual_threshold and record.INFO['MQ'] > MQ:
                            found_positions_mix.update({absolute_positon: record.REF})
                    except KeyError as e:
                        # raise type(e)(str(e) + f' \n#### error in {filename}\n#### see value "{absolute_positon}"').with_traceback(sys.exc_info()[2])
                        # all_mal.append(type(e)(str(e) + f' \n#### error in {filename}\n#### see value "{absolute_positon}"').with_traceback(sys.exc_info()[2]))
                        os.remove(filename)
                        mal.append(os.path.basename(filename))
                        pass
                return filename, mal, found_positions, found_positions_mix
            except (ZeroDivisionError, ValueError, UnboundLocalError, TypeError) as e:
                # raise type(e)(str(e) + f' \n#### error in {filename}\n#### see value "{record.CHROM}:{record.POS}"').with_traceback(sys.exc_info()[2])
                # all_mal.append(type(e)(str(e) + f' \n#### error in {filename}\n#### see value "{record.CHROM}:{record.POS}"').with_traceback(sys.exc_info()[2]))
                os.remove(filename)
                mal.append(os.path.basename(filename))
                return filename, mal, {'': ''}, {'': ''}
        except (SyntaxError, AttributeError) as e:
            # print(type(e)(str(e) + f'\n### VCF SyntaxError {filename} File Removed'))
            # all_mal.append(f'VCF SyntaxError {filename} File Removed')
            os.remove(filename)
            mal.append(os.path.basename(filename))
            return filename, mal, {'': ''}, {'': ''}
        
    def tree_to_svg(self, tree_file, svg_file):

        def get_label(leaf):
            return leaf.name
        tree = Phylo.read(tree_file, 'newick')
        tree.ladderize()
        Phylo.draw(tree, label_func=get_label, do_show=False)
        pylab.axis('off')
        pylab.savefig(svg_file, format='svg', bbox_inches='tight', dpi=500)

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

    def decide_snps(self, filename):
        all_positions = self.all_positions
        sample_row=[]
        sample_map_qualities = {}
        mal = []
        just_name = os.path.basename(filename)
        just_name = just_name.replace('.vcf', '')
        just_name = re.sub('\..*', '*', just_name) # if after the .vcf is removed there is stilll a "." in the name it is assumed the name did not get changed
        sample_row.append(just_name)
        # for each line in vcf
        vcf_reader = vcf.Reader(open(filename, 'r'))
        sample_dict = {}
        try:
            for record in vcf_reader:
                try:
                    # Freebayes VCFs place MQ values are placed into a list.  GATK as a float
                    record.INFO['MQ'] = record.INFO['MQ'][0]
                except TypeError:
                    pass
                except KeyError:
                    pass
                record_position = str(record.CHROM) + ":" + str(record.POS)
                try:
                    record_qual = int(record.QUAL)
                except TypeError:
                    record_qual = 0 
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
                    if str(record.ALT[0]) != "None" and len(record.ALT[0]) == 1 and record.INFO['AC'][0] == 2 and record_qual > self.N_threshold:
                        sample_dict.update({record_position: str(record.ALT[0])})
                    elif str(record.ALT[0]) != "None" and len(record.ALT[0]) == 1 and record.INFO['AC'][0] == 1 and record_qual > self.N_threshold:
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
                    elif str(record.ALT[0]) != "None" and record_qual <= 50:
                        sample_dict.update({record_position: record.REF[0]})
                    elif str(record.ALT[0]) != "None" and record_qual <= self.N_threshold:
                        sample_dict.update({record_position: "N"})
                    elif str(record.ALT[0]) != "None": #Insurance -- Will still report on a possible SNP even if missed with above statement
                        sample_dict.update({record_position: str(record.REF[0])})
                    elif str(record.ALT[0]) == "None":
                        sample_dict.update({record_position: "-"})

            # merge dictionaries and order
            merge_dict={}
            merge_dict.update(all_positions) #abs_pos:REF
            merge_dict.update(sample_dict) # abs_pos:ALT replacing all_positions, because keys must be unique
            sample_df = pd.DataFrame(merge_dict, index=[just_name])
        except ValueError:
            print(f'### FILE REMOVED: ValueError: {filename}')
            os.remove(filename)
            mal.append(os.path.basename(filename))
            sample_df = None
        return sample_df, just_name, sample_map_qualities, mal

    def get_snps(self, directory):

        if self.filter_finder:
            print(f'Finding filters: {directory}')
            filter_finder = Filter_Finder(f'{self.cwd}/{directory}', self.excel_path, directory, self.debug)
            filter_finder.filter_finder()

        all_positions = {}
        df_list = []
        collect_mal_list = []
        sample_path_name = f'{self.sample_path_name}/{directory}'
        vcf_list = glob.glob(f'{sample_path_name}/*.vcf') #reset because some could have been deleted, removing from self.list in pool function doesn't work
        # print(f'\nIn {directory} finding SNP position in {len(vcf_list):,} sample...')
        find_initial_positions = self.find_initial_positions
        for filename in vcf_list:
            try:
                filename, mal, found_positions, found_positions_mix = find_initial_positions(filename)
                all_positions.update(found_positions)
            except TypeError:
                mal = [f'TypeError for {filename}']
                #if "found_positions" is not a dictionary it is error and should not be added.  See log for malformed VCF files
                pass
            collect_mal_list = collect_mal_list + mal
        # print(f'{len(all_positions):,} SNP positions found')
        # order before adding to file to match with ordering of individual samples below
        # all_positions is abs_pos:REF
        self.all_positions = OrderedDict(sorted(all_positions.items()))
        ref_positions_df = pd.DataFrame(self.all_positions, index=['root'])
        vcf_list = glob.glob(f'{sample_path_name}/*.vcf') #reset because some could have been deleted
        all_map_qualities = {}
        # print(f'Deciding which SNPs to keep from {len(vcf_list):,} samples...')
        decide_snps = self.decide_snps
        for filename in vcf_list:
            sample_df, just_name, sample_map_qualities, mal = decide_snps(filename)
            df_list.append(sample_df)
            all_map_qualities.update({just_name: sample_map_qualities})
            collect_mal_list = collect_mal_list + mal
        all_sample_df = pd.concat(df_list)
        # print(" Done")
        ##### ALL POSITIONS HAVE NOW BEEN SELECTED FOR EACH SAMPLE #####
        # SELECT PARISOMONY INFORMATIVE SNPSs - removes columns where all fields are the same
        # pars_df = self.get_parsimonious_pos(all_sample_df) #not necessary for all unless few samples and -a option
        prefilter_df = pd.concat([ref_positions_df, all_sample_df], join='inner') #add reference to top row
        all_mq_df = pd.DataFrame.from_dict(all_map_qualities)
        mq_averages = all_mq_df.mean(axis=1).astype(int)
        mq_averages.to_json(f'{sample_path_name}/average_mq.json', orient='split')
        ## For debugging
        prefilter_df.to_json(f'{sample_path_name}/prefilter_df.json')

        # prefilter_df.to_csv('/Users/tstuber/Desktop/script_test_files/suis1_step2/test/text.txt', sep='\t')
        self.gather_and_filter(alignment=prefilter_df, excel_path=self.excel_path, directory=directory)
        if os.path.isfile(f'{self.sample_path_name}/{directory}/{directory}-{self.st}.tre'):
            self.build_tables(directory)

        ### Remove unwanted files
        vcf_list = glob.glob(f'{sample_path_name}/*.vcf')
        for each_file in vcf_list:
            os.remove(each_file)
        json_list = glob.glob(f'{sample_path_name}/*.json')
        for each_file in json_list:
            os.remove(each_file)
        return collect_mal_list

    def get_position_list(self, sheet_names, excel_path, usecols):
        exclusion_list=[]
        try:
            filter_to_all = pd.read_excel(excel_path, header=1, usecols=[usecols])
            for value in filter_to_all.values:
                value = str(value[0])
                if "-" not in value.split(":")[-1]:
                    exclusion_list.append(value)
                elif "-" in value:
                    try:
                        chrom, sequence_range = value.split(":")
                    except ValueError as e:
                        raise type(e)(str(e) + f' \n#### error in {excel_path}\n#### see value "{value}"').with_traceback(sys.exc_info()[2])
                    value = sequence_range.split("-")
                    for position in range(int(value[0].replace(',', '')), int(value[1].replace(',', '')) + 1):
                        exclusion_list.append(chrom + ":" + str(position))
            return exclusion_list
        except ValueError:
            exclusion_list = []
            return exclusion_list

    def get_parsimonious_pos(self, in_df):
        try:
            ref_series = in_df.loc['root']
            in_df = in_df.drop(['root']) #in all_vcf root needs to be removed
        except KeyError:
            pass
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
        else:
            sample_path_name = f'{self.sample_path_name}/{directory}'
        group = directory
        raxml = self.raxml
        st = self.st
        if excel_path is None or self.no_filters:
            filtered_all_df = alignment
            ref_series = filtered_all_df.loc['root']
            sheet_names = None
        elif excel_path:
            #filter positions to be removed from all
            xl = pd.ExcelFile(excel_path)
            sheet_names = xl.sheet_names
            exclusion_list_all = self.get_position_list(sheet_names, excel_path, 0) #Use the first column to filter "all" postions
            exclusion_list_group = self.get_position_list(sheet_names, excel_path, directory) #Use the first column to filter "all" postions
            exclusion_list = exclusion_list_all + exclusion_list_group
            filtered_all_df = alignment.drop(columns=exclusion_list, errors='ignore')  #filters for all applied
            # ref_series = filtered_all_df.loc['root']
            # print(f'{filtered_all_df.shape} Table size after filtering')
        parsimonious_df = self.get_parsimonious_pos(filtered_all_df)
        parsimonious_df.to_json(f'{sample_path_name}/{group}_filtered.json', orient='split')
        alignment_file = f'{sample_path_name}/{group}_alignment-{st}.fasta'
        self.df_to_fasta(parsimonious_df, alignment_file)
        samples_number, columns = parsimonious_df.shape
        if samples_number < 4:
            with open(f'{sample_path_name}/TOO_FEW_SAMPLES_TO_BUILD_TREE', 'w') as message_out:
                print(f'check sample numbers', file=message_out)
        else:
            os.system(f'{raxml} -s {alignment_file} -n raxml -m GTRCATI -o root -w {sample_path_name} -p 456123 -T 4 > /dev/null 2>&1')
            try:
                os.rename(f'{sample_path_name}/RAxML_bestTree.raxml', f'{sample_path_name}/{group}-{st}.tre')
                # self.tree_to_svg(f'{sample_path_name}/{group}-{st}.tre', f'{sample_path_name}/{group}-{st}.svg')
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

    def group_vcfs(self, excel_path):
        sample_path_name = self.sample_path_name
        find_initial_positions = self.find_initial_positions
        group_names=[]
        xl = pd.ExcelFile(excel_path)
        sheet_names = xl.sheet_names
        ws = pd.read_excel(excel_path, sheet_name=sheet_names[0])
        defining_snps = ws.iloc[0]
        defsnp_iterator = iter(defining_snps.iteritems())
        next(defsnp_iterator)
        defining_snps={}
        inverted_defining_snps={}
        for abs_pos, group in defsnp_iterator:
            group_names.append(group)
            if '!' in abs_pos:
                inverted_defining_snps[abs_pos.replace('!', '')] = group
            else:
                defining_snps[abs_pos] = group #Make defining snp/group dict
        vcf_list = glob.glob(f'{sample_path_name}/*.vcf')
        pre_vcf_list_size = len(vcf_list)
        for filename in vcf_list:
            if os.stat(filename).st_size == 0:
                print(f'\n#### Error, File empty {filename}, removed')
                os.remove(filename)
                self.all_mal.append(os.path.basename(filename))
        vcf_list = glob.glob(f'{sample_path_name}/*.vcf')
        print(f'{len(vcf_list):,} VCF files')
        count = 0
        samples_groups_dict = {}
        gather_mal_list = []
        def bin_and_html_table(filename, mal, found_positions, found_positions_mix, count):
            sample_groups_list = []
            gather_mal_list = []
            tablename = os.path.basename(filename)
            if mal:
                samples_groups_dict[tablename] = ['<font color="red">Rejected VCF file</font>']
                gather_mal_list = gather_mal_list + mal
            else:
                try:
                    defining_snp = False
                    for abs_position in list(defining_snps.keys() & (found_positions.keys() | found_positions_mix.keys())): #absolute positions in set union of two list
                        # print(f'{defining_snps[abs_position]}: {filename}')
                        sys.stdout.write("\r%i VCF files moved" % count)
                        sys.stdout.flush()
                        group = defining_snps[abs_position]
                        group_directory = f'{sample_path_name}/{group}'
                        sample_groups_list.append(group)
                        if len(list(defining_snps.keys() & found_positions_mix.keys())) > 0:
                            tablename = f'{os.path.basename(filename)} <font color="red">[[MIXED]]</font>'
                        if not os.path.exists(group_directory):
                            os.makedirs(group_directory)
                        shutil.copy(filename, group_directory)
                        defining_snp = True
                    if not set(inverted_defining_snps.keys()).intersection(found_positions.keys() | found_positions_mix.keys()):
                        for abs_position in list(inverted_defining_snps.keys()):
                            sys.stdout.write("\r%i VCF files moved" % count)
                            sys.stdout.flush()
                            group = inverted_defining_snps[abs_position]
                            group_directory = f'{sample_path_name}/{group}'
                            sample_groups_list.append(group)
                            if not os.path.exists(group_directory):
                                os.makedirs(group_directory)
                            shutil.copy(filename, group_directory)
                            defining_snp = True
                    if defining_snp:
                        samples_groups_dict[tablename] = sorted(sample_groups_list)
                    else:
                        samples_groups_dict[tablename] = ['<font color="red">No defining SNP</font>']
                except TypeError:
                    message = f'File TypeError'
                    print(f'{message}: {filename}')
                    samples_groups_dict[tablename] = [f'{message}: {filename}']
                    pass
            return samples_groups_dict, gather_mal_list
        if self.debug:
            for filename in vcf_list:
                filename, mal, found_positions, found_positions_mix = find_initial_positions(filename)
                count += 1
                samples_groups_dict, gather_mal_list = bin_and_html_table(filename, mal, found_positions, found_positions_mix, count)
        else:
            with futures.ProcessPoolExecutor(max_workers=self.cpu_count_half) as pool: #ProcessPoolExecutor ThreadPoolExecutor
                for filename, mal, found_positions, found_positions_mix in pool.map(find_initial_positions, vcf_list):
                    count += 1
                    samples_groups_dict, gather_mal_list = bin_and_html_table(filename, mal, found_positions, found_positions_mix, count)
        self.all_mal = self.all_mal + gather_mal_list
        post_vcf_list_size = len(glob.glob(f'{sample_path_name}/*.vcf'))
        difference = pre_vcf_list_size - post_vcf_list_size
        print(f'\n{difference} corrupt and removed')
        self.htmlsummary['fixed_vcf_number_removed'] = self.htmlsummary['fixed_vcf_number_removed'] + difference
        return samples_groups_dict

    def build_tables(self, group):
        sample_path_name = self.sample_path_name
        st = self.st
        alignment_json = f'{sample_path_name}/{group}/{group}_filtered.json'
        mq = f'{sample_path_name}/{group}/average_mq.json'
        tree = f'{sample_path_name}/{group}/{group}-{st}.tre'
        alignment = pd.read_json(alignment_json, orient='split')
        os.remove(alignment_json)
        mq_series = pd.read_json(mq, typ='series', orient='split')
        with open(tree, 'rt') as tree_file: #must be the single line newick format.  Not Nexus which will be mutliline often with formating
            for line in tree_file:
                line = re.sub('[:,]', '\n', line)
                line = re.sub('[)(]', '', line)
                line = re.sub('[0-9].*\.[0-9].*\n', '', line)
                line = re.sub('root\n', '', line)
        sample_order = line.split('\n')
        sample_order = list(filter(None, sample_order))
        sample_order.insert(0, 'root')
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
        #row1 = row1.drop('root')

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
        #row2 = row2.drop('root')
        tree_order = tree_order.append([row1])
        tree_order = tree_order.append([row2])
        #In pandas=0.18.1 even this does not work:
        #    abc = row1.to_frame()
        #    abc = abc.T --> tree_order.shape (5, 18), abc.shape (1, 18)
        #    tree_order.append(abc)
        #Continue to get error: "*** ValueError: all the input arrays must have same number of dimensions"

        tree_order = tree_order.T
        tree_order = tree_order.sort_values(['snp_from_top', 'snp_per_column'], ascending=[True, False])
        tree_order = tree_order.T

        # map quality to dataframe
        mqdf = mq_series.to_frame(name='MQ')
        mqdf = mqdf.T

        # remove snp_per_column and snp_from_top rows
        cascade_order = tree_order[:-2]
        # cascading table
        cascade_order_mq = pd.concat([cascade_order, mqdf], join='inner')   
        ###Break
        max_size=10000 #max columns allowed in tables
        count=0
        chunk_start=0
        chunck_end=0
        column_count = cascade_order_mq.shape[1]
        if column_count > max_size:
            while column_count > max_size:
                count += 1
                # print(f'{column_count} columns > {max_size}, cascade table break {count}')
                chunck_end += max_size
                df = cascade_order_mq.iloc[:, chunk_start:chunck_end]
                df.to_json(f'{sample_path_name}/{group}/{group}_cascade_order_mq{count}.json', orient='split')
                self.excel_formatter(f'{sample_path_name}/{group}/{group}_cascade_order_mq{count}.json', f'{sample_path_name}/{group}/{group}_cascade_table{count}-{st}.xlsx', group)
                os.remove(f'{sample_path_name}/{group}/{group}_cascade_order_mq{count}.json')
                chunk_start += max_size
                column_count -= max_size
            count += 1
            # print(f'Last break {column_count} columns, cascade table break {count}')
            df = cascade_order_mq.iloc[:, chunk_start:]
            df.to_json(f'{sample_path_name}/{group}/{group}_cascade_order_mq{count}.json', orient='split')
            self.excel_formatter(f'{sample_path_name}/{group}/{group}_cascade_order_mq{count}.json', f'{sample_path_name}/{group}/{group}_cascade_table{count}-{st}.xlsx', group)
            os.remove(f'{sample_path_name}/{group}/{group}_cascade_order_mq{count}.json')
        else: # no break needed
            cascade_order_mq.to_json(f'{sample_path_name}/{group}/{group}_cascade_order_mq.json', orient='split')
            self.excel_formatter(f'{sample_path_name}/{group}/{group}_cascade_order_mq.json', f'{sample_path_name}/{group}/{group}_cascade_table-{st}.xlsx', group)
            os.remove(f'{sample_path_name}/{group}/{group}_cascade_order_mq.json')
        # sorted position table
        sort_df = cascade_order.T
        sort_df['abs_value'] = sort_df.index
        sort_df[['chrom','pos']] = sort_df['abs_value'].str.split(':', expand=True)
        sort_df = sort_df.drop(['abs_value', 'chrom'], axis=1)
        sort_df.pos = sort_df.pos.astype(int)
        sort_df = sort_df.sort_values(by=['pos'])
        sort_df = sort_df.drop(['pos'], axis=1)
        sort_df = sort_df.T
        sort_order_mq = pd.concat([sort_df, mqdf], join='inner')
        ###Break
        count=0
        chunk_start=0
        chunck_end=0
        column_count = sort_order_mq.shape[1]
        if column_count > max_size:
            while column_count > max_size:
                count += 1
                # print(f'{column_count} columns > {max_size}, sorted table break {count}')
                chunck_end += max_size
                df = sort_order_mq.iloc[:, chunk_start:chunck_end]
                df.to_json(f'{sample_path_name}/{group}/{group}_sort_order_mq_{count}.json', orient='split')
                self.excel_formatter(f'{sample_path_name}/{group}/{group}_sort_order_mq_{count}.json', f'{sample_path_name}/{group}/{group}_sort_table_{count}-{st}.xlsx', group)
                os.remove(f'{sample_path_name}/{group}/{group}_sort_order_mq_{count}.json')
                chunk_start += max_size
                column_count -= max_size
            count += 1
            # print(f'Last break {column_count} columns, sorted table break {count}')
            df = sort_order_mq.iloc[:, chunk_start:]
            df.to_json(f'{sample_path_name}/{group}/{group}_sort_order_mq_{count}.json', orient='split')
            self.excel_formatter(f'{sample_path_name}/{group}/{group}_sort_order_mq_{count}.json', f'{sample_path_name}/{group}/{group}_sort_table_{count}-{st}.xlsx', group)
            os.remove(f'{sample_path_name}/{group}/{group}_sort_order_mq_{count}.json')
        else: # no break needed
            sort_order_mq.to_json(f'{sample_path_name}/{group}/{group}_sort_order_mq.json', orient='split')
            self.excel_formatter(f'{sample_path_name}/{group}/{group}_sort_order_mq.json', f'{sample_path_name}/{group}/{group}_sort_table-{st}.xlsx', group)
            os.remove(f'{sample_path_name}/{group}/{group}_sort_order_mq.json')

    def excel_formatter(self, df_json, write_to, group):
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
    Current working directory used by default, but can specify working directory with -w.
    Directory must contain VCF files with file extension ".vcf"

    Usage: 
        # See available reference options:
        vSNP_step2.py -t

        # Run with a specific reference option:
        vSNP_step2.py -r Brucella_suis1

        # Find reference from VCF chrom column:
        vSNP_step2.py

        # If VCF chrom does not cross-reference with reference an all_vcf is made:
        vSNP_step2.py

        # An all VCF table can also be created along with using a reference option:
        vSNP_step2.py -ar Brucella_suis1

    Dependencies:
    - Reference options are grouped and accessed via named directories.  New directories are added using, $ path_adder.py.  In vSNP's installed package dependency paths are stored in, "dependency_paths.txt".  Directory/reference options are shown using -t option.
        Seven files can be included:
            Excel: (see template_define_filter.xlsx) with defining SNPs and filter positions.   <Required for grouping>
            Excel: metadata.xlsx  3 column file: VCF file name, updated file name, representative (optional boolean).  File name must contain "meta" somewhere in its name.  <Optional>
            Excel: remove_from_analysis.xlsx 1 column file: removes files based on name minus .vcf extension.  File name must contain "remove" somewhere in its name.  <Optional>
            FASTA (.fasta): used by vSNP_step1.py as reference.  <Required, unless explicitely given with -r option.  See senario 1>
            GBK (.gbk): used to annotate VCF files and tables.  <Optional>
            GFF (.gff): used by IGV to show annotation.  <Optional>
            IGV file: .genome IGV file mapping FASTA and GFF.  <Optional>

        - "template.xlsx" is availabe in "dependencies" directory.
          Definition: absolute position - combination of reference and position.  This is the combination of VCF CHROM and POS, CHROM:POS. 
          There are 3 aspects to Excel file:
            1) top row: absolute position to define a group
            2) second row: group name.  First column must contain "-All" in the name, other naming restrictions: no special characters or spaces.  Dashes and underscores are allowed in group name.
            3) third row and below: positions to filter from table and trees.
                Must be in CHROM:POS format
                POS numbers are allowed to have commas: both 1000 and 1,000 are accepted
                POS number ranges can be be used.  Note: CHROM:POS number range conforms to IGV
                
    - vSNP comes installed with Mycobacterium_AF2122.  For additional reference options see: https://github.com/USDA-VS/vSNP_dependencies.git

    ---------------------------------------------------------


    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-r', '--reference', action='store', dest='reference', default=None, required=False, help='provide a valid reference, see -t output')
    parser.add_argument('-t', '--table', action='store_true', dest='table', required=False, help='see valid reference types available')
    parser.add_argument('-a', '--all', action='store_true', dest='all', required=False, help='create table with all isolates')
    parser.add_argument('-s', '--subset', action='store_true', dest='subset', required=False, help='create trees with a subset of sample that represent the whole')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='turn off map.pooling of samples')
    parser.add_argument('-n', '--no_filters', action='store_true', dest='no_filters', default=False, help='run without applying filters')
    parser.add_argument('-f', '--filter_finder', action='store_true', dest='filter_finder', default=False, help='write possible positions to filter to text file')
    parser.add_argument('-w', '--cwd', action='store', dest='working_directory', required=False, default='.', help='Optional: path to VCF files')
    parser.add_argument('-g', '--gbk', action='store', dest='gbk', required=False, default=None, help='Optional: provide full path to gbk file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')
    args = parser.parse_args()
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    working_directory=args.working_directory

    if args.table:
        reference_options = Ref_Options(None)
        reference_options.print_options()
        sys.exit(0)

    startTime = datetime.now()
    starting_vcf_number = len(glob.glob(f'{working_directory}/*.vcf'))
    if args.reference == None:
        reference_type = Reference_Chromosome(working_directory)
        reference = reference_type.get_reference()
        print(f'Reference type found: {reference}\n')
        file_management = File_Management(excel_genotype_codes=None, working_directory=working_directory)
    else:
        reference = args.reference
    try:
        reference_options = Ref_Options(reference)
        reference_options.files_in_directory()
        if reference_options.remove:
            remove_from_analysis = Remove_From_Analysis(working_directory=working_directory, excel_remove=reference_options.remove)
            remove_from_analysis.remove_files()
        if reference_options.metadata:
            file_management = File_Management(excel_genotype_codes=reference_options.metadata, working_directory=working_directory)
            #Remove samples not in subset
            if args.subset:
                file_management.subset()
            file_management.change_names()
            file_management.backup_vcfs()
        else:
            file_management = File_Management(excel_genotype_codes=None, working_directory=working_directory)
            if args.subset:
                file_management.subset()
            file_management.backup_vcfs()
    except AttributeError:
        print(f'Check your reference optin name against `vsnp_path_adder.py -s` list')
        pass
    print('Fixing VCF files...')
    vcf_list=glob.glob(f'{args.working_directory}/*.vcf')
    all_mal = file_management.fix_vcfs(vcf_list=vcf_list, debug=args.debug)
    new_vcf_list = glob.glob(f'{args.working_directory}/*.vcf')
    difference_between_org_and_fixed = len(vcf_list) - len(new_vcf_list)

    print(f'Starting {starting_vcf_number:,} VCF files')

    snp_alignment = Get_Snps(reference=reference, working_directory=working_directory, gbk=args.gbk, filter_finder=args.filter_finder, debug=args.debug, no_filters=args.no_filters, all_vcf=args.all, table=args.table)
    snp_alignment.htmlsummary["starttime"] = startTime
    snp_alignment.htmlsummary["reference"] = reference
    snp_alignment.htmlsummary['all_vcf_boolen'] = args.all
    snp_alignment.htmlsummary['subset_boolen'] = args.subset
    snp_alignment.htmlsummary['starting_vcf_number'] = starting_vcf_number
    snp_alignment.all_mal = all_mal + snp_alignment.all_mal
    snp_alignment.htmlsummary['fixed_vcf_number_removed'] = difference_between_org_and_fixed
    try:
        snp_alignment.htmlsummary["remove_from_analysis_count"] = remove_from_analysis.removed_file_count
    except NameError:
        pass
    snp_alignment.htmlsummary['keep_count'] = file_management.keep_count    

    if snp_alignment.excel_path:
        samples_groups_dict = snp_alignment.group_vcfs(snp_alignment.excel_path)
        snp_alignment.htmlsummary["samples_groups_dict"] = samples_groups_dict
    if args.all or args.subset or snp_alignment.excel_path is None:
        os.makedirs(f'{working_directory}/all_vcf')
        for each in glob.glob(f'{working_directory}/*.vcf'):
            shutil.copy(each, f'{working_directory}/all_vcf')
    directory_list = [dd for dd in os.listdir(working_directory) if os.path.isdir(f'{working_directory}/{dd}')]

    count = 0
    print(f'\n\n{len(directory_list):,} groups')
    if snp_alignment.debug:
        for directory in directory_list:
            print(f'\n{directory}')
            collect_mal_list = snp_alignment.get_snps(directory)
            snp_alignment.all_mal = snp_alignment.all_mal + collect_mal_list
            count += 1
            sys.stdout.write("\r%i groups completed" % count)
            sys.stdout.flush()
    else:
        with futures.ProcessPoolExecutor(max_workers=snp_alignment.cpu_count_half) as pool: #ProcessPoolExecutor ThreadPoolExecutor
            for collect_mal_list in pool.map(snp_alignment.get_snps, directory_list):
                snp_alignment.all_mal = snp_alignment.all_mal + collect_mal_list
                count += 1
                sys.stdout.write("\r%i groups completed" % count)
                sys.stdout.flush()
                
    # os.remove(f'{snp_alignment.sample_path_name}/average_mq.json')
    snp_alignment.htmlsummary["malformed"] = snp_alignment.all_mal
    runtime = datetime.now() - startTime
    snp_alignment.htmlsummary["endtime"] = datetime.now()
    snp_alignment.htmlsummary["runtime"] = runtime
    html_step2_summary = Step2_Summary()
    html_step2_summary.html_step2_summary(snp_alignment.htmlsummary, working_directory=working_directory)

    try:
        if reference_options.excel:
            dependents_dir = f'{working_directory}/dependents'
            os.makedirs(dependents_dir)
            try:
                shutil.copy(reference_options.excel, dependents_dir)
                shutil.copy(reference_options.remove, dependents_dir)
            except (OSError, TypeError):
                pass
            snp_alignment.zipit(dependents_dir, dependents_dir)
    except AttributeError:
        with open(f'{snp_alignment.sample_path_name}/NO_DEPENDENCY_FILES', 'w') as message_out:
            print(f'ran without dependency files', file=message_out)

    for each_file in glob.glob(f'{snp_alignment.cwd}/*.vcf'):
        os.remove(each_file)

    print(f'\n\n Total Runtime: {runtime}:  \n')