#!/usr/bin/env python

__version__ = "2.03"

import os
import sys
import glob
import vcf
import numpy as np
import pandas as pd
from pathlib import Path
import argparse
import textwrap
from concurrent import futures
import multiprocessing
multiprocessing.set_start_method('spawn', True)
from concurrent import futures


class Filter_Finder:
    '''

    '''

    def __init__(self, working_directory, defining_snp_filter_file, group, debug=False):

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

        if group:
            self.group = group
        else:
            self.group = self.cwd.split('/')[-1] #this "group" will be used to match to columns in defining_snp_filter_file
        self.cpu_count_half = int(multiprocessing.cpu_count() / 2)
        self.vcf_list = glob.glob(f'{working_directory}/*vcf')
        self.defining_snp_filter_file = defining_snp_filter_file
        self.debug = debug

    def find_filter_dict(self, each_vcf):
        dict_qual = {}
        dict_map = {}
        try:
            vcf_reader = vcf.Reader(open(each_vcf, 'r'))
            for record in vcf_reader:
                try:
                    # Freebayes VCFs place MQ values are placed into a list.  GATK as a float
                    record.INFO['MQ'] = record.INFO['MQ'][0]
                except TypeError:
                    pass
                except KeyError:
                    pass
                absolute_positon = str(record.CHROM) + ":" + str(record.POS)
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
        except ValueError:
            print(f'Fix VCF: {each_vcf} Try running: from vSNP_step2 import Get_Snps, snp_alignment = Get_Snps(), snp_alignment.fix_vcfs()  Then rerun.')
        except SyntaxError:
            print(f'VCF file is corrupt and therefore not included in the analysis: {each_vcf}')
            return {}, {}


    def filter_finder(self, ):
        #write to files
        vcf_list = self.vcf_list
        find_filter_dict = self.find_filter_dict
        positions_to_filter = f'{self.group}_positions_to_filter.txt'
        positions_to_filter_details = f'{self.group}_positions_to_filter_details.txt'
        good_snps = f'{self.group}_good_snps_details.txt'
        write_out_positions = open(f'{self.cwd}/{positions_to_filter}', 'w')
        write_out_details = open(f'{self.cwd}/{positions_to_filter_details}', 'w')
        write_out_good_snps = open(f'{self.cwd}/{good_snps}', 'w')

        #calculate mean/max qual and map at all possible positions
        dd_qual = {}
        dd_map = {}
        if self.debug:
            for each_vcf in vcf_list:
                # print("working on: %s" % each_vcf)
                dict_qual, dict_map = find_filter_dict(each_vcf)
                keys = set(dd_qual).union(dict_qual)
                no = []
                #make position (key) and qual/maps list (value)
                dd_qual = dict((k, dd_qual.get(k, no) + dict_qual.get(k, no)) for k in keys)
                keys = set(dd_map).union(dict_map)
                no = []
                dd_map = dict((k, dd_map.get(k, no) + dict_map.get(k, no)) for k in keys)
        else:
            with futures.ProcessPoolExecutor(max_workers=self.cpu_count_half) as pool: #ProcessPoolExecutor ThreadPoolExecutor
                for dict_qual, dict_map in pool.map(find_filter_dict, vcf_list, chunksize=8):
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

        def print_possible_position_to_filter(all_maybe_filter):
            for absolute_positon in all_maybe_filter:
                ave_qual_value = ave_qual[absolute_positon]
                max_qual_value = max_qual[absolute_positon]
                ave_map_value = ave_map[absolute_positon]
                max_map_value = max_map[absolute_positon]
                # print("%s, max_qual_value: %s, ave_qual_value: %s, max_map_value: %s, ave_map_value: %s" % (absolute_positon, max_qual_value, ave_qual_value, max_map_value, ave_map_value))
                if max_qual_value < 1300 and ave_qual_value < 700 or ave_map_value < 56:
                    print("%s, max_qual_value: %s, ave_qual_value: %s, max_map_value: %s, ave_map_value: %s" % (absolute_positon, max_qual_value, ave_qual_value, max_map_value, ave_map_value), file=write_out_details)
                    print(absolute_positon, file=write_out_positions)
                else:
                    print("%s, max_qual_value: %s, ave_qual_value: %s, max_map_value: %s, ave_map_value: %s" % (absolute_positon, max_qual_value, ave_qual_value, max_map_value, ave_map_value), file=write_out_good_snps)
            write_out_positions.close()
            write_out_details.close()
            write_out_good_snps.close()

        # Removing those already from all positions to filter
        if self.defining_snp_filter_file:
            defining_snp_filter_file = self.defining_snp_filter_file
            exclusion_list=[]
            try:
                xl_filters = pd.read_excel(defining_snp_filter_file, header=1, usecols=[self.group])
                for value in xl_filters.values:
                    value = str(value[0])
                    if "-" not in value.split(":")[-1]:
                        exclusion_list.append(value)
                    elif "-" in value:
                        try:
                            chrom, sequence_range = value.split(":")
                        except ValueError as e:
                            raise type(e)(str(e) + f' \n#### error in {defining_snp_filter_file}\n#### see value "{value}"').with_traceback(sys.exc_info()[2])
                        value = sequence_range.split("-")
                        for position in range(int(value[0].replace(',', '')), int(value[1].replace(',', '')) + 1):
                            exclusion_list.append(chrom + ":" + str(position))
                len(exclusion_list)  
                all_maybe_filter_minus_already_in_excel = [x for x in all_maybe_filter if x not in exclusion_list] #just keep those not already in the excel filter file
                print_possible_position_to_filter(all_maybe_filter_minus_already_in_excel)
            except (KeyError, ValueError) as e:  # directory or group specified does not match any column names in defining_snp/filter file
                print(f'Error: {e}')
                print_possible_position_to_filter(all_maybe_filter)


if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-w', '--cwd', action='store', dest='working_directory', required=False, default='.', help='Working directoy to be added as path to dependency files.')
    parser.add_argument('-f', '--filters', action='store', dest='defining_snp_filter_file', required=False, default=None, help='Defining SNPs/Filter file')
    parser.add_argument('-g', '--group', action='store', dest='group', required=False, default=None, help='Group/directory name')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='turn off map.pooling of samples')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    filter_finder = Filter_Finder(args.working_directory, args.defining_snp_filter_file, args.group, args.debug)
    filter_finder.filter_finder()