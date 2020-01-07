#!/usr/bin/env python

__version__ = "2.0.0"

import os
import re
import glob
import shutil
import sys
import time
import zipfile
import pandas as pd
import argparse
import textwrap
from concurrent import futures
import multiprocessing
multiprocessing.set_start_method('spawn', True)

class File_Management:

    def __init__(self, excel_genotype_codes=None, working_directory='.', extension="vcf", debug=False):
        if excel_genotype_codes:
            df = pd.read_excel(excel_genotype_codes, index_col=0, usecols=[0, 1], names=['file_name', 'genotype_code'])
            excel_dict = df.to_dict('dict')['genotype_code']
            self.excel_dict = excel_dict
        self.excel_genotype_codes = excel_genotype_codes
        self.extension = extension
        self.change_name_file = f'{working_directory}/name_change_log.txt'
        self.working_directory = working_directory
        self.starting_files = f'{working_directory}/starting_files'
        self.keep_count = 0
        self.cpu_count_half = int(multiprocessing.cpu_count() / 2)
        self.debug = debug

    def get_genotype_code(self, sample):
        sample_list = [key for key in self.excel_dict.keys() if str(key) == str(sample)] #force key to string encase a int
        if len(sample_list) == 1:
            try:
                genotype_name = self.excel_dict[sample]
                genotype_name = re.sub("/", "_", genotype_name)
                genotype_name = re.sub("\.", "_", genotype_name)
                genotype_name = re.sub("\*", "_", genotype_name)
                genotype_name = re.sub("\?", "_", genotype_name)
                genotype_name = re.sub("\(", "_", genotype_name)
                genotype_name = re.sub("\)", "_", genotype_name)
                genotype_name = re.sub("\[", "_", genotype_name)
                genotype_name = re.sub("\]", "_", genotype_name)
                genotype_name = re.sub(" ", "_", genotype_name)
                genotype_name = re.sub("{", "_", genotype_name)
                genotype_name = re.sub("}", "_", genotype_name)
                genotype_name = re.sub("-_", "_", genotype_name)
                genotype_name = re.sub("_-", "_", genotype_name)
                genotype_name = re.sub("--", "_", genotype_name)
                genotype_name = re.sub("_$", "", genotype_name)
                genotype_name = re.sub("-$", "", genotype_name)
                genotype_name = re.sub("\'", "", genotype_name)
                genotype_name = re.sub(",", "", genotype_name)
                found = True
                return (genotype_name, found)
            except KeyError:
                found = False
                genotype_name = f'{sample} name not found in Excel file\nSample name minus everything left of first occurring _ or . must match text in column 1 exactly'
                return (genotype_name, found)
            except TypeError:
                found = False
                genotype_name = f'{sample} Sample found but no corresponding genotype code'
                return (genotype_name, found)
        else:
            found = False
            genotype_name = 'Multiple names found for sample'
            return (genotype_name, found)

    def subset(self,):
        '''
        Uses 3rd column of Excel file used to rename isolates.  If *anything* character is in this 3rd column the row's sample will be use in the subset if available.
        '''
        print(f'*** Only using a subset of samples ***')
        print(f'Using 3rd column of {self.excel_genotype_codes} to subset samples')
        input_list = glob.glob(f'{self.working_directory}/*{self.extension}')
        xl_filters = pd.read_excel(self.excel_genotype_codes, header=0, index_col=0)
        keep_list = xl_filters.iloc[:,1:2].dropna()
        keep_list = keep_list.index.tolist()
        not_in_subset = 'not_in_subset'
        os.makedirs(not_in_subset)
        keep_count = self.keep_count
        for sample in input_list:
            time_test = time.time() - os.path.getmtime(sample) < (1 * 24 * 60 * 60) # 1day * (24*60*60)sec in day
            sample_base = os.path.basename(sample)
            sample_name = re.sub("[._].*", "", sample_base)
            if sample_name in keep_list:
                print(f'Subset selecting, keeping: {sample_name}')
                keep_count += 1
            elif time_test:
                print(f'Subset selecting, keeping: {sample_name}')
                keep_count += 1
            else:
                shutil.move(sample, not_in_subset)
        self.zipit(not_in_subset, not_in_subset)
        self.keep_count = keep_count
        print("")

    def change_names(self):
        # self.backup_vcfs()
        change_name_file = self.change_name_file
        get_genotype_code = self.get_genotype_code
        extension = self.extension
        working_directory = self.working_directory
        starting_files = self.starting_files
        num_named_changed = 0
        num_not_changed = 0
        input_list = glob.glob(f'{working_directory}/*{extension}')
        with open(change_name_file, 'w') as name_log:
            for sample in input_list:
                sample_base = os.path.basename(sample)
                sample_name = re.sub("[._].*", "", sample_base)
                genotype_name, found = get_genotype_code(sample_name)
                if found:
                    print(f'Renamed: {sample_base}, --> {genotype_name}.{extension}', file=name_log)
                    num_named_changed += 1
                    os.rename(sample, f'{working_directory}/{genotype_name}.{extension}')
                else:
                    num_not_changed += 1
                    print(f'Name Not Changed: {sample_base}: {genotype_name}', file=name_log)
        print('Name changes:')
        print(f'\t {num_named_changed:,} file names updated')
        print(f'\t {num_not_changed:,} not updated\n')

    def backup_vcfs(self):
        starting_files = self.starting_files
        os.makedirs(starting_files)
        all_starting_files = glob.glob(f'{self.working_directory}/*{self.extension}')
        for i in all_starting_files:
            shutil.copy(i, starting_files)
        if os.path.isfile(self.change_name_file):
            shutil.move(self.change_name_file, starting_files)
        self.zipit(starting_files, starting_files) # zip starting files directory

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

    def fix_each_vcf(self, each_vcf):
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
            except IndexError:
                mal.append("##### IndexError: Deleting corrupt VCF file: " + each_vcf)
                os.remove(temp_file)
                os.remove(each_vcf)
                mal.append(os.path.basename(each_vcf))
            except UnicodeDecodeError:
                mal.append("##### UnicodeDecodeError: Deleting corrupt VCF file: " + each_vcf)
                os.remove(temp_file)
                os.remove(each_vcf)
                mal.append(os.path.basename(each_vcf))
        if os.stat(each_vcf).st_size == 0:
            mal.append("##### Empty file: Deleting corrupt VCF file: " + each_vcf)
            os.remove(temp_file)
            os.remove(each_vcf)
            mal.append(os.path.basename(each_vcf))

        write_out.close()
        os.rename(temp_file, each_vcf)
        # revert timestamp to original allows elites to properly sort on file modification time
        os.utime(each_vcf, times=(initial_file_time_stats.st_mtime, initial_file_time_stats.st_mtime))
        return mal

    def fix_vcfs(self, vcf_list=[], debug=False):
        all_mal = []
        fix_each_vcf = self.fix_each_vcf
        if debug:
            for filename in vcf_list:
                mal = fix_each_vcf(filename)
                if mal:
                    all_mal.append(mal)
        else:
            with futures.ProcessPoolExecutor(max_workers=self.cpu_count_half) as pool: #ProcessPoolExecutor ThreadPoolExecutor
                for mal in pool.map(fix_each_vcf, vcf_list):
                    if mal:
                        all_mal.append(mal)
        return all_mal

if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    
   
    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-g', '--codes', action='store', dest='excel_genotype_codes', required=True, help='Excel file containing genotype code\nColumn 1: to match sample name minus everything left of first occurring _ or .\nColumn 2: New name')
    parser.add_argument('-e', '--extension', action='store', dest='extension', required=False, default="vcf", help='File extension type to be renamed')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()

    manage = File_Management(excel_genotype_codes=args.excel_genotype_codes, extension=args.extension)
    manage.change_names()
    # update_names.backup_vcfs()
    # update_names.backup_and_change_names()