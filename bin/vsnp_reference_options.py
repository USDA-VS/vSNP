#!/usr/bin/env python

__version__ = "0.2.01"

import os
import sys
import glob
import re
import argparse
import textwrap

class Ref_Options:
    '''
    The default location for dependencies is in the "dependencies" of installed directory.  Additional directories can be added.  Paths of parent directory of dependency types are to be used.  Therefore when the following dependency types "species1" and "species2" are added the full working directory up to parent is provided as path.  Paths can be added manually or with path_adder.
    Example:
        /home/user/parent/species1
        /home/user/parent/species2
        Add: /home/user/parent to include species1 and species2 as reference options
    '''

    def __init__(self, select_ref=None):
        all_ref_options = []
        script_path = os.path.dirname(os.path.realpath(__file__))
        self.script_path = script_path
        #don't use just the script path dependencies, but gather external dependency paths too
        ref_options_file = os.path.abspath(f'{script_path}/../dependencies/reference_options_paths.txt')
        self.ref_options_file = ref_options_file
        with open(f'{ref_options_file}', 'r') as dep_paths:
            dependency_paths = [line.strip() for line in dep_paths]
        #the additional dependency paths point to more reference options
        for path in dependency_paths:
            ref_options = glob.glob(f'{path}/*')
            all_ref_options = all_ref_options + ref_options
        all_ref_options = [x for x in all_ref_options if os.path.isdir(x)] #only capture directories
        self.all_ref_options = all_ref_options
        for option in all_ref_options:
            if select_ref == option.split('/')[-1]:
                #if the asked for reference is a match grab files
                self.path = option
                excel = glob.glob(f'{option}/*xlsx')
                excel = [efile for efile in excel if not re.search('~\$*', efile)] #ignore opened files
                #there are 3 excel files.  only 1 excel file for variable "excel".  it must be the non-*meta* and non-*remove* file
                excel = [efile for efile in excel if not re.search('.*remove.*', efile)]
                remove = glob.glob(f'{option}/*remove*xlsx')
                remove = [efile for efile in remove if not re.search('~\$*', efile)] #ignore opened files
                excel = [efile for efile in excel if not re.search('.*meta.*', efile)]
                metadata = glob.glob(f'{option}/*meta*xlsx')
                metadata = [efile for efile in metadata if not re.search('~\$*', efile)] #ignore opened files
                fasta = glob.glob(f'{option}/*fasta')
                gbk = glob.glob(f'{option}/*gbk')
                #check that multiple files are not found for a single variable.  Each variable must point to just one file.
                if len(excel) > 1:
                    print(f'\n\n##### Exiting script {select_ref} contains more than one an Excel at {option}\n')
                    sys.exit(0)
                if len(metadata) > 1:
                    print(f'\n\n##### Exiting script {select_ref} contains more than one metadata file at {option}\n')
                    sys.exit(0)
                if len(remove) > 1:
                    print(f'\n\n##### Exiting script {select_ref} contains more than one remove file at {option}\n')
                    sys.exit(0)
                if len(fasta) > 1:
                    print(f'\n\n##### Exiting script {select_ref} contains more than one an FASTA at {option}\n')
                    sys.exit(0)
                if len(gbk) > 1:
                    print(f'\n\n##### Exiting script {select_ref} contains more than one an GBK at {option}n')
                    print(f'Concatenate multiple gbks to a single file with "\\\\" separator between appends')
                    sys.exit(0)

                #if not a single index 0 in list then mark as False, not having file
                try:
                    self.excel = excel[0]
                except IndexError:
                    self.excel = False
                try:
                    self.metadata = metadata[0]
                except IndexError:
                    self.metadata = False
                try:
                    self.remove = remove[0]
                except IndexError:
                    self.remove = False
                try:
                    self.fasta = fasta[0]
                except IndexError:
                    self.fasta = False
                try:
                    self.gbk = gbk[0]
                except IndexError:
                    self.gbk = False
                return
            
    def files_in_directory(self):
        all_files = glob.glob(f'{self.path}/*')
        for each_file in all_files:
            print(f'{each_file}')
        print("")

    def print_options(self):
        each_reference_option=[]
        print("\nReference option files available:")
        print(f'Path are listed here: {self.ref_options_file}')
        for option in self.all_ref_options:
            each_reference_option.append((f'\t{os.path.split(option)[-1]}'))
        for each in sorted(each_reference_option):
            print(each)
        print("\nSee vsnp_path_adder.py -h for more information\n")
        
        
if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-s', '--select_ref', action='store', dest='select_ref', required=True, help='Required: single read, R1 when Illumina read')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    select_ref = args.select_ref
    ro = Ref_Options(select_ref)
    