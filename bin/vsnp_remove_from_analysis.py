#!/usr/bin/env python

__version__ = "0.2.01"

import os
import sys
import pandas as pd
import glob
import argparse
import textwrap

class Remove_From_Analysis:
    '''
    ''' 
    def __init__(self, working_directory='.', excel_remove=None, extension="vcf"):
        if working_directory == '.':
            working_directory = os.getcwd()
        if list(working_directory)[0] == '/':
            print(f'working directory: {working_directory}')
        else:
            print(f'##### PROVIDE A FULL PATH')
            print(f'directory given: "{working_directory}"')
            sys.exit(0)

        df = pd.read_excel(excel_remove, index_col=0, usecols=[0], header=None)
        remove_list=[]
        for each_sample in df.index:
            remove_list.append(f'{working_directory}/{each_sample}*.{extension}')
        self.excel_remove = excel_remove
        self.remove_list = remove_list

    def remove_files(self):
        num_files_removed = 0
        print(f'Removing samples listed in {self.excel_remove}')
        for each_sample in self.remove_list:
            glob_list = glob.glob(each_sample)
            for item in glob_list:
                num_files_removed += 1
                # print(f'\tRemoving: {item}')
                os.remove(item)
        print(f'{num_files_removed:,} files removed from the analysis\n')
        self.removed_file_count = num_files_removed


if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-r', '--codes', action='store', dest='excel_remove', required=True, help='Excel file containing samples to remove from analysis\nColumn 1: to match sample name minus extension.\nNo header allowed')
    parser.add_argument('-w', '--cwd', action='store', dest='working_directory', required=False, default='.', help='Optional: path to VCF files')
    parser.add_argument('-e', '--extension', action='store', dest='extension', required=False, default="vcf", help='File extension type to be renamed')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()

    remove_from_analysis = Remove_From_Analysis(working_directory=args.working_directory, excel_remove=args.excel_remove, extension=args.extension)
    remove_from_analysis.remove_files()