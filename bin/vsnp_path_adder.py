#!/usr/bin/env python

__version__ = "2.0.0"

import os
import sys
import glob
import re
import argparse
import textwrap

class Add_Path:
    '''
    Dependencies:
        Defining SNP/Filter file, xlsx - REQUIRED
        FASTA - REQUIRED
        GBK - OPTIONAL
        GFF - OPTIONAL
        metadata.xlsx (2 column: filename and updated name) - OPTIONAL
    '''

    def __init__(self, added_directory):
        if added_directory == '.':
            added_directory = os.getcwd()
        else:
            added_directory = added_directory
        try: 
            if list(added_directory)[0] == '/':
                print(f'\nDirectory Given: {added_directory}\n')
                self.valid_path = True
            else:
                self.valid_path = False
        except TypeError:
            # If None is passed error because not iterable
            pass
        script_path = os.path.dirname(os.path.realpath(__file__))
        self.script_path = script_path
        abspath = os.path.abspath(f'{script_path}/../dependencies/reference_options_paths.txt')
        self.abspath = abspath
        self.added_directory = added_directory

    #insure no blank lines are in file
    def remove_blank_lines(self):
        with open(f'{self.abspath}',"r") as open_file:
            all_lines=open_file.readlines()
        with open(f'{self.abspath}',"w") as write_back:  
            [write_back.write(line) for line in all_lines if line.strip() ] 

    def add_to_path(self,):
        #append new path to file
        with open(f'{self.abspath}', 'a') as dep_paths:
            print(f'{self.added_directory}', file=dep_paths)

        # Remove duplicate lines from file
        with open(self.abspath, "r") as infile:
            all_lines = [line.rstrip() for line in infile]
        unique_lines = set(all_lines)
        with open(self.abspath, "w") as outfile:
            for each_line in unique_lines:
                print(f'{each_line}', file=outfile)

    def show_options(self,):
        # Show reference options
        with open(f'{self.abspath}', 'r') as dep_paths:
            reference_options_paths = [line.rstrip() for line in dep_paths]
        for path in reference_options_paths:
            print(f'Path: {path}')
            ref_options = glob.glob(f'{path}/*')
            ref_options = [x for x in ref_options if os.path.isdir(x)]
            if ref_options:
                for option in sorted(ref_options):
                    print(f'\t{os.path.basename(option)}')
                print(f'\n')
            elif not os.path.exists(path):
                print(f'\tPath does not exist\n')
            else:
                print(f'\tPath is empty\n')
             

if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------
    Using no arguments or -s option show the same output.
    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-d', '--cwd', action='store', dest='directory', required=False, help='Absolute directory path to be added to find reference option files.')
    parser.add_argument('-s', '--show', action='store_true', dest='show', required=False, help='Show available directories.')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')
    args = parser.parse_args()

    if args.directory:
        directory = args.directory
        add_path = Add_Path(directory)
        add_path.remove_blank_lines()
        if add_path.valid_path:
            add_path.add_to_path()
            add_path.show_options()
            print(f'\nPaths files can be manually edited here: {add_path.abspath}\n')
        else:
            print(f'##### PROVIDE AN ABSOLUTE PATH')
            print(f'Directory Given: "{add_path.added_directory}"')
            print(f'\nPaths files can be manually edited here: {add_path.abspath}\n')
    else:
        print(f'\nNO NEW DIRECTORY ADDED\nONLY SHOWING WHAT IS CURRENTLY AVAILABLE')
        add_path = Add_Path(None)
        add_path.remove_blank_lines()
        if os.stat(f'{add_path.abspath}').st_size == 0:
            print(f'\nPaths have not been added.  File is empty. \nAdd path using -d option.\nSee -h for more options\n\t{add_path.abspath}\n')
        else:
            add_path.show_options()
            print(f'\nPaths files can be manually edited here: {add_path.abspath}\n')