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
    The default location for dependencies is in the "dependencies" of installed directory.  Additional directories can be added by adding the full directory path to dependency_paths.txt.  Paths of parent directory of dependency types are to be used.  Therefore when the following dependency types "species1" and "species2" are added the full working directory up to parent is provided as path.  Paths can be added manually or with path_adder.
    Example:
        /home/user/parent/species1
        /home/user/parent/species2
        Add: /home/user/parent to include species1 and species2 as reference options

    Dependencies:
        Defining SNP/Filter file, xlsx - REQUIRED
        FASTA - REQUIRED
        GBK - OPTIONAL
        GFF - OPTIONAL
        metadata.xlsx (2 column: filename and updated name) - OPTIONAL
    '''

    def __init__(self, working_directory):

        if args.working_directory == '.':
            working_directory = os.getcwd()
        else:
            working_directory = args.working_directory
        if list(working_directory)[0] == '/':
            print(f'working directory: {working_directory}')
        else:
            print(f'##### PROVIDE A FULL PATH')
            print(f'directory given: "{working_directory}"')
            sys.exit(0)
        script_path = os.path.dirname(os.path.realpath(__file__))
        self.script_path = script_path
        with open(f'{script_path}/../dependencies/dependency_paths.txt', 'a') as dep_paths:
            print(f'{working_directory}', file=dep_paths)

        with open(f'{script_path}/../dependencies/dependency_paths.txt', 'r') as dep_paths:
            dependency_paths = [line.strip() for line in dep_paths]
        ref_options = glob.glob(f'{script_path}/dependencies/*')
        for path in dependency_paths:
            more_ref_options = glob.glob(f'{path}/*')
            ref_options = ref_options + more_ref_options
        print(f'Current ref paths:')
        for option in ref_options:
            print(f'\t{option}')        

if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-w', '--cwd', action='store', dest='working_directory', required=False, default='.', help='Working directoy to be added as path to dependency files... aka add a new reference')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    working_directory = args.working_directory
    Add_Path(working_directory)