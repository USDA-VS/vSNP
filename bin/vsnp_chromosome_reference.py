#!/usr/bin/env python

__version__ = "2.0.0"

import os
import sys
import glob
import re
import argparse
import textwrap
import vcf

class Reference_Chromosome:
    '''
    Single VCF file will be chosen, and chromosome read (first column).  Chromosome will be looked up in dictionary and return the reference type.

    '''

    def __init__(self, working_directory):

        species_cross_reference={}
        species_cross_reference[""] = ""
        species_cross_reference["NC_006932.1"] = "Brucella_abortus1"
        species_cross_reference["NC_006933.1"] = "Brucella_abortus1"
        species_cross_reference["NZ_CP007682.1"] = "Brucella_abortus3"
        species_cross_reference["NZ_CP007683.1"] = "Brucella_abortus3"
        species_cross_reference["NC_010103.1"] = "Brucella_canis"
        species_cross_reference["NC_010104.1"] = "Brucella_canis"
        species_cross_reference["Bceti1Cudo"] = "Brucella_ceti1"
        species_cross_reference["NC_022905.1"] = "Brucella_ceti2"
        species_cross_reference["NC_022906.1"] = "Brucella_ceti2"
        species_cross_reference["NC_003317.1"] = "Brucella_melitensis-bv1"
        species_cross_reference["NC_003318.1"] = "Brucella_melitensis-bv1"
        species_cross_reference["NZ_CP018508.1"] = "Brucella_melitensis-bv1b"
        species_cross_reference["NZ_CP018509.1"] = "Brucella_melitensis-bv1b"
        species_cross_reference["NC_012441.1"] = "Brucella_melitensis-bv2"
        species_cross_reference["NC_012442.1"] = "Brucella_melitensis-bv2"
        species_cross_reference["NZ_CP007760.1"] = "Brucella_melitensis-bv3"
        species_cross_reference["NZ_CP007761.1"] = "Brucella_melitensis-bv3"
        species_cross_reference["KN046827.1"] = "Brucella_neotomae"
        species_cross_reference["NC_009504.1"] = "Brucella_ovis"
        species_cross_reference["NC_009505.1"] = "Brucella_ovis"
        species_cross_reference["NC_017251.1"] = "Brucella_suis1"
        species_cross_reference["NC_017250.1"] = "Brucella_suis1"
        species_cross_reference["NC_010169.1"] = "Brucella_suis2"
        species_cross_reference["NC_010167.1"] = "Brucella_suis2"
        species_cross_reference["NZ_CP007719.1"] = "Brucella_suis3"
        species_cross_reference["NZ_CP007718.1"] = "Brucella_suis3"
        species_cross_reference["B-REF-BS4-40"] = "Brucella_suis4"
        species_cross_reference["NC_002945.4"] = "Mycobacterium_AF2122"
        species_cross_reference["NC_000962.3"] = "Mycobacterium_H37"
        species_cross_reference["CP033688.1"] = "para-CP033688"
        species_cross_reference["NC_002944.2"] = "para-NC002944"
        species_cross_reference["Gua1_1407_2018"] = "guatemala"
        self.species_cross_reference = species_cross_reference

        if working_directory == '.':
            working_directory = os.getcwd()

        if list(working_directory)[0] == '/':
            print(f'working directory: {working_directory}')
        else:
            print(f'##### PROVIDE A FULL PATH')
            print(f'directory given: "{working_directory}"')
            sys.exit(0)
        self.working_directory = working_directory

    def get_reference(self):
        vcf_file = glob.glob(f'{self.working_directory}/*vcf')[0]
        print(f'Looking for chrom in: {vcf_file}')
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))
        try:
            for record in vcf_reader:
                for key, value in self.species_cross_reference.items():
                    if record.CHROM == key:
                        return (value)
        except ValueError:
            #if corrupt VCF go to the next file
            pass


if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-w', '--cwd', action='store', dest='working_directory', required=False, default='.', help='Optional: path to VCF files')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    working_directory = args.working_directory
    reference_type = Reference_Chromosome(working_directory)
    reference = reference_type.get_reference()
    print(f'Reference type found: {reference}')