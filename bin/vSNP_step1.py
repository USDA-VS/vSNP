#!/usr/bin/env python

__version__ = "2.0.0"

import os
import sys
import shutil
import glob
import pandas as pd
import argparse
import textwrap

from vsnp_fastq_quality import FASTQ_Quality
from vsnp_alignment_vcf import Align_Reads
from vsnp_reference_options import Ref_Options
from vsnp_best_reference import Best_Reference


parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

---------------------------------------------------------
Three Senario Options:
    Senario 1: Provide a FASTA
    vSNP_step1.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz -r *fasta
    vSNP_step1.py -r1 *fastq.gz -r *fasta

    Senario 2: Provide a Reference Option
    Run -t option to see table of reference options: vSNP_step1.py -t
    vSNP_step1.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz -r Mycobacterium_AF2122

    Senario 3: Find Best Reference.  Only for TB complex, paraTB, and Brucella.
    Run without -r option
    vSNP_step1.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz

FASTA, gbk, and gff files for multi-chromosome genomes must be concatenated to single file

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
    - vSNP comes installed with Mycobacterium_AF2122.  For additional reference options see: https://github.com/USDA-VS/vSNP_dependencies.git


'''), epilog='''---------------------------------------------------------''')

    
parser.add_argument('-r1', '--read1', action='store', dest='read1', required=False, help='Required: single read, R1 when Illumina read')
parser.add_argument('-r2', '--read2', action='store', dest='read2', required=False, default=None, help='Optional: R2 Illumina read')
parser.add_argument('-r', '--reference', action='store', dest='reference', required=False, default=None, help="Optional: Provide reference option or FASTA file.  If neither are given, no -r option, then a TB/Brucella/paraTB best reference are searched")
parser.add_argument('-g', '--gbk', action='store', dest='gbk', required=False, default=None, help='Optional: gbk to annotate VCF file')
parser.add_argument('-t', '--table', action='store_true', dest='table', help='See reference options')
parser.add_argument('-group', '--group', action='store', dest='group', required=False, default=None, help="Optional: group output via best_reference.py, ie specify TB or Bruc which initials spoligo or MLST generation")
parser.add_argument('-skip_assembly', '--skip_assembly', action='store_true', dest='skip_assembly', help='Optional: skip assembly of unmapped reads')
parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

args = parser.parse_args()
read1 = args.read1
read2 = args.read2
gbk = args.gbk
species = None #species output via best_reference.py, ie specify directory name with dependents to be reported in Excel stats file.  When no species from best_reference.py FASTA name is used.
group = args.group
skip_assembly = args.skip_assembly

if args.table:
    reference_options = Ref_Options(None)
    reference_options.print_options()
    sys.exit(0)

if not read1:
    print(f'Read1 required')
    sys.exit(1)

#IF -r REFERENCE FASTA PROVIDED USE IT, ie .fasta
cwd = os.getcwd()
try:
    if args.reference.lower().endswith(('.fasta', '.fa', '.fas')):
        try:
            shutil.copy2(args.reference, cwd)
            reference = f'{cwd}/{os.path.basename(args.reference)}'
        except shutil.SameFileError:
            reference = f'{cwd}/{os.path.basename(args.reference)}'
    elif args.reference:
        #IF -r SPECIFIES A REFERENCE OPTION USE IT, ie directory name
        reference_options = Ref_Options(args.reference)
        try:
            shutil.copy2(reference_options.fasta, cwd)
            reference = f'{cwd}/{os.path.basename(reference_options.fasta)}'
        except AttributeError:
            print(f'{args.reference} reference type not available, see available options:')
    else:
        print(f'Bad argument, -r did not trigger an AttributeError')
#IF NO REFERENCE PROVIDED SEEK A "BEST REFERNCE", for tb, brucella and para
except AttributeError:
    best_ref = Best_Reference(read1, read2)
    print(f'{read1} Finding Best Reference')
    best_ref.best_reference()
    print(f'### Best Reference: {best_ref.group}, {best_ref.species}')
    reference_options = Ref_Options(best_ref.species)
    species = best_ref.species
    group = best_ref.group
    try:
        shutil.copy2(reference_options.fasta, cwd)
        reference = f'{cwd}/{os.path.basename(reference_options.fasta)}'
    except AttributeError:
        print(f'Best Reference not found: {best_ref.group}, {best_ref.species}')
        fq = FASTQ_Quality(read1, read2)
        print(f'\nGetting FASTQ quality for stats file...')
        fq.get_quality()
        print(f'Printing FASTQ stats to excel')
        fq.excel_fastq_only_stats_out(f'{fq.sample_name}_fastq_stats.xlsx')
        print(f'{args.reference} reference type not available, see options:')
        reference_options.print_options()
        sys.exit(0)  

align_reads = Align_Reads(read1, read2, reference, gbk, species, group, skip_assembly)
align_reads.align()
print(f'{os.path.basename(read1)} done')