#!/usr/bin/env python

import os
import sys
import shutil
import multiprocessing
import argparse
import textwrap
import subprocess
import gzip
import glob
import csv
import json
import time
import regex
import re
import numpy as np
import pandas as pd
import zipfile
import xlsxwriter
import xlrd
import pysam
import vcf
import smtplib
from multiprocessing import Pool
from dask import delayed
from itertools import repeat as itertools_repeat
from collections import Iterable
from numpy import mean
from email.utils import formatdate
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders
from distutils.dir_util import copy_tree
from datetime import datetime
from concurrent import futures
from collections import OrderedDict
from collections import Counter
from collections import defaultdict
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO

from parameters import Get_Specie_Parameters
import functions

root_dir = str(os.getcwd())

parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

---------------------------------------------------------
vSNP --> get SNPs, group SNPs, verify SNPs

vSNP is called on a working directory containing FASTQ or VCF files.

See documentation at: https://usda-vs.github.io/snp_analysis/

        Step 1: FASTQs --> VCF

        Step 2: VCFs --> Tables & Trees

'''), epilog='''---------------------------------------------------------''')

#universal
parser.add_argument('-s', '--species', action='store', dest='species', help='OPTIONAL: Used to FORCE species type <see options above>')
parser.add_argument('-d', '--debug', action='store_true', dest='debug_call', help='debug, run without pool.map')
parser.add_argument('-g', '--get', action='store_true', dest='get', help='get, get to the core functions for debugging')
parser.add_argument('-n', '--no_annotation', action='store_true', dest='no_annotation', help='no_annotation, run without annotation')
parser.add_argument('-a', '--all_vcf', action='store_true', dest='all_vcf', help='make tree using all VCFs')
parser.add_argument('-o', '--only_all_vcf', action='store_true', dest='only_all_vcf', help='make tree using all VCFs')
parser.add_argument('-e', '--elite', action='store_true', dest='elite', help='create a tree with on elite sample representation')
parser.add_argument('-f', '--filter', action='store_true', dest='filter_finder', help='Find possible positions to filter')
parser.add_argument('-p', '--processor', action='store', dest='processor', help='max processor usage')
parser.add_argument('-q', '--quiet', action='store_true', dest='quiet', help='[**APHIS only**] prevent stats going to cumlative collection')
parser.add_argument('-m', '--email', action='store', dest='email', help='[**APHIS only**, specify own SMTP address for functionality] email options: all, s, tod, jess, suelee, chris, email_address')
parser.add_argument('-u', '--upload', action='store_true', dest='upload', help='[**APHIS only**, specify own storage for functionality] upload files to the bioinfo drive')
parser.add_argument('-t', '--table', action='store_true', dest='table', help='print reference/species table')
parser.add_argument('-i', '--ignore_filters', action='store_true', dest='ignore_filters', help='print reference/species table')
parser.add_argument('-l', '--label', action='store', dest='label', help='[Step 2 table filter.  Provide annotated label. Ex: "ppe|pgrs|repeat".  greedy/case-insensitive')
args = parser.parse_args()

if args.table:
    pretty_table = functions.reference_table()
    print(pretty_table)
    sys.exit()

if args.only_all_vcf:
    args.all_vcf = True
print("\nSET ARGUMENTS: ")
print(args)
arg_options = {
    "species": args.species,
    "debug_call": args.debug_call,
    "get": args.get,
    "no_annotation": args.no_annotation,
    "all_vcf": args.all_vcf,
    "only_all_vcf": args.only_all_vcf,
    "elite": args.elite,
    "filter_finder": args.filter_finder,
    "processor": args.processor,
    "quiet": args.quiet,
    "upload": args.upload,
    "ignore_filters": args.ignore_filters,
    "label": args.label,
}
print("")

email_dict = {}
email_dict["all"] = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov, patrick.m.camp@aphis.usda.gov, David.T.Farrell@aphis.usda.gov, Robin.L.Swanson@aphis.usda.gov, Doris.M.Bravo@aphis.usda.gov, eto3@cdc.gov, kristina.lantz@aphis.usda.gov, Tyler.Thacker@aphis.usda.gov"
email_dict["tod"] =  "tod.p.stuber@aphis.usda.gov"
email_dict["jess"] =  "Jessica.A.Hicks@aphis.usda.gov"
email_dict["suelee"] =  "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Tyler.Thacker@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov, Doris.M.Bravo@aphis.usda.gov, kristina.lantz@aphis.usda.gov, patrick.m.camp@aphis.usda.gov"
email_dict["suelee-"] =  "tod.p.stuber@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov"
email_dict["tyler-"] =  "tod.p.stuber@aphis.usda.gov, Tyler.Thacker@aphis.usda.gov"
email_dict["chris"] =  "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov, eto3@cdc.gov, kristina.lantz@aphis.usda.gov, Tyler.Thacker@aphis.usda.gov, patrick.m.camp@aphis.usda.gov"
email_dict["chris-"] =  "tod.p.stuber@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov"
email_dict["kris"] =  "kristina.lantz@aphis.usda.gov, tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov, mary.k.smith@aphis.usda.gov, patrick.m.camp@aphis.usda.gov"
email_dict["doris"] =  "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, doris.m.bravo@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov, kristina.lantz@aphis.usda.gov, patrick.m.camp@aphis.usda.gov"

arg_options['email_list'] = email_dict.get(args.email, None)

################################################################################################################################################

all_file_types = glob.glob('*.*')
all_file_types_count = len([x for x in all_file_types if not re.match(r'.*log', x)]) #don't include immediately made .log files
fastq_check = len(glob.glob('*fastq.gz'))
vcf_check = len(glob.glob('*vcf'))

# Check that there are either FASTQs or VCFs, not both
if fastq_check > 0:
    fastq_check = True
if vcf_check > 0:
    vcf_check = True
if fastq_check and vcf_check:
    print("\n#####You have a mix of FASTQ and VCF files.  This is not allowed\n\n")
    sys.exit(0)

arg_options['root_dir'] = root_dir

if arg_options['processor']:
    cpu_count = int(args.processor)
else:
    cpu_count = int(multiprocessing.cpu_count())
limited_cpu_count = int(cpu_count / 4)
if limited_cpu_count == 0:
    limited_cpu_count = 1

arg_options['cpu_count'] = cpu_count
arg_options['limited_cpu_count'] = limited_cpu_count

# Check that there an equal number of both R1 and R2 reads
if fastq_check:
    # Pair check
    pair_check = len(glob.glob('*_R2*fastq.gz'))
    if pair_check > 0:
        R1 = glob.glob('*_R1*fastq.gz')
        R2 = glob.glob('*_R2*fastq.gz')
    else:
        R1 = glob.glob('*fastq.gz')
        R2 = None

    if arg_options['all_vcf'] or arg_options['elite'] or arg_options['upload'] or arg_options['filter_finder']:
        print("#####Incorrect use of options when running loop/script 1")
        sys.exit(0)

    print("\n--> RUNNING LOOP/SCRIPT 1\n")
    #Enter script 1 -->
    functions.run_loop(arg_options)
    print("See files, vSNP has finished alignments")
elif vcf_check:
    #fix files
    malformed = []
    vcf_list = glob.glob('*vcf')
    print("Fixing files...\n")
    if arg_options['debug_call'] and not arg_options['get']:
        for each_vcf in vcf_list:
            print(each_vcf)
            mal = functions.fix_vcf(each_vcf, arg_options)
            malformed = list(mal)
    else:
        with futures.ProcessPoolExecutor() as pool:
            mal = pool.map(functions.fix_vcf, vcf_list, itertools_repeat(arg_options))
            malformed = malformed + list(mal)
    malformed = [x for x in malformed if x] # remove blanks
    print("done fixing")
    arg_options['malformed'] = malformed

    if not arg_options['species']:
        species = functions.get_species(arg_options)
        if species is None:
            print("\nEXITED\n##### Unable to find a species corresponding to CHROM found in VCF files")
            sys.exit(0)
        arg_options['species'] = species
        print("species %s" % species)
    vcfs_count = len(glob.glob('*vcf'))
    if (all_file_types_count != vcfs_count):
        print("\n##### You have more than just VCF files in your directory.  Only VCF files are allowed if running script 2\n\n")
        sys.exit(0)
    else:
        if arg_options['quiet']:
            print("#####Incorrect use of options when running script 2, when running step 2 -q cannot be used")
            sys.exit(0)
        else:
            if arg_options['species']:
                print("\n--> RUNNING SCRIPT 2\n")
                #Enter script 2 -->
                functions.run_script2(arg_options)
            else:
                print("#####Based on VCF CHROM id (reference used to build VCF) a matching species cannot be found neither was there a -s option given")
                sys.exit(0)
else:
    print("#####Error determining file type.")
    sys.exit(0)
