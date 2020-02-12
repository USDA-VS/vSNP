#!/usr/bin/env python

__version__ = "2.03"

import os
import sys
import shutil
import gzip
import glob
import time
import re
import numpy as np
import pandas as pd
import zipfile
import pysam
import vcf
import humanize
import argparse
import textwrap
from numpy import mean
from datetime import datetime
from Bio import SeqIO


class FASTQ_Container:
    #Provide nested dot notation to object for each read with stats, fq.read1.fastq --> 'sample_S25_L001_R1.fastq.gz'
    def __init__(self, sample_name, fastq, file_size, total_read_count, sampling_size, length_mean, read_average, reads_gt_q30):
        self.sample_name = sample_name
        self.fastq = fastq
        self.file_size = file_size
        self.total_read_count = total_read_count
        self.sampling_size = sampling_size
        self.length_mean = length_mean
        self.read_average = read_average
        self.reads_gt_q30 = reads_gt_q30


class FASTQ_Quality:
    '''
    Files must be .gz zipped
    Paired reads must contain _R1 or _R2 in file name

    Along with the FASTQ file and sample name 5 FASTQ attributes are obtained:
        file_size --> FASTQ file size (human readable)
        total_read_count --> Total read count within each FASTQ file
        sampling_size  --> Random analyzed read count
        length_mean --> Average read length
        read_average --> Average read quality
        reads_gt_q30 --> Read counts with an average quality greater than 30

    Note: when calculating percent of reads above Q30 use reads_gt_q30/sampling_size

    After ran object will contain nested dot notation for each read, fq.read1.fastq --> 'sample_S25_L001_R1.fastq.gz'
    '''
    def __init__(self, read1, read2=None, sampling_number=10000):
        quality_key = {'!':'0', '"':'1', '#':'2', '$':'3', '%':'4', '&':'5', "'":'6', '(':'7', ')':'8', '*':'9', '+':'10', ',':'11', '-':'12', '.':'13', '/':'14', '0':'15', '1':'16', '2':'17', '3':'18', '4':'19', '5':'20', '6':'21', '7':'22', '8':'23', '9':'24', ':':'25', ';':'26', '<':'27', '=':'28', '>':'29', '?':'30', '@':'31', 'A':'32', 'B':'33', 'C':'34', 'D':'35', 'E':'36', 'F':'37', 'G':'38', 'H':'39', 'I':'40', 'J':'41', 'K':'42', 'L':'43', 'M':'44', 'N':'45', 'O':'46', 'P':'47', 'Q':'48', 'R':'49', 'S':'50', 'T':'51', 'U':'52', 'V':'53', 'W':'54', 'X':'55', 'Y':'56', 'Z':'57', '_':'1', ']':'1', '[':'1', '\\':'1', '\n':'1', '`':'1', 'a':'1', 'b':'1', 'c':'1', 'd':'1', 'e':'1', 'f':'1', 'g':'1', 'h':'1', 'i':'1', 'j':'1', 'k':'1', 'l':'1', 'm':'1', 'n':'1', 'o':'1', 'p':'1', 'q':'1', 'r':'1', 's':'1', 't':'1', 'u':'1', 'v':'1', 'w':'1', 'x':'1', 'y':'1', 'z':'1', ' ':'1'}
        self.fastq_list = []
        for read in [read1, read2]:
            if read: #append if not None
                self.fastq_list.append(read)
        if read2:
            self.paired = True
        self.sample_name = re.sub('_.*', '', os.path.basename(read1))
        self.sampling_number = sampling_number
        self.root_dir = str(os.getcwd())
        self.quality_key = quality_key

    
    def fastq_stats(self, fastq):
        #Determine read type
        read1_pattern = re.compile('.*_R1.*')
        read2_pattern = re.compile('.*_R2.*')
        if read1_pattern.match(fastq):
            read = "read1"
        elif read2_pattern.match(fastq):
            read = "read2"
        else:
            read = "read1"
        file_size = humanize.naturalsize(os.path.getsize(fastq))
        df = pd.read_csv(gzip.open(fastq, "rt"), header=None, sep='^') #basically set sep to None
        #Starting at row 3, keep every 4 row
        #Random sample specified number of rows
        total_read_count = int(len(df.index)/4)
        sampling_size = int(self.sampling_number)
        if sampling_size > total_read_count:
            sampling_size = total_read_count
        df = df.iloc[3::4].sample(sampling_size)
        dict_mean={}
        list_length=[]
        for index, row in df.iterrows():
            base_qualities=[]
            for base in list(row.array[0]):
                base_qualities.append(int(self.quality_key[base]))
            dict_mean[index] = np.mean(base_qualities)
            list_length.append(len(row.array[0]))
        length_mean = np.mean(list_length)
        df_mean = pd.DataFrame.from_dict(dict_mean, orient='index', columns=['ave'])
        read_average = df_mean['ave'].mean()
        reads_gt_q30 = len(df_mean[df_mean['ave'] >= 30])
        print(f'{fastq} --> File Size: {file_size}, Total Reads: {total_read_count:,}, Random Analyzed Reads: {sampling_size:,}, Mean Read Length: {length_mean:.1f}, Mean Read Quality: {read_average:.1f}, Reads Passing Q30: {reads_gt_q30/sampling_size:0.1%}')
        container = FASTQ_Container(self.sample_name, fastq, file_size, total_read_count, sampling_size, length_mean, read_average, reads_gt_q30)
        return read, container

    def get_quality(self):

        for myfastq in self.fastq_list:
            read, container = self.fastq_stats(myfastq)
            setattr(self, read, container)

    def excel_fastq_only_stats_out(self, out_file):
        #stats to excel
        df = pd.DataFrame(index=[self.sample_name], columns=['Reference', 'Read1 FASTQ', 'Read1 File Size', 'Read1 Total Reads', 'Read1 Mean Read Length', 'Read1 Mean Read Quality', 'Read1 Reads Passing Q30', 'Read2 FASTQ', 'Read2 File Size', 'Read2 Total Reads', 'Read2 Mean Read Length', 'Read2 Mean Read Quality', 'Read2 Reads Passing Q30', 'Total Reads', 'All Mapped Reads', 'Reference with Coverage', 'Average Depth of Coverage', 'Unmapped Reads', 'Unmapped Assembled Contigs', 'Good SNP Count'])
        df.at[self.sample_name, 'Reference'] = 'N/A'
        df.at[self.sample_name, 'Read1 FASTQ'] = f'{self.read1.fastq}'
        df.at[self.sample_name, 'Read1 File Size'] = f'{self.read1.file_size}'
        df.at[self.sample_name, 'Read1 Total Reads'] = f'{self.read1.total_read_count:,}'
        df.at[self.sample_name, 'Read1 Mean Read Length'] = f'{self.read1.length_mean:.1f}'
        df.at[self.sample_name, 'Read1 Mean Read Quality'] = f'{self.read1.read_average:.1f}'
        df.at[self.sample_name, 'Read1 Reads Passing Q30'] = f'{self.read1.reads_gt_q30/self.read1.sampling_size:0.1%}'
        total_reads = self.read1.total_read_count
        try:
            df.at[self.sample_name, 'Read2 FASTQ'] = f'{self.read2.fastq}'
            df.at[self.sample_name, 'Read2 File Size'] = f'{self.read2.file_size}'
            df.at[self.sample_name, 'Read2 Total Reads'] = f'{self.read2.total_read_count:,}'
            df.at[self.sample_name, 'Read2 Mean Read Length'] = f'{self.read2.length_mean:.1f}'
            df.at[self.sample_name, 'Read2 Mean Read Quality'] = f'{self.read2.read_average:.1f}'
            df.at[self.sample_name, 'Read2 Reads Passing Q30'] = f'{self.read2.reads_gt_q30/self.read2.sampling_size:0.1%}'
            total_reads = self.read1.total_read_count + self.read2.total_read_count
        except AttributeError:
            #when no read2
            df.at[self.sample_name, 'Read2 FASTQ'] = 'NA'
            df.at[self.sample_name, 'Read2 File Size'] = 'NA'
            df.at[self.sample_name, 'Read2 Total Reads'] = 'NA'
            df.at[self.sample_name, 'Read2 Mean Read Length'] = 'NA'
            df.at[self.sample_name, 'Read2 Mean Read Quality'] = 'NA'
            df.at[self.sample_name, 'Read2 Reads Passing Q30'] = 'NA'
        df.at[self.sample_name, 'Total Reads'] = f'{total_reads:,}'
        # df.at[fq.sample_name, 'All Mapped Reads'] = f'{self.allbam_mapped_reads:,}'
        # df.at[fq.sample_name, 'Reference with Coverage'] = f'{self.genome_coverage}'
        # df.at[fq.sample_name, 'Average Depth of Coverage'] = f'{self.ave_coverage:.1f}X'
        # df.at[fq.sample_name, 'Unmapped Reads'] = f'{self.unmapped_reads:,}'
        # df.at[fq.sample_name, 'Unmapped Assembled Contigs'] = f'{self.abyss_contig_count:,}'
        # df.at[fq.sample_name, 'Good SNP Count'] = f'{self.good_snp_count:,}'
        df.index.name = 'sample'
        df.to_excel(out_file)


if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
   
    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-r1', '--read1', action='store', dest='read1', required=True, help='Required: single read, R1 if Illumina read')
    parser.add_argument('-r2', '--read2', action='store', dest='read2', required=False, default=None, help='Required: R2 Illumina read')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    read1 = args.read1
    read2 = args.read2

    fq = FASTQ_Quality(read1, read2)
    fq.get_quality()
    fq.excel_fastq_only_stats_out(f'{fq.sample_name}_fastq_stats.xlsx')