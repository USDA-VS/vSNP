#!/usr/bin/env python

import os
import sys
import re
import allel
import argparse
import textwrap
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
from Bio import SeqIO


class FileName:
    '''
    return the file and name as object.  name is string before first "."
    usage: FileName("myfile.txt")
    '''
    def __init__(self, file):
        self.file = file
        self.name = re.sub('[.].*', '', file)


class FastaLineBreaks:
    '''
    return FASTA file with sequence with line breaks every 60 characters
    '''
    def __init__(self, nobreak_file):
        records = SeqIO.parse(nobreak_file, "fasta")
        SeqIO.write(records, nobreak_file + "-temp", "fasta")
        os.replace(nobreak_file + "-temp", nobreak_file)

class Merge_VCF_to_FASTA:
    '''
    Setup optional parameters
    '''
    def __init__(self, qual_threshold, map_threshold, ambiguity_NOT):
        self.qual_threshold = qual_threshold
        self.map_threshold = map_threshold
        self.ambiguity_NOT = ambiguity_NOT

    def merge_vcf_into_fasta(self, fasta_file, vcf_file):
        '''

        '''
        fasta = FileName(fasta_file)
        vcf = FileName(vcf_file)

        output_file_name = vcf.name + "_merged_into_" + fasta.name + ".fasta"
        fasta_out = open(output_file_name, 'w')
        print(f'Threshold cutoff, not applying calls below {self.qual_threshold} QUAL')
        print(f'Threshold cutoff, not applying calls below {self.map_threshold} MQ')
        if self.ambiguity_NOT:
            print(f'AC=1 calls will NOT be reported with IUPAC ambiguity nucleotide codes.')
        else:
            print(f'AC=1 calls will be reported with IUPAC ambiguity nucleotide codes')
        

        for seq_record in SeqIO.parse(fasta.file, "fasta"):
            qual_threshold = int(self.qual_threshold)
            map_threshold = int(self.map_threshold)
            ambiguity_NOT = self.ambiguity_NOT

            chrom_seq_dict = {}
            record = 0
            chrom = seq_record.id
            sequence = seq_record.seq
            for base in sequence:
                record += 1 
                chrom_seq_dict[chrom + "-" + str(record)] = base
            fasta_file_df = pd.DataFrame.from_dict(chrom_seq_dict, orient='index')
            fasta_file_df.columns = ["REF"]
            fasta_file_df = fasta_file_df.rename(columns={"REF": "ALT"})
            fasta_file_df.index.name = "Absolute_Pos"

            vcf_df = allel.vcf_to_dataframe(vcf.file, fields=['variants/CHROM', 'variants/POS', 'variants/QUAL', 'variants/REF', 'variants/ALT', 'variants/AC', 'variants/DP', 'variants/MQ'], alt_number=1)
            df1 = vcf_df[((vcf_df.QUAL >= qual_threshold) & (vcf_df.ALT.str.len() == 1) & (vcf_df.MQ >= map_threshold) | (vcf_df.REF == "N"))]
            if not ambiguity_NOT:  #aka applying ambiguity
                df1.ALT = np.where(df1.AC.eq(1), df1.REF + df1.ALT, df1.ALT)  #cat to ALT column when AC=1
                df1.ALT = df1.ALT.replace({"AG": "R", "CT": "Y", "GC": "S", "AT": "W", "GT": "K", "AC": "M", "GA": "R", "TC": "Y", "CG": "S", "TA": "W", "TG": "K", "CA": "M"})
            df1.ALT = np.where(df1.REF.eq("N"), df1.REF, df1.ALT) #move the Ns to the ALT column
            df1["Absolute_Pos"] = df1["CHROM"].map(str) + "-" + df1["POS"].map(str)
            vcf_merge_read = df1[["Absolute_Pos", "ALT"]]
            vcf_merge_read = vcf_merge_read.set_index('Absolute_Pos')

            fasta_file_df.update(vcf_merge_read)  #merge vcf changes to fasta dataframe
            print(">{}".format(vcf.name + "_merged_into_" + fasta.name), file=fasta_out)
            print("{}".format(''.join(list(fasta_file_df.to_dict()['ALT'].values()))), file=fasta_out)

        print(f'\nFile written to: {output_file_name}\n')
        fasta_out.close()
        FastaLineBreaks(output_file_name)


if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Merge VCF changes into FASTA sequence.
    Usage: merge_vcf_into_fasta.py -r <fasta_file> -v <vcf_file>

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-f', '--fasta_file', action='store', dest='fasta_file', required=True, help='REQUIRED: In file to be processed')
    parser.add_argument('-v', '--vcf_file', action='store', dest='vcf_file', required=True, help='REQUIRED: In file to be processed')
    parser.add_argument('-q', '--qual_threshold', action='store', dest='qual_threshold', default=150, help='OPTIONAL: QUAL values below this value will not be applied.')
    parser.add_argument('-m', '--map_threshold', action='store', dest='map_threshold', default=56, help='OPTIONAL: Map Quality values below this value will not be applied.')
    parser.add_argument('-a', '--ambiguity_NOT', action='store_true', dest='ambiguity_NOT', help='OPTIONAL: when -a used ambiguity will not be applied.  By default it will be applied.')

    args = parser.parse_args()
    print("\nSET ARGUMENTS: ")
    print(args)
    print("\n")

    merge_vcf = Merge_VCF_to_FASTA(args.qual_threshold, args.map_threshold, args.ambiguity_NOT)
    merge_vcf.merge_vcf_into_fasta(args.fasta_file, args.vcf_file)

# Created February 2020 by Tod Stuber
