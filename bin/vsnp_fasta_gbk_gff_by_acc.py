#!/usr/bin/env python

__version__ = "2.03"

import os
import argparse
import textwrap
from Bio import SeqIO
from Bio import Entrez

class Downloader:

    def __init__(self, accession):
        self.entrezDbName = 'nucleotide'
        self.email = 'mickey_mouse@gmail.com'
        self.accession = accession

    def gbk(self):
        Entrez.email = self.email
        print(f"Downloading {self.accession} gbk")
        entryData = Entrez.efetch(db=self.entrezDbName, id=self.accession, retmode="text", rettype='gb')
        writeFile = self.accession + ".gbk"
        local_file=open(writeFile,"w")
        local_file.write(entryData.read())
        entryData.close()
        local_file.close()

    def gff(self):
        Entrez.email = self.email
        print(f"Downloading {self.accession} gff3")
        entryData = Entrez.efetch(db=self.entrezDbName, id=self.accession, retmode="text", rettype='gff3')
        writeFile = self.accession + ".gff"
        local_file=open(writeFile,"w")
        local_file.write(entryData.read())
        entryData.close()
        local_file.close()

    def fasta(self):
        Entrez.email = self.email
        print(f"Downloading {self.accession} FASTA")
        entryData = Entrez.efetch(db=self.entrezDbName, id=self.accession, retmode="text", rettype='fasta')
        writeFile = self.accession + ".fasta"
        local_file=open(writeFile,"w")
        local_file.write(entryData.read())
        entryData.close()
        local_file.close()

        handle = open(writeFile, "r")
        for record in SeqIO.parse(handle, "fasta"):
            print(f"{record.description}")
            print(f'Sequence length: {len(record.seq):,}\n')
        handle.close()
        return record.description
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------

    Usage:
        fasta_gbk_gff_by_acc.py -a NC_002945 -f
        fasta_gbk_gff_by_acc.py -a NC_006932 -fg
        fasta_gbk_gff_by_acc.py -a NC_006933 -fbg
        fasta_gbk_gff_by_acc.py -a CP023243 -fbg
        fasta_gbk_gff_by_acc.py -a NZ_CP023243 -fbg
        **bad request: fasta_gbk_gff_by_acc.py -a NZ_AQME0000000 -fbg # must be complete chromosome

    vSNP requires multi-chromosome genomes to be concatenated to single file

    Search genomes: https://www.ncbi.nlm.nih.gov/genome

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-a', '--accession', action='store', dest='accession', help='NCBI chromosome number')
    parser.add_argument('-f', '--fasta', action='store_true', dest='fasta', help='get FASTA file')
    parser.add_argument('-b', '--gbk', action='store_true', dest='gbk', help='get gbk file')
    parser.add_argument('-g', '--gff', action='store_true', dest='gff', help='get gff file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    download = Downloader(args.accession)
    if args.fasta:
        download.fasta()
    if args.gbk:
        download.gbk()
    if args.gff:
        download.gff()
    
