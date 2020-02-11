#!/usr/bin/env python

__version__ = "2.0.0"

import os
import gzip
import re
import regex
import argparse
import textwrap
from collections import OrderedDict
import multiprocessing
multiprocessing.set_start_method('spawn', True)
from concurrent import futures
from dask import delayed
from Bio.SeqIO.QualityIO import FastqGeneralIterator


class Spoligo:

    def __init__(self, read1, read2, debug=None):
        self.debug = debug
        self.cpu_count_half = int(multiprocessing.cpu_count() / 2)
        self.fastq_list = []
        for read in [read1, read2]:
            if read: #append if not None
                self.fastq_list.append(read)
        self.paired = False
        if read2:
            self.paired = True
        real_path = os.path.dirname(os.path.realpath(__file__))
        self.spoligo_db = real_path + "/../dependencies/spoligotype_db.txt" 
        self.sample_name = re.sub('[_.].*', '', os.path.basename(read1))
        spoligo_dictionary = {}
        spoligo_dictionary["spacer01"] = ["TGATCCAGAGCCGGCGACCCTCTAT", "ATAGAGGGTCGCCGGCTCTGGATCA"]
        spoligo_dictionary["spacer02"] = ["CAAAAGCTGTCGCCCAAGCATGAGG", "CCTCATGCTTGGGCGACAGCTTTTG"]
        spoligo_dictionary["spacer03"] = ["CCGTGCTTCCAGTGATCGCCTTCTA", "TAGAAGGCGATCACTGGAAGCACGG"]
        spoligo_dictionary["spacer04"] = ["ACGTCATACGCCGACCAATCATCAG", "CTGATGATTGGTCGGCGTATGACGT"]
        spoligo_dictionary["spacer05"] = ["TTTTCTGACCACTTGTGCGGGATTA", "TAATCCCGCACAAGTGGTCAGAAAA"]
        spoligo_dictionary["spacer06"] = ["CGTCGTCATTTCCGGCTTCAATTTC", "GAAATTGAAGCCGGAAATGACGACG"]
        spoligo_dictionary["spacer07"] = ["GAGGAGAGCGAGTACTCGGGGCTGC", "GCAGCCCCGAGTACTCGCTCTCCTC"]
        spoligo_dictionary["spacer08"] = ["CGTGAAACCGCCCCCAGCCTCGCCG", "CGGCGAGGCTGGGGGCGGTTTCACG"]
        spoligo_dictionary["spacer09"] = ["ACTCGGAATCCCATGTGCTGACAGC", "GCTGTCAGCACATGGGATTCCGAGT"]
        spoligo_dictionary["spacer10"] = ["TCGACACCCGCTCTAGTTGACTTCC", "GGAAGTCAACTAGAGCGGGTGTCGA"]
        spoligo_dictionary["spacer11"] = ["GTGAGCAACGGCGGCGGCAACCTGG", "CCAGGTTGCCGCCGCCGTTGCTCAC"]
        spoligo_dictionary["spacer12"] = ["ATATCTGCTGCCCGCCCGGGGAGAT", "ATCTCCCCGGGCGGGCAGCAGATAT"]
        spoligo_dictionary["spacer13"] = ["GACCATCATTGCCATTCCCTCTCCC", "GGGAGAGGGAATGGCAATGATGGTC"]
        spoligo_dictionary["spacer14"] = ["GGTGTGATGCGGATGGTCGGCTCGG", "CCGAGCCGACCATCCGCATCACACC"]
        spoligo_dictionary["spacer15"] = ["CTTGAATAACGCGCAGTGAATTTCG", "CGAAATTCACTGCGCGTTATTCAAG"]
        spoligo_dictionary["spacer16"] = ["CGAGTTCCCGTCAGCGTCGTAAATC", "GATTTACGACGCTGACGGGAACTCG"]
        spoligo_dictionary["spacer17"] = ["GCGCCGGCCCGCGCGGATGACTCCG", "CGGAGTCATCCGCGCGGGCCGGCGC"]
        spoligo_dictionary["spacer18"] = ["CATGGACCCGGGCGAGCTGCAGATG", "CATCTGCAGCTCGCCCGGGTCCATG"]
        spoligo_dictionary["spacer19"] = ["TAACTGGCTTGGCGCTGATCCTGGT", "ACCAGGATCAGCGCCAAGCCAGTTA"]
        spoligo_dictionary["spacer20"] = ["TTGACCTCGCCAGGAGAGAAGATCA", "TGATCTTCTCTCCTGGCGAGGTCAA"]
        spoligo_dictionary["spacer21"] = ["TCGATGTCGATGTCCCAATCGTCGA", "TCGACGATTGGGACATCGACATCGA"]
        spoligo_dictionary["spacer22"] = ["ACCGCAGACGGCACGATTGAGACAA", "TTGTCTCAATCGTGCCGTCTGCGGT"]
        spoligo_dictionary["spacer23"] = ["AGCATCGCTGATGCGGTCCAGCTCG", "CGAGCTGGACCGCATCAGCGATGCT"]
        spoligo_dictionary["spacer24"] = ["CCGCCTGCTGGGTGAGACGTGCTCG", "CGAGCACGTCTCACCCAGCAGGCGG"]
        spoligo_dictionary["spacer25"] = ["GATCAGCGACCACCGCACCCTGTCA", "TGACAGGGTGCGGTGGTCGCTGATC"]
        spoligo_dictionary["spacer26"] = ["CTTCAGCACCACCATCATCCGGCGC", "GCGCCGGATGATGGTGGTGCTGAAG"]
        spoligo_dictionary["spacer27"] = ["GGATTCGTGATCTCTTCCCGCGGAT", "ATCCGCGGGAAGAGATCACGAATCC"]
        spoligo_dictionary["spacer28"] = ["TGCCCCGGCGTTTAGCGATCACAAC", "GTTGTGATCGCTAAACGCCGGGGCA"]
        spoligo_dictionary["spacer29"] = ["AAATACAGGCTCCACGACACGACCA", "TGGTCGTGTCGTGGAGCCTGTATTT"]
        spoligo_dictionary["spacer30"] = ["GGTTGCCCCGCGCCCTTTTCCAGCC", "GGCTGGAAAAGGGCGCGGGGCAACC"]
        spoligo_dictionary["spacer31"] = ["TCAGACAGGTTCGCGTCGATCAAGT", "ACTTGATCGACGCGAACCTGTCTGA"]
        spoligo_dictionary["spacer32"] = ["GACCAAATAGGTATCGGCGTGTTCA", "TGAACACGCCGATACCTATTTGGTC"]
        spoligo_dictionary["spacer33"] = ["GACATGACGGCGGTGCCGCACTTGA", "TCAAGTGCGGCACCGCCGTCATGTC"]
        spoligo_dictionary["spacer34"] = ["AAGTCACCTCGCCCACACCGTCGAA", "TTCGACGGTGTGGGCGAGGTGACTT"]
        spoligo_dictionary["spacer35"] = ["TCCGTACGCTCGAAACGCTTCCAAC", "GTTGGAAGCGTTTCGAGCGTACGGA"]
        spoligo_dictionary["spacer36"] = ["CGAAATCCAGCACCACATCCGCAGC", "GCTGCGGATGTGGTGCTGGATTTCG"]
        spoligo_dictionary["spacer37"] = ["CGCGAACTCGTCCACAGTCCCCCTT", "AAGGGGGACTGTGGACGAGTTCGCG"]
        spoligo_dictionary["spacer38"] = ["CGTGGATGGCGGATGCGTTGTGCGC", "GCGCACAACGCATCCGCCATCCACG"]
        spoligo_dictionary["spacer39"] = ["GACGATGGCCAGTAAATCGGCGTGG", "CCACGCCGATTTACTGGCCATCGTC"]
        spoligo_dictionary["spacer40"] = ["CGCCATCTGTGCCTCATACAGGTCC", "GGACCTGTATGAGGCACAGATGGCG"]
        spoligo_dictionary["spacer41"] = ["GGAGCTTTCCGGCTTCTATCAGGTA", "TACCTGATAGAAGCCGGAAAGCTCC"]
        spoligo_dictionary["spacer42"] = ["ATGGTGGGACATGGACGAGCGCGAC", "GTCGCGCTCGTCCATGTCCCACCAT"]
        spoligo_dictionary["spacer43"] = ["CGCAGAATCGCACCGGGTGCGGGAG", "CTCCCGCACCCGGTGCGATTCTGCG"]
        self.spoligo_dictionary = spoligo_dictionary

    def finding_sp(self, spacer_sequence):
        # spacer_id, spacer_sequence = spacer_id_and_spacer_sequence
        total_count = 0
        total_finds = 0
        #if total < 6: # doesn't make a big different.  Might as well get full counts
        #total += sum(seq.count(x) for x in (v)) #v=list of for and rev spacer
        total_finds = [len(regex.findall("(" + spacer + "){s<=1}", self.seq_string)) for spacer in spacer_sequence]
        for number in total_finds:
            total_count += number
        return (total_count)

    def binary_to_octal(self, binary):
        #binary_len = len(binary)
        i = 0
        ie = 1
        octal = ""
        while ie < 43:
            ie = i + 3
            # print(binary[i:ie])
            region = binary[i:ie]
            region_len = len(region)
            i += 3
            if int(region[0]) == 1:
                if region_len < 2: # for the lone spacer 43.  When present needs to be 1 not 4.
                    oct = 1
                else:
                    oct = 4
            else:
                oct = 0
            try:
                if int(region[1]) == 1:
                    oct += 2
                if int(region[2]) == 1:
                    oct += 1
            except IndexError:
                pass
            octal = octal + str(oct)
        return(octal)

    def spoligo(self):

        octal = None
        sbcode = None
        db_binarycode = None
        sample_binary = None

        seq_string = ""
        count_summary = {}
        sequence_list = []
        try:
            for fastq in self.fastq_list:
                with gzip.open(fastq, "rt") as in_handle:
                    # all 3, title and seq and qual, were needed
                    for title, seq, qual in FastqGeneralIterator(in_handle):
                        sequence_list.append(seq)
        except TypeError:
            # TypeError if not paired
            pass

        if len(seq) > 99:
            #Three 10bp sequences dispersed across repeat region, forward and reverse
            capture_spacer_sequence = re.compile(".*TTTCCGTCCC.*|.*GGGACGGAAA.*|.*TCTCGGGGTT.*|.*AACCCCGAGA.*|.*TGGGTCTGAC.*|.*GTCAGACCCA.*")
            sequence_list = list(filter(capture_spacer_sequence.match, sequence_list))
            seq_string = "".join(sequence_list)
        else:
            #if < 100 then search all reads, not just those with repeat regions.
            seq_string = "".join(sequence_list)
        self.seq_string = seq_string
        finding_sp = self.finding_sp
        # if self.debug:
        #     for spacer_id_and_spacer_sequence in self.spoligo_dictionary.items():
        #         total_count, spacer_id = finding_sp(spacer_id_and_spacer_sequence)
        #         count_summary.update({spacer_id: total_count})
        # else:
        #     spoligo_dictionary = self.spoligo_dictionary
        #     cpu_count_half = self.cpu_count_half
        #     with futures.ProcessPoolExecutor(max_workers=cpu_count_half) as pool:
        #         for total_count, spacer_id in pool.map(finding_sp, spoligo_dictionary.items()):
        #             count_summary.update({spacer_id: total_count})
        for spacer_id, spacer_sequence in self.spoligo_dictionary.items():
            count = delayed(self.finding_sp)(spacer_sequence)
            count_summary.update({spacer_id: count})
        pull = delayed(count_summary)
        count_summary = pull.compute()
        count_summary = OrderedDict(sorted(count_summary.items()))
        spoligo_binary_dictionary = {}
        for k, v in count_summary.items():
            if v > 4:
                spoligo_binary_dictionary.update({k: 1})
            else:
                spoligo_binary_dictionary.update({k: 0})
        spoligo_binary_dictionary = OrderedDict(sorted(spoligo_binary_dictionary.items()))
        spoligo_binary_list = []
        for v in spoligo_binary_dictionary.values():
            spoligo_binary_list.append(v)
        sample_binary = ''.join(str(e) for e in spoligo_binary_list)  #sample_binary correct
        self.sample_binary = sample_binary
        self.octal = self.binary_to_octal(sample_binary)
        write_out = open("spoligo.txt", 'w')
        found = False
        with open(self.spoligo_db) as spoligo_db_file: # put into dictionary or list
            for line in spoligo_db_file:
                line = line.rstrip()
                sbcode = line.split()[1]
                db_binarycode = line.split()[2]
                if sample_binary == db_binarycode:
                    found = True
                    self.sbcode = sbcode
        if sample_binary == '0000000000000000000000000000000000000000000':
            self.sbcode = "spoligo not found, binary all zeros, see spoligo file"
            # print("CHECK SAMPLE!  NO SPACERS FOUND.  LIKELY NOT TB COMPLEX.  ALTHOUGH SB2277 IS A ZERO STRING BINARY\n")
            print("CHECK SAMPLE!  NO SPACERS FOUND.  LIKELY NOT TB COMPLEX.  ALTHOUGH SB2277 IS A ZERO STRING BINARY", file=write_out)
            print("\nOne mismatch allowed spacers search against both R1 and R2 reads.\n", file=write_out)
            for k, v in count_summary.items():
                print(k, v, file=write_out)
        elif found:
            print(f'{self.octal}, {self.sbcode}, {sample_binary}', file=write_out)
            print("One mismatch allowed spacer search against both R1 and R2 reads.\n", file=write_out)
            for k, v in count_summary.items():
                print(k, v, file=write_out)
            print(f'sample_binary: {sample_binary}', file=write_out)
        else:
            self.sbcode = "Not Found"
            print(f'{self.octal}, {self.sbcode}, {sample_binary}', file=write_out)
            print("\nSPOLIGO SB NUMBER NOT FOUND\n", file=write_out)
            print("\nOne mismatch allowed spacer search against both R1 and R2 reads.\n", file=write_out)
            for k, v in count_summary.items():
                print(k, v, file=write_out)

if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
   
    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-r1', '--read1', action='store', dest='read1', required=True, help='Required: single read')
    parser.add_argument('-r2', '--read2', action='store', dest='read2', required=False, default=None, help='Optional: paired read')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='turn off map.pooling of samples')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    read1 = args.read1
    read2 = args.read2

    spoligo = Spoligo(read1, read2, args.debug)
    spoligo.spoligo()
    
    print("Done")