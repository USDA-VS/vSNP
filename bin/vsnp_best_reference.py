#!/usr/bin/env python

__version__ = "2.03"

import os
import re
import gzip
import argparse
import textwrap
from collections import OrderedDict
from concurrent import futures
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import multiprocessing
multiprocessing.set_start_method('spawn', True)


class Best_Reference:

    def __init__(self, read1, read2):
        self.fastq_list = []
        for read in [read1, read2]:
            if read: #append if not None
                self.fastq_list.append(read)
        self.paired = False
        if read2:
            self.paired = True
        self.sample_name = re.sub('[_.].*', '', os.path.basename(read1))

        oligo_dictionary = {}
        oligo_dictionary["01_ab1"] = "AATTGTCGGATAGCCTGGCGATAACGACGC"
        oligo_dictionary["02_ab3"] = "CACACGCGGGCCGGAACTGCCGCAAATGAC"
        oligo_dictionary["03_ab5"] = "GCTGAAGCGGCAGACCGGCAGAACGAATAT"
        oligo_dictionary["04_mel"] = "TGTCGCGCGTCAAGCGGCGTGAAATCTCTG"
        oligo_dictionary["05_suis1"] = "TGCGTTGCCGTGAAGCTTAATTCGGCTGAT"
        oligo_dictionary["06_suis2"] = "GGCAATCATGCGCAGGGCTTTGCATTCGTC"
        oligo_dictionary["07_suis3"] = "CAAGGCAGATGCACATAATCCGGCGACCCG"
        oligo_dictionary["08_ceti1"] = "GTGAATATAGGGTGAATTGATCTTCAGCCG"
        oligo_dictionary["09_ceti2"] = "TTACAAGCAGGCCTATGAGCGCGGCGTGAA"
        oligo_dictionary["10_canis4"] = "CTGCTACATAAAGCACCCGGCGACCGAGTT"
        oligo_dictionary["11_canis"] = "ATCGTTTTGCGGCATATCGCTGACCACAGC"
        oligo_dictionary["12_ovis"] = "CACTCAATCTTCTCTACGGGCGTGGTATCC"
        oligo_dictionary["13_ether2"] = "CGAAATCGTGGTGAAGGACGGGACCGAACC"
        oligo_dictionary["14_63B1"] = "CCTGTTTAAAAGAATCGTCGGAACCGCTCT"
        oligo_dictionary["15_16M0"] = "TCCCGCCGCCATGCCGCCGAAAGTCGCCGT"
        oligo_dictionary["16_mel1b"] = "TCTGTCCAAACCCCGTGACCGAACAATAGA" #added 2018-01-30
        oligo_dictionary["17_tb157"] = "CTCTTCGTATACCGTTCCGTCGTCACCATGGTCCT"
        oligo_dictionary["18_tb7"] = "TCACGCAGCCAACGATATTCGTGTACCGCGACGGT"
        oligo_dictionary["19_tbbov"] = "CTGGGCGACCCGGCCGACCTGCACACCGCGCATCA"
        oligo_dictionary["20_tb5"] = "CCGTGGTGGCGTATCGGGCCCCTGGATCGCGCCCT"
        oligo_dictionary["21_tb2"] = "ATGTCTGCGTAAAGAAGTTCCATGTCCGGGAAGTA"
        oligo_dictionary["22_tb3"] = "GAAGACCTTGATGCCGATCTGGGTGTCGATCTTGA"
        oligo_dictionary["23_tb4"] = "CGGTGTTGAAGGGTCCCCCGTTCCAGAAGCCGGTG"
        oligo_dictionary["24_tb6"] = "ACGGTGATTCGGGTGGTCGACACCGATGGTTCAGA"
        oligo_dictionary["25_para"] = "CCTTTCTTGAAGGGTGTTCG"
        oligo_dictionary["26_para_sheep"] = "CGTGGTGGCGACGGCGGCGGGCCTGTCTAT"
        oligo_dictionary["27_para_cattle"] = "TCTCCTCGGTCGGTGATTCGGGGGCGCGGT"
        brucella_identifications = {}
        brucella_identifications["1111111111111111"] = "odd" #Unexpected findings
        brucella_identifications["0111111111111111"] = "Brucella_abortus1" #Brucella abortus bv 1, 2 or 4
        brucella_identifications["1011111111111111"] = "Brucella_abortus3" #Brucella abortus bv 3
        brucella_identifications["1101111111111111"] = "Brucella_abortus1" #Brucella abortus bv 5, 6 or 9
        brucella_identifications["1110111111111101"] = "Brucella_melitensis-bv1"
        brucella_identifications["0000010101101101"] = "Brucella_melitensis-bv1"
        brucella_identifications["1110111111111100"] = "Brucella_melitensis-bv1b" #added 2018-01-30
        brucella_identifications["0000010101101100"] = "Brucella_melitensis-bv1b" #added 2018-01-30
        brucella_identifications["1110111111111011"] = "Brucella_melitensis-bv2"
        brucella_identifications["0000010101101001"] = "Brucella_melitensis-bv2"
        brucella_identifications["0100010101101001"] = "Brucella_melitensis-bv2"
        brucella_identifications["1110011111101011"] = "Brucella_melitensis-bv2"
        brucella_identifications["1110111111110111"] = "Brucella_melitensis-bv3"
        brucella_identifications["1110011111100111"] = "Brucella_melitensis-bv3"
        brucella_identifications["1111011111111111"] = "Brucella_suis1"
        brucella_identifications["1111101111111111"] = "Brucella_suis2"
        brucella_identifications["1111110111111101"] = "Brucella_suis3"
        brucella_identifications["1111111011111111"] = "Brucella_ceti1"
        brucella_identifications["1111111001111111"] = "Brucella_ceti1"
        brucella_identifications["1111111101111111"] = "Brucella_ceti2"
        brucella_identifications["1111111110111101"] = "Brucella_suis4"
        brucella_identifications["1111111110011101"] = "Brucella_canis"
        brucella_identifications["1111111111101111"] = "Brucella_ovis"
        bovis_identifications = {}
        bovis_identifications["11101111"] = "Mycobacterium_H37" #tb1
        bovis_identifications["11101101"] = "Mycobacterium_H37" #tb1
        bovis_identifications["01100111"] = "Mycobacterium_H37" #tb2
        bovis_identifications["01101011"] = "Mycobacterium_H37" #tb3
        bovis_identifications["11101011"] = "Mycobacterium_H37" #tb3
        bovis_identifications["01101111"] = "Mycobacterium_H37" #tb4a
        bovis_identifications["01101101"] = "Mycobacterium_H37" #tb4b
        bovis_identifications["11101101"] = "Mycobacterium_H37" #tb4b
        bovis_identifications["01101111"] = "Mycobacterium_H37" #tb4b
        bovis_identifications["11111111"] = "Mycobacterium_H37" #tb5
        bovis_identifications["11001111"] = "Mycobacterium_H37" #tb6
        bovis_identifications["10101110"] = "Mycobacterium_H37" #tb7
        bovis_identifications["11001110"] = "Mycobacterium_AF2122" #bovis
        bovis_identifications["11011110"] = "Mycobacterium_AF2122" #bovis
        bovis_identifications["11001100"] = "Mycobacterium_AF2122"  #bovis
        para_identifications = {}
        para_identifications["110"] = "para-CP033688"
        para_identifications["101"] = "para-NC002944"

        self.oligo_dictionary = oligo_dictionary
        self.brucella_identifications = brucella_identifications
        self.bovis_identifications = bovis_identifications
        self.para_identifications = para_identifications

    def best_reference(self):

        write_out = open("best_reference.txt", 'w')
        count_summary = {}

        for value in self.oligo_dictionary.values():
            returned_value, count = self.finding_best_ref(value)
            for key, value in self.oligo_dictionary.items():
                if returned_value == value:
                    count_summary.update({key: count})
                    count_summary = OrderedDict(sorted(count_summary.items()))

        # https://github.com/Microsoft/ptvsd/issues/1165 for Error in atexit._run_exitfuncs: OSError: handle is closed in VS Code
        # with futures.ProcessPoolExecutor() as pool:  #ProcessPoolExecutor ThreadPoolExecutor
        #     for returned_value, count in pool.map(self.finding_best_ref, self.oligo_dictionary.values()):
        #         for key, value in self.oligo_dictionary.items():
        #             if returned_value == value:
        #                 count_summary.update({key: count})
        #                 count_summary = OrderedDict(sorted(count_summary.items()))        
        count_list = []
        for v in count_summary.values():
            count_list.append(v)
        brucella_sum = sum(count_list[:16])
        bovis_sum = sum(count_list[16:24])
        para_sum = sum(count_list[24:])
        print("Best reference Brucella counts:", file=write_out)
        for i in count_list[:16]:
            print(i, end=',', file=write_out)
        print("\nBest reference TB counts:", file=write_out)
        for i in count_list[16:24]:
            print(i, end=',', file=write_out)

        print("\nBest reference Para counts:", file=write_out)
        for i in count_list[24:]:
            print(i, end=',', file=write_out)

        #Binary dictionary
        binary_dictionary = {}
        for k, v in count_summary.items():
            if v > 1:
                binary_dictionary.update({k: 1})
            else:
                binary_dictionary.update({k: 0})
        binary_dictionary = OrderedDict(sorted(binary_dictionary.items()))

        binary_list = []
        for v in binary_dictionary.values():
            binary_list.append(v)
        brucella_binary = binary_list[:16]
        brucella_string = ''.join(str(e) for e in brucella_binary)
        bovis_binary = binary_list[16:24]
        bovis_string = ''.join(str(e) for e in bovis_binary)
        para_binary = binary_list[24:]
        para_string = ''.join(str(e) for e in para_binary)

        if brucella_sum > 3:
            self.group = "Brucella"
            if brucella_string in self.brucella_identifications:
                self.species = self.brucella_identifications[brucella_string]
                print("\n\nBrucella group, species %s" % self.brucella_identifications[brucella_string], file=write_out)
            else:
                print("\n\nBrucella group, but no match", file=write_out)
                self.species = "No Findings"
        elif bovis_sum > 3:
            self.group = "TB"
            if bovis_string in self.bovis_identifications:
                print("\n\nTB group, species %s" % self.bovis_identifications[bovis_string], file=write_out)
                self.species = self.bovis_identifications[bovis_string]
            else:
                print("\n\nTB group, but no match", file=write_out)
                self.species = "No Findings"
        elif para_sum >= 1:
            self.group = "paraTB"
            if para_string in self.para_identifications:
                print("\n\nPara group, species %s" % self.para_identifications[para_string], file=write_out)
                self.species = self.para_identifications[para_string]
            else:
                print("\n\nM. paratuberculosis group, but no match", file=write_out)
                self.species = "No Findings"
        else:
            self.group = "No Findings"
            self.species = "No Findings"
            print("Unable to find a best reference species or group")
        write_out.close()

    def finding_best_ref(self, value):
        count = 0
        for fastq in self.fastq_list:
            with gzip.open(fastq, 'rt') as in_handle:
                # all 3, title and seq and qual, were needed
                for title, seq, qual in FastqGeneralIterator(in_handle):
                    count += seq.count(value)
        return(value, count)

if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Use oligos to determine species.  Most often if the absents of a single oligo from a set specific for either brucella or bovis will confer species type.  Some species will the absents of more than one oligo.  Oligo findings are translated to binary patterns.
   
    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-r1', '--read1', action='store', dest='read1', required=True, help='Required: single read')
    parser.add_argument('-r2', '--read2', action='store', dest='read2', required=False, default=None, help='Optional: paired read')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    read1 = args.read1
    read2 = args.read2

    best_ref = Best_Reference(read1, read2)
    best_ref.best_reference()
    print(f'### {best_ref.group}, {best_ref.species}')