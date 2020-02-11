#!/usr/bin/env python

__version__ = "2.0.02"

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

from vsnp_fastq_quality import FASTQ_Quality
from vsnp_spoligotype import Spoligo
from vsnp_bruc_mlst import Bruc_MLST
from vsnp_chromosome_reference import Reference_Chromosome
from vsnp_group_reporter import GroupReporter

def print_options():
    script_path = os.path.dirname(os.path.realpath(__file__))
    ref_options = glob.glob(f'{script_path}/dependencies/*')
    print(f'\nReference Options:')
    for option in ref_options:
        if os.path.isdir(option):
            print(f'\t{os.path.basename(option)}')
    print("")

class Align_Reads:

    def __init__(self, read1, read2, reference, gbk, species=None, group=None, skip_assembly=False):
        self.startTime = datetime.now()
        start_time = datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
        print("\n########################")
        print(f'Start time: {start_time}')
        self.fastq_list = []
        for read in [read1, read2]:
            if read: #append if not None
                self.fastq_list.append(read)
        self.paired = False
        self.read1 = self.fastq_list[0]
        self.read2 = None
        if read2:
            self.paired = True
            self.read2 = self.fastq_list[1]
        self.reference = reference
        self.sample_name = re.sub('[_.].*', '', os.path.basename(read1))
        self.root_dir = str(os.getcwd())
        self.gbk = gbk
        self.species = species
        self.group = group
        self.skip_assembly = skip_assembly

    def add_zero_coverage(self, sample_name, sample_reference, nodupbam, hapall, zero_coverage_vcf):
        coverage_dict = {}
        coverage_list = pysam.depth(nodupbam, split_lines=True)
        for line in coverage_list:
            chrom, position, depth = line.split('\t')
            coverage_dict[chrom + "-" + position] = depth
        coverage_df = pd.DataFrame.from_dict(coverage_dict, orient='index', columns=["depth"])
        zero_dict = {}
        for record in SeqIO.parse(sample_reference, "fasta"):
            chrom = record.id
            total_len = len(record.seq)
            for pos in list(range(1, total_len + 1)):
                zero_dict[str(chrom) + "-" + str(pos)] = 0
        zero_df = pd.DataFrame.from_dict(zero_dict, orient='index', columns=["depth"])
        #df with depth_x and depth_y columns, depth_y index is NaN
        coverage_df = zero_df.merge(coverage_df, left_index=True, right_index=True, how='outer')
        #depth_x "0" column no longer needed
        coverage_df = coverage_df.drop(columns=['depth_x'])
        coverage_df = coverage_df.rename(columns={'depth_y': 'depth'})
        #covert the NaN to 0 coverage
        coverage_df = coverage_df.fillna(0)
        coverage_df['depth'] = coverage_df['depth'].apply(int)
        total_length = len(coverage_df)
        ave_coverage = coverage_df['depth'].mean()
        zero_df = coverage_df[coverage_df['depth'] == 0]
        total_zero_coverage = len(zero_df)
        print("\tPositions with no coverage: {:,}" .format(total_zero_coverage))
        total_coverage = total_length - total_zero_coverage
        genome_coverage = "{:.2%}".format(total_coverage / total_length)
        vcf_df = pd.read_csv(hapall, sep='\t', header=None, names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"], comment='#')
        good_snp_count = len(vcf_df[(vcf_df['ALT'].str.len() == 1) & (vcf_df['REF'].str.len() == 1) & (vcf_df['QUAL'] > 150)])
        if total_zero_coverage > 0:
            header_out = open('v_header.csv', 'w+')
            with open(hapall) as fff:
                for line in fff:
                    if re.search('^#', line):
                        print(line.strip(), file=header_out)
            header_out.close()
            vcf_df_snp = vcf_df[vcf_df['REF'].str.len() == 1]
            vcf_df_snp = vcf_df_snp[vcf_df_snp['ALT'].str.len() == 1]
            vcf_df_snp['ABS_VALUE'] = vcf_df_snp['CHROM'].map(str) + '-' + vcf_df_snp['POS'].map(str)
            vcf_df_snp = vcf_df_snp.set_index('ABS_VALUE')
            cat_df = pd.concat([vcf_df_snp, zero_df], axis=1, sort=False)
            cat_df = cat_df.drop(columns=['CHROM', 'POS', 'depth'])
            cat_df[['ID', 'ALT', 'QUAL', 'FILTER', 'INFO']] = cat_df[['ID', 'ALT', 'QUAL', 'FILTER', 'INFO']].fillna('.')
            cat_df['REF'] = cat_df['REF'].fillna('N')
            cat_df['FORMAT'] = cat_df['FORMAT'].fillna('GT')
            cat_df['Sample'] = cat_df['Sample'].fillna('./.')
            cat_df['temp'] = cat_df.index.str.rsplit('-', n=1)
            cat_df[['CHROM', 'POS']] = pd.DataFrame(cat_df.temp.values.tolist(), index=cat_df.index)
            cat_df = cat_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Sample']]
            cat_df['POS'] = cat_df['POS'].astype(int)
            cat_df = cat_df.sort_values(['CHROM', 'POS'])
            cat_df.to_csv('v_annotated_body.csv', sep='\t', header=False, index=False)
            cat_files = ['v_header.csv', 'v_annotated_body.csv']
            with open(zero_coverage_vcf, "wb") as outfile:
                for cf in cat_files:
                    with open(cf, "rb") as infile:
                        outfile.write(infile.read())
            os.remove('v_header.csv')
            os.remove('v_annotated_body.csv')
        else:
            shutil.copyfile(hapall, zero_coverage_vcf)
        return (zero_coverage_vcf, good_snp_count, ave_coverage, genome_coverage)

    def align(self):

        fq = FASTQ_Quality(self.read1, self.read2)
        print(f'\n########################')
        fq.get_quality()

        sample_name = self.sample_name
        fastq_list = self.fastq_list
        reference = self.reference
        print(f'### {sample_name} Making indexes...')
        os.system(f'samtools faidx {reference}')
        os.system(f'picard CreateSequenceDictionary REFERENCE={reference} OUTPUT={reference.rsplit(".", 1)[0]}.dict 2> /dev/null')
        os.system(f'bwa index {reference} 2> /dev/null')
        samfile = f'{sample_name}.sam'
        print(f'### {sample_name} BWA Aligning reads...')
        if self.paired:
            os.system(f'bwa mem -M -R "@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tPI:250" -t 8 {reference} {fastq_list[0]} {fastq_list[1]} > {samfile}')
        else:
            os.system(f'bwa mem -M -R "@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tPI:250" -t 8 {reference} {fastq_list[0]} > {samfile}')
        all_bamfile = f'{sample_name}_all.bam'
        os.system(f'samtools view -Sb {samfile} -o {all_bamfile}')
        sorted_bamfile = f'{sample_name}_sorted.bam'
        print(f'### {sample_name} Sorting BAM...')
        os.system(f'samtools sort {all_bamfile} -o {sorted_bamfile}')
        os.system(f'samtools index {sorted_bamfile}')
        rmdup_bamfile = f'{sample_name}.bam'
        print(f'### {sample_name} Removing Duplicates...')
        os.system(f'picard MarkDuplicates INPUT={sorted_bamfile} OUTPUT={rmdup_bamfile} ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=dup_metrics.csv 2> /dev/null')
        os.system(f'samtools index {rmdup_bamfile}')
        chrom_ranges = open("chrom_ranges.txt", 'w')
        for record in SeqIO.parse(reference, "fasta"):
            chrom = record.id
            total_len = len(record.seq)
            min_number = 0
            step = 100000
            if step < total_len:
                for chunk in range(min_number, total_len, step)[1:]:
                    print("{}:{}-{}".format(chrom, min_number, chunk), file=chrom_ranges)
                    min_number = chunk
            print("{}:{}-{}".format(chrom, min_number, total_len), file=chrom_ranges)
        chrom_ranges.close()
        unfiltered_hapall = f'{sample_name}_unfiltered_hapall.vcf'
        print(f'### {sample_name} Calling SNPs with Freebayes...')
        os.system(f'freebayes-parallel chrom_ranges.txt 8 -E -1 -e 1 -u --strict-vcf -f {reference} {rmdup_bamfile} > {unfiltered_hapall}')
        # "fix" MQ notation in VCF to match GATK output
        mapfix_hapall = f'{sample_name}_mapfix_hapall.vcf'
        write_fix = open(mapfix_hapall, 'w+')
        with open(unfiltered_hapall, 'r') as unfiltered:
            for line in unfiltered:
                line = line.strip()
                new_line = re.sub(r';MQM=', r';MQ=', line)
                new_line = re.sub(r'ID=MQM,', r'ID=MQ,', new_line)
                print(new_line, file=write_fix)
            write_fix.close()
        # remove clearly poor positions
        filtered_hapall = f'{sample_name}_filtered_hapall.vcf'
        os.system(f'vcffilter -f "QUAL > 20" {mapfix_hapall} > {filtered_hapall}')
        if self.skip_assembly:
            print(f'Skipping assembly of unmapped reads')
            self.abyss_contig_count = 'Assembly of unmapped reads not done'
        else:
            print(f'### {sample_name} Assembling Unmapped Reads...')
            abyss_contig_count = 0
            if self.paired:
                unaligned_read1 = f'{sample_name}_R1_unaligned.fastq'
                unaligned_read2 = f'{sample_name}_R2_unaligned.fastq'
                os.system(f'samtools fastq -f4 -1 {unaligned_read1} -2 {unaligned_read2} --reference {reference} --threads 8 {rmdup_bamfile} 2> /dev/null' )
                abyss_out = f'{sample_name}_unaligned_contigs.fasta'
                try:
                    os.system(f'ABYSS --out {abyss_out} --coverage 5 --kmer 64 {unaligned_read1} {unaligned_read2} &> /dev/null')
                    with open(abyss_out) as f:
                        for line in f:
                            abyss_contig_count += line.count(">")
                except FileNotFoundError:
                    abyss_contig_count = 0
            else:
                unaligned_read = f'{sample_name}_unaligned.fastq'
                os.system(f'samtools fastq -f4 -0 {unaligned_read} --reference {reference} --threads 8 {rmdup_bamfile} 2> /dev/null' )
                abyss_out = f'{sample_name}_unaligned_contigs.fasta'
                try:
                    os.system(f'ABYSS --out {abyss_out} --coverage 5 --kmer 64 {unaligned_read} &> /dev/null')
                    with open(abyss_out) as f:
                        for line in f:
                            abyss_contig_count += line.count(">")
                except FileNotFoundError:
                    abyss_contig_count = 0
            self.abyss_contig_count = abyss_contig_count
        os.makedirs('unmapped_reads')
        unaligned_files = glob.glob('*unaligned*')
        for file in unaligned_files:
            newZip = zipfile.ZipFile(f'{file}.gz', 'w')
            newZip.write(file, compress_type=zipfile.ZIP_DEFLATED)
            newZip.close()
            os.remove(file)
            shutil.move(f'{file}.gz', 'unmapped_reads')

        # Full bam stats
        stat_out = open("stat_align.txt", 'w')
        stat_out.write(os.popen("samtools idxstats {} " .format(rmdup_bamfile)).read())
        stat_out.close()
        with open("stat_align.txt", 'r') as f:
            first_line = f.readline()
            first_line = first_line.rstrip()
            first_line = re.split(':|\t', first_line)
            self.allbam_mapped_reads = int(first_line[2])
        with open("stat_align.txt") as f:
            for line_num, line in enumerate(f):
                if line_num == 1:
                    dup_line_two = line
                    dup_line_two = dup_line_two.split()
                    self.unmapped_reads = int(dup_line_two[3])
        
        os.remove(f'{reference.rsplit(".", 1)[0]}.dict')
        os.remove(f'{reference}.pac')
        os.remove(f'{reference}.bwt')
        os.remove(f'{reference}.sa')
        os.remove(f'{reference}.ann')
        os.remove(f'{reference}.amb')
        os.remove(f'{sample_name}.sam')
        os.remove(f'{sample_name}_all.bam')
        os.remove(f'{sample_name}_sorted.bam')
        os.remove(f'{sample_name}_sorted.bam.bai')
        os.remove(f'dup_metrics.csv')
        os.remove(f'chrom_ranges.txt')
        os.remove(f'{sample_name}_unfiltered_hapall.vcf')
        os.remove(f'{sample_name}_mapfix_hapall.vcf')
        os.remove(f'stat_align.txt')

        zero_coverage_vcf = f'{sample_name}_zc.vcf'
        print(f'### {sample_name} Adding Zero Coverage...')
        zero_coverage_vcf, self.good_snp_count, self.ave_coverage, self.genome_coverage = self.add_zero_coverage(sample_name, reference, rmdup_bamfile, filtered_hapall, zero_coverage_vcf)

        if self.gbk:
            print(f'### {sample_name} Annotating VCF file...')
            annotated_vcf = f'{sample_name}_annotated.vcf'
            gbk_file = self.gbk
            annotation_dict = {}
            for gbk in gbk_file:
                print("gbk file: %s" % gbk)
                write_out = open('temp.csv', 'w+')
                gbk_dict = SeqIO.to_dict(SeqIO.parse(gbk, "genbank"))
                gbk_chrome = list(gbk_dict.keys())[0]
                for key, value in gbk_dict.items():
                    for feature in value.features:
                        if "CDS" in feature.type or "rRNA" in feature.type:
                            myproduct = None
                            mylocus = None
                            mygene = None
                            try:
                                myproduct = feature.qualifiers['product'][0]
                            except KeyError:
                                pass
                            try:
                                mylocus = feature.qualifiers['locus_tag'][0]
                            except KeyError:
                                pass
                            try:
                                mygene = feature.qualifiers['gene'][0]
                            except KeyError:
                                pass
                            print(key, int(feature.location.start), int(feature.location.end), mylocus, myproduct, mygene, sep='\t', file=write_out)
                write_out.close()
                df = pd.read_csv('temp.csv', sep='\t', names=["chrom", "start", "stop", "locus", "product", "gene"])
                df = df.sort_values(['start', 'gene'], ascending=[True, False])
                df = df.drop_duplicates('start')
                pro = df.reset_index(drop=True)
                pro.index = pd.IntervalIndex.from_arrays(pro['start'], pro['stop'], closed='both')
                annotation_dict[gbk_chrome] = pro

            header_out = open('v_header.csv', 'w+')
            with open(zero_coverage_vcf) as fff:
                for line in fff:
                    if re.search('^#', line):
                        print(line.strip(), file=header_out)
            header_out.close()

            vcf_df = pd.read_csv(zero_coverage_vcf, sep='\t', header=None, names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"], comment='#')
            vcf_df['ABS_VALUE'] = vcf_df['CHROM'].map(str) + '-' + vcf_df['POS'].map(str)
            vcf_df = vcf_df.set_index('ABS_VALUE')

            annotate_condense_dict = {}
            for gbk_chrome, pro in annotation_dict.items():
                matching_chrom_df = vcf_df[vcf_df['CHROM'] == gbk_chrome]
                for index, row in matching_chrom_df.iterrows():
                    pos = row.POS
                    try:
                        aaa = pro.iloc[pro.index.get_loc(int(pos))][['chrom', 'locus', 'product', 'gene']]
                        try:
                            chrom, name, locus, tag = aaa.values[0]
                            annotate_condense_dict[str(chrom) + "-" + str(pos)] = "{}, {}, {}".format(name, locus, tag)
                        except ValueError:
                            # if only one annotation entire chromosome (such with flu) then having [0] fails
                            chrom, name, locus, tag = aaa.values
                            annotate_condense_dict[str(chrom) + "-" + str(pos)] = "{}, {}, {}".format(name, locus, tag)
                    except KeyError:
                        annotate_condense_dict[str(gbk_chrome) + "-" + str(pos)] = "No annotated product"

            annotate_df = pd.DataFrame.from_dict(annotate_condense_dict, orient='index', columns=["ID"])
            annotate_df.index.name = 'ABS_VALUE'
            vcf_df.drop(['ID'], axis=1, inplace=True)

            vcf_df = vcf_df.merge(annotate_df, how='left', left_index=True, right_index=True)
            vcf_df = vcf_df[["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"]]
            vcf_df.to_csv('v_annotated_body.csv', sep='\t', header=False, index=False)

            cat_files = ['v_header.csv', 'v_annotated_body.csv']
            with open(annotated_vcf, "wb") as outfile:
                for cf in cat_files:
                    with open(cf, "rb") as infile:
                        outfile.write(infile.read())

        if self.group == 'TB':
            print(f'Getting spoligo...')
            spoligo = Spoligo(self.read1, self.read2)
            spoligo.spoligo()
        if self.group == 'Brucella':
            print(f'Getting Brucella mlst...')
            mlst = Bruc_MLST(self.read1, self.read2)
            mlst.find_mlst()

        #Capture program versions for step 1
        try:
            verison_out = open("versions.txt", 'w')
            print(os.popen('conda list bwa | grep -v "^#"; \
                conda list abyss | grep -v "^#"; \
                conda list picard | grep -v "^#"; \
                conda list samtools | grep -v "^#"; \
                conda list freebayes | grep -v "^#"; \
                conda list biopython | grep -v "^#"').read(), file=verison_out)
            print(f'Dependent source: {self.sample_name} {self.fastq_list} {self.reference} {self.gbk }', file=verison_out)
            verison_out.close()
        except:
            print(f'##### {self.sample_name} DID NOT CAPTURE VERSION INFO')
            pass
           
        unaligned_files = glob.glob('*unaligned*')
        for file in unaligned_files:
            newZip = zipfile.ZipFile(f'{file}.gz', 'w')
            newZip.write(file, compress_type=zipfile.ZIP_DEFLATED)
            newZip.close()
            os.remove(file)
            shutil.move(f'{file}.gz', 'unmapped_reads')
            
        reference_type = Reference_Chromosome(self.root_dir)
        ref_option = reference_type.get_reference()
        print(f'Reference type found: {ref_option}\n')
        if ref_option:
            group_reporter = GroupReporter(zero_coverage_vcf, ref_option)
            group_list = group_reporter.get_groups()
            group_string = ", ".join(str(x) for x in group_list)
        else:
            group_string = 'Unable to cross-reference chrom, check vsnp_chromosome_reference.py'

        os.makedirs('alignment')
        files_grabbed = [glob.glob(e) for e in ['*bam', '*bam.bai', '*vcf', '*fasta', '*fasta.fai']]
        for file_list in files_grabbed:
            for afile in file_list:
                shutil.move(afile, 'alignment')

        os.makedirs('zips')
        for afile in self.fastq_list:
            shutil.move(afile, 'zips')

        runtime = (datetime.now() - self.startTime)
        print(f'\n{self.sample_name} alignment/VCF runtime: {runtime}\n')

        #stats to excel
        df = pd.DataFrame(index=[fq.sample_name], columns=['Reference', 'Read1 FASTQ', 'Read1 File Size', 'Read1 Total Reads', 'Read1 Mean Read Length', 'Read1 Mean Read Quality', 'Read1 Reads Passing Q30', 'Read2 FASTQ', 'Read2 File Size', 'Read2 Total Reads', 'Read2 Mean Read Length', 'Read2 Mean Read Quality', 'Read2 Reads Passing Q30', 'Total Reads', 'All Mapped Reads', 'Reference with Coverage', 'Average Depth of Coverage', 'Unmapped Reads', 'Unmapped Assembled Contigs', 'Good SNP Count', 'Group Placements', 'Bovis SB Code', 'Brucella MLST'])
        if self.species:
            reference = self.species
        else:
            reference = chrom
        df.at[fq.sample_name, 'Reference'] = f'{reference}'
        df.at[fq.sample_name, 'Read1 FASTQ'] = f'{fq.read1.fastq}'
        df.at[fq.sample_name, 'Read1 File Size'] = f'{fq.read1.file_size}'
        df.at[fq.sample_name, 'Read1 Total Reads'] = f'{fq.read1.total_read_count:,}'
        df.at[fq.sample_name, 'Read1 Mean Read Length'] = f'{fq.read1.length_mean:.1f}'
        df.at[fq.sample_name, 'Read1 Mean Read Quality'] = f'{fq.read1.read_average:.1f}'
        df.at[fq.sample_name, 'Read1 Reads Passing Q30'] = f'{fq.read1.reads_gt_q30/fq.read1.sampling_size:0.1%}'
        total_reads = fq.read1.total_read_count
        try:
            df.at[fq.sample_name, 'Read2 FASTQ'] = f'{fq.read2.fastq}'
            df.at[fq.sample_name, 'Read2 File Size'] = f'{fq.read2.file_size}'
            df.at[fq.sample_name, 'Read2 Total Reads'] = f'{fq.read2.total_read_count:,}'
            df.at[fq.sample_name, 'Read2 Mean Read Length'] = f'{fq.read2.length_mean:.1f}'
            df.at[fq.sample_name, 'Read2 Mean Read Quality'] = f'{fq.read2.read_average:.1f}'
            df.at[fq.sample_name, 'Read2 Reads Passing Q30'] = f'{fq.read2.reads_gt_q30/fq.read2.sampling_size:0.1%}'
            total_reads = fq.read1.total_read_count + fq.read2.total_read_count
        except AttributeError:
            #when no read2
            df.at[fq.sample_name, 'Read2 FASTQ'] = 'No Findings'
            df.at[fq.sample_name, 'Read2 File Size'] = 'No Findings'
            df.at[fq.sample_name, 'Read2 Total Reads'] = 'No Findings'
            df.at[fq.sample_name, 'Read2 Mean Read Length'] = 'No Findings'
            df.at[fq.sample_name, 'Read2 Mean Read Quality'] = 'No Findings'
            df.at[fq.sample_name, 'Read2 Reads Passing Q30'] = 'No Findings'
        df.at[fq.sample_name, 'Total Reads'] = f'{total_reads:,}'
        df.at[fq.sample_name, 'All Mapped Reads'] = f'{self.allbam_mapped_reads:,}'
        df.at[fq.sample_name, 'Reference with Coverage'] = f'{self.genome_coverage}'
        df.at[fq.sample_name, 'Average Depth of Coverage'] = f'{self.ave_coverage:.1f}X'
        df.at[fq.sample_name, 'Unmapped Reads'] = f'{self.unmapped_reads:,}'
        try:
            df.at[fq.sample_name, 'Unmapped Assembled Contigs'] = f'{self.abyss_contig_count:,}'
        except ValueError:
            df.at[fq.sample_name, 'Unmapped Assembled Contigs'] = f'Count not given'
        df.at[fq.sample_name, 'Good SNP Count'] = f'{self.good_snp_count:,}'
        df.at[fq.sample_name, 'Group Placements'] = f'{group_string}'
        if self.group == 'TB':
            df.at[fq.sample_name, 'Bovis SB Code'] = f'{spoligo.sbcode}'
            df.at[fq.sample_name, 'Brucella MLST'] = f'Not Applicable'
        elif self.group == 'Brucella':
            df.at[fq.sample_name, 'Bovis SB Code'] = f'Not Applicable'
            df.at[fq.sample_name, 'Brucella MLST'] = f'{mlst.mlst_type}'
        else:
            df.at[fq.sample_name, 'TB Spoligo Octal Code'] = f'Not Applicable'
            df.at[fq.sample_name, 'Bovis SB Code'] = f'Not Applicable'
            df.at[fq.sample_name, 'Brucella MLST'] = f'Not Applicable'           
        df.index.name = 'sample'
        df.to_excel(f'{fq.sample_name}_alignment_vcf_stats.xlsx')

if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------

    Usage:
        alignment_vcf.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz -r *fasta
        alignment_vcf.py -r1 *fastq.gz -r *fasta
        alignment_vcf.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz -r *fasta -g *gbk

    Dependencies:
        FASTA: reference
        GBK: used to annotate VCF files and tables
    ---------------------------------------------------------

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-r1', '--read1', action='store', dest='read1', required=True, help='Required: single read, R1 when Illumina read')
    parser.add_argument('-r2', '--read2', action='store', dest='read2', required=False, default=None, help='Optional: R2 Illumina read')
    parser.add_argument('-r', '--reference', action='store', dest='reference', required=True, default=None, help="Optional: Provide reference option or FASTA file.  If neither are given, no -r option, then a TB/Brucella/paraTB best reference are searched")
    parser.add_argument('-g', '--gbk', action='store', dest='gbk', required=False, default=None, help='Optional: gbk to annotate VCF file')
    parser.add_argument('-skip_assembly', '--skip_assembly', action='store_true', dest='skip_assembly', help='skip assembly of unmapped reads')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    read1 = args.read1
    read2 = args.read2
    reference = args.reference
    gbk = args.gbk
    species = None #species output via best_reference.py, ie directory name with dependents is reported in Excel stats file
    group = None #group output via best_reference.py, ie TB or Bruc which initials spoligo or MLST generation
    skip_assembly = args.skip_assembly
    
    align_reads = Align_Reads(read1, read2, reference, gbk, species, group, skip_assembly)
    align_reads.align()