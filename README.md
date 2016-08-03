Bacterial SNP Pipeline
======================

INTRODUCTION
============

The scripts contained in this directory are meant to provide a low and high resolution SNP analysis of closely related isolates.  Individuals should customize scripts to cater to their lab’s throughput.  There are two main scripts: processZips.sh and vcftofasta.sh.  processZips.sh outputs BAM and VCFs from Illumina paired-end data.  

- processZips.sh is called on a working directory containing paired files with “R1” and “R2” designations.
- VCFs output from processZips.sh are collected into a single working directory and vcftofasta.sh is called outputting alignment files to be viewed in a tree viewer of choice, and SNP tables which provide a high resolution view of SNP data.

Currently the script is specifically designed to analyze Mycobacterium tuberculosis complex and Brucella specie isolates.  The scripts are not meant to be installed and go.  A basic understanding of shell scripting is required.

CONTENTS
========

processZips.sh is straight forward once the dependencies below are installed.  The script is ready once the paths to Java programs have been manually changed in the script.  Simply place your paired-end Illumina files alone in your working directory and call the script with the required argument.

spoligoSpacerfinder.sh will output the spoligo octal code if reads are from a TB complex isolate and processZips.sh bovis is called.

bruc_oligo_identifier.sh will identify Brucella reads as Brucella abortus bv 1, 2 or 4, Brucella abortus bv 3, Brucella abortus bv 5, 6 or 9, Brucella melitensis, Brucella suis bv1, Brucella suis bv2, Brucella suis bv3, Brucella suis bv4, Brucella ceti group 1, Brucella ceti group 2, Brucella canis, or Brucella ovis and then starts processZips.sh with the correct arguement/reference.

loopFiles.sh is used to help speed up the analysis when multiple isolates are needing to be placed through the pipeline.  If a list of fastq.gz file are placed in directory and named properly the isolates will be organized into their individual folders and bruc_oligo_identifier and processZips will be called.  Paired files must have R1 and R2 designations.  They must have matching names.  Names are determined by taking any character before "_" or ".".  Therefore files 01-3941_S1_L001_R1_001.fastq.gz and 01-3941_S1_L001_R2_001.fastq.gz will be placed into a directory named “01-3941”.  If the file name begins with a “B” the sample is ran as a Brucella species.  If no “B” is present as 01-3941 the isolate is ran as TB complex.

email_loopFiles.sh is used to help communicate results to interested individuals.  It calls loopFiles.sh, which then calls bruc_oligo_identifier.sh, which then calls processZips.sh.  Please remember to change the e-mail addresses.

vcftofasta.sh takes VCFs created with the same reference and creates alignment files and SNP tables.  The “ready-mem.vcf” output by processZips.sh must be used.  This script is less straight forward than the previous scripts.  

- It is dependent on SNP positions to define groups, subgroups and clades.  
- It relies heavily on regular expressions to move files around.  Therefore some thought must go in to how files are named.  Basically if your file has a unique name right of the first “_” or “.” in its name you have a good start.  Also be sure your name is not to simple that it will find more than one dependent file.  For example, if a file is named, RB51.vcf and you have more than two RB51 being used in sample names, such as: another sample named Ref-RB51 then there will be conflicts.
- Filter files can be used to filter regions of the genome that consistently produce poor SNP results.
- Coverage files output by processZip.sh placed into their own directory can be used to factor in depth-of-coverage into the SNP analysis.

All four of these items listed can be ignored when first running the script.  The script will run and complete.  Each item can then be added into the analysis based on your individual needs. 


SYSTEM REQUIREMENTS
===================

Hardware Requirements—

The first script ran on 20 isolates from a typical MiSeq run using 2 x 250 chemistry outputting 3 mb genomes with approximately 150X average depth-of-coverage on an Apple Mac Pro with 2 x 6-core Xeon processors and 48 GB of memory will complete in 1.5 hours.  Following the first script a collection of 1000 VCFs will analyze through vcftofasta.sh in 2 hours.  Although more processing and memory is favorable, scripts have been tested successfully using a Macbook laptop with an Intel Core 2 Duo processor and 4 GB of memory.

Software Requirements—
- Place all files in SNP_analysis into the shell’s path.
- Xcode developer tools (tested on OS X 10.8 and 10.9)
- BWA, http://bio-bwa.sourceforge.net/bwa.shtml ver. 0.7.5 or greater
- Samtools, http://samtools.sourceforge.net/samtools.shtml ver 0.1.19 or greater
- Picard, http://picard.sourceforge.net/command-line-overview.shtml ver 1.100 or greater
- GATK, http://www.broadinstitute.org/gatk/ ver 3.1 or greater
- IGV tools http://www.broadinstitute.org/igv/igvtools
- bamtools https://github.com/pezmaster31/bamtools
- ABySS http://www.bcgsc.ca/platform/bioinfo/software/abyss ver 1.3.4 or greater
- RStudio IDE and R libraries ggplot2 and gsalib
- RAxML http://sco.h-its.org/exelixis/web/software/raxml/index.html

If this is unfamiliar the link below is a great place to start:
http://gatkforums.broadinstitute.org/discussion/2899/howto-install-all-software-packages-required-to-follow-the-gatk-best-practices

Script File Dependents—
See script_dependents directory in SNP_analysis.

Reference in fasta format
VCF formatted file containing high quality SNPs, for optimizing BAM file.
Defining SNP positions (One can use what is available or determine their own SNP positions)
Filter files (At the beginning of vcftofasta.sh in the argument controls the filters can be turned on or off, see the B. abortus filter files for sampling the contents.  Files must be named Group, Subgroup or Clade.  Capitalization must be followed)

OUTPUT
======

- BAM file of mapped reads
- VCF files
	- ready-mem.vcf made using UnifiedGenotyper
	- hapreadyAll.vcf made using HaplotypeCaller
	- hapreadyOnlySNPs.vcf same as above but only SNPs
- Unmapped reads and assembled contigs
- Spoligotype octal code if TB isolate
- Identification if Brucella species
- Index files

INTERPRETATION OF RESULTS
=========================

- Once a SNP occurs and establishes in a population it does not revert back
- Observations of homoplasy are rare
- Group, subgroup and clade clusters only show informative SNPs for the isolates within that cluster
- SNPs observed in a single isolate are less informative than SNPs seen in multiple isolates and therefore established in a population

CITATIONS
======
Investigation of the cause of geographic disparities in IDEXX ELISA sensitivity in serum samples from Mycobacterium bovis-infected cattle.
Trost B, Stuber T, Surujballi O, Nelson J, Robbe-Austerman S, Smith NH, Desautels L, Tikoo SK, Griebel P.
Sci Rep. 2016 Mar 7;6:22763. doi: 10.1038/srep22763.
PMID: 26949166

Descriptive Epidemiology and Whole Genome Sequencing Analysis for an Outbreak of Bovine Tuberculosis in Beef Cattle and White-Tailed Deer in Northwestern Minnesota.
Glaser L, Carstensen M, Shaw S, Robbe-Austerman S, Wunschmann A, Grear D, Stuber T, Thomsen B.
PLoS One. 2016 Jan 19;11(1):e0145735. doi: 10.1371/journal.pone.0145735. eCollection 2016.
PMID: 26785113

Identification of Source of Brucella suis Infection in Human by Using Whole-Genome Sequencing, United States and Tonga.
Quance C, Robbe-Austerman S, Stuber T, Brignole T, DeBess EE, Boyd L, LeaMaster B, Tiller R, Draper J, Humphrey S, Erdman MM.
Emerg Infect Dis. 2016 Jan;22(1):79-82. doi: 10.3201/eid2201.150843.
PMID: 26689610

Anatomical distribution of Mycobacterium bovis genotypes in experimentally infected white-tailed deer.
Thacker TC, Palmer MV, Robbe-Austerman S, Stuber TP, Waters WR.
Vet Microbiol. 2015 Oct 22;180(1-2):75-81. doi: 10.1016/j.vetmic.2015.07.006. Epub 2015 Jul 26.
PMID: 26243696

CONTACT
=======
Tod Stuber,
tod.p.stuber@usda.gov,
515-343-6935

Please contact me with questions, comments, or suggestions.
