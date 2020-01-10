# Overview

vSNP is a 2-step process.
1. `vSNP_step1.py` takes single or paired FASTQ files and either a FASTA or directory name (see path_adder.py below) containing a FASTA.  For Mycobacterium tuberculosis complex and Brucella species if no FASTA or directory name is provide a "best reference" is used.
2. `vSNP_step2.py` takes a working directory containing VCF files, output from step 1 (specifically *_zc.vcf files), and builds tables and trees.

Other than a FASTA file no other dependency file is required to run vSNP, but more are recommended.  The main, additional dependency file recommended is an Excel file containing defining SNPs and filters.  A template is provided in the script installation `dependencies/template_defining_filter.xlsx`.

`path_adder.py` is used to direct vSNP to a directory containing dependencies for a specific reference.

![](./dependencies/directory_screen_shot.png)

For example, running `path_adder.py -w /full/path/to/vSNP_dependencies` (see below) all subdirectories are accessible using the `-r` option.  `vSNP_step1.py -r1 *_R1*gz -r2 *R2*gz -r Mycobacterium_AF2122` will use NC_002945v4.fasta to align reads.  It is recommended the path to dependencies be placed on shared storage available to both compute resources and subject matter expert.

For more run detail use `-h` option.

# Quick Setup

It is expected user is familiar with command-line basics, knows how to add directories to their PATH, and has installed conda.

Download vSNP from https://github.com/USDA-VS/<br>
`git clone https://github.com/USDA-VS/vSNP.git`<br>
Add to your PATH.

Change directory to vSNP and create vsnp conda environment.<br>
`conda env create`<br>
`conda activate vsnp`

Run `vSNP_step1.py -h` to see usage details.

Run `vSNP_step1.py -t` to see "Mycobacterium_AF2122" available.

Download test files:<br>
`git clone https://github.com/USDA-VS/fastq_data_set-tb_complex`<br>
Place each sample is in their own directory and run `vSNP_step1.py -r1 *_R1*gz -r2 *_R2*gz` on each directory.

cp *zc.vcf output from step 1 into their own working directory.<br>
Run `vSNP_step2.py` on this working directory.

Additional reference options can be added from GitHub repository.<br>
`git clone https://github.com/USDA-VS/vSNP_dependencies.git`<br>
`path_adder.py -w </full/path/to/vSNP_dependencies>`

# Archived version:
https://github.com/USDA-VS/vSNP_archive/tree/master/vSNP_version1