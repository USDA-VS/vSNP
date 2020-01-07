# Overview

vSNP is a 2-step process.
1. `vSNP_step1.py` takes single or paired FASTQ files and either a FASTA or directory name (see path_adder.py below) containing a FASTA.
2. `vSNP_step2.py` takes a working directory containing VCF files, output from step 1 (specifically *_zc.vcf files), and builds tables and trees.

Other than a FASTA file no other dependency file is required to run vSNP, but more are recommended.  The main, additional dependency file recommended is an Excel file containing defining SNPs and filters.  A template is provided in the script installation `dependencies/template_defining_filter.xlsx`.

`path_adder.py` is used to direct vSNP to a directory containing dependencies for a specific reference.

![](./dependencies/directory_screen_shot.png)

For example, running `path_adder.py -w /path/to/vsnp_dependencies` all subdirectories are accessible using the `-r` option.  `vSNP_step1.py -r1 *_R1*gz -r2 *R2*gz -r Mycobacterium_AF2122` will use NC_002945v4.fasta to align reads.  It is recommended the path to dependencies be placed on shared storage available to both compute resources and subject matter expert.

For more run detail use the `-h` option.

# Quick Setup

Installation expects user is familiar with the command-line, knows how to add directories to their PATH, and has installed conda.

Download vSNP from https://github.com/USDA-VS/<br>
`git clone https://github.com/USDA-VS/temp_vsnp.git`

Change directory to vSNP and create vsnp conda environment.<br>
`conda env create`<br>
`conda activate vsnp`

Download test files:<br>
`git clone https://github.com/USDA-VS/fastq_data_set-tb_complex`

Run `vSNP_step1.py -t` to see "Mycobacterium_AF2122" is available

Run `vSNP_step1.py` on test files.

cp *zc.vcf output from step 1 into their own working directory<br>
Run `vSNP_step2.py`