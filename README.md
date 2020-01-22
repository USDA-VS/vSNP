# Overview

## vSNP -- validate SNPs

### vSNP developed to rapidly call, validate, and compare SNPs from FASTQ files in a timely manner utilizing large data sets.

### Features:
1. Output VCF file containing SNPs (no indels) and positions containing no coverage.
2. Ability to group isolates.
3. Ability to filter positions.
4. Visualize trees and corresponding SNP calls in spreadsheet.
5. Ability for run-detail to be provided and optimized by multiple users.

vSNP is a 2-step process.
1. `vSNP_step1.py` takes single or paired FASTQ files and either a FASTA or directory name containing a FASTA (see vsnp_path_adder.py below).  For Mycobacterium tuberculosis complex and Brucella species if no FASTA or directory name is provide a "best reference" is selected.  For each sample step 1 is only ran one time to generate the zc.vcf.  Each time step 1 is ran the VCF file is added to a directory containing all other VCF files ran from past runs on the same reference. (Only VCF files generated from the same reference can be compared together.)
2. `vSNP_step2.py` builds tables and trees when ran on a directory containing the collection of zc.vcf files output from step 1.  Step 2 is designed to be ran on this complete and growing collection of VCF files, able to handle thousands of VCF files, outputting detailed comparisons in minutes.

Except for a FASTA file no other dependency file is required to run vSNP, however other files are recommended.  The main, additional dependency file recommended is an Excel file containing defining SNPs and filters.  A template is provided at `dependencies/template_defining_filter.xlsx`.

`vsnp_path_adder.py` is used to direct vSNP to a directory containing reference options.

![](./dependencies/directory_screen_shot.png)

For example, after running `vsnp_path_adder.py -d /full/path/to/vSNP_referenece_options` (see setup) all subdirectories are accessible using the `-r` option.  `vSNP_step1.py -r1 *_R1*gz -r2 *R2*gz -r Mycobacterium_AF2122` will use NC_002945v4.fasta to align reads.  It is recommended the path to reference options be placed on shared storage available to both compute resources and subject matter expert.

For more run detail use `-h` option.

# Quick Setup

It is expected the setup user is familiar with the command-line and can installed conda.

Follow conda install instructions:<br>
https://bioconda.github.io/user/install.html#set-up-channels<br>
Channel setup is important.<br>

Do not work in base.  If needed make new environmnet.<br>
`conda create --name myenv`

Install:<br>
`conda install vsnp`

Run `vSNP_step1.py -h` to see usage details.

Download and add reference options.  It is recommended options be placed on storage accessible to both compute resources and subject matter expert.<br>
`git clone https://github.com/USDA-VS/vSNP_reference_options.git`

Use `vsnp_path_adder.py` to add options.  See `vsnp_path_adder.py -h` for help.

Download test files:<br>
`git clone https://github.com/USDA-VS/fastq_data_set-tb_complex`<br>
Place each sample is in their own directory and run `vSNP_step1.py -r1 *_R1*gz -r2 *_R2*gz` on each directory.

cp *zc.vcf output from step 1 into their own directory.<br>
Run `vSNP_step2.py` on this directory.

As with reference options, output from step 1 and 2 must be available on storage accessible to subject matter expert.

Recommended directory structure for step 1 and 2 output:

![](./dependencies/step1_screenshot.png)

![](./dependencies/step2_screenshot.png)

# Archived version:
https://github.com/USDA-VS/vSNP_archive/tree/master/vSNP_version1