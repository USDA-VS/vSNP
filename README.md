# Overview

## vSNP -- validate SNPs

### Whole genome sequencing for disease tracing and outbreak investigations is routinely required for high consequence diseases.  vSNP is an accreditation-friendly and robust tool, designed for easy error correction and validation of SNP calls. vSNP rapidly generates annotated SNP tables and corresponding  phylogenetic trees that can be easily sized and scaled for reporting purposes.   It is able to process large scale datasets, and it can easily accommodate multiple references efficiently.

### Features:
* Allows for the creation of customized groups by specified SNP position(s).
* Functionality to filter desired SNP positions by group or for all groups.
* Two types of output:
  1. A spreadsheet, or SNP table containing SNP calls sorted in evolutionary order with annotation and genome position information.
  2. Corresponding phylogenetic trees
* Ability for run-detail to be provided and optimized by multiple users.


vSNP is a 2-step process for efficiency.

1. `vSNP_step1.py` takes single or paired FASTQ files and either a FASTA or directory name containing a FASTA (see vsnp_path_adder.py below) as a reference. For Mycobacterium tuberculosis complex and Brucella species, if no FASTA file or directory name is provided a "best reference" is automatically selected.


   NOTE: It is only necessary to process each set of raw data once for each reference.  The VCF file from Step 1 can be saved for future analyses with additional samples with Step 2.  It is recommended to place VCF files generated with a single reference into a directory to compile a dataset for future analyses. Only VCF files generated from the same reference can be compared in Step 2.  See recommended directory structure for step 1 and 2 output below.
2. `vSNP_step2.py` builds SNP tables and corresponding phylogenetic trees when ran on a directory containing a collection of zc.vcf files output from step 1. Step 2 is able to handle large datasets with thousands of VCF files, outputting detailed comparisons in minutes.

Except for a FASTA file, no other dependency file is required to run vSNP.  However, other files are recommended. The additional dependency file that is most recommended is an Excel file containing defining SNPs and filters. This will allow for the creation of custom groups.  A template is provided at `dependencies/template_defining_filter.xlsx`.

`vsnp_path_adder.py` is used to direct vSNP to a directory containing reference options.

![](./dependencies/directory_screen_shot.png)

For example, after running the following all subdirectories are accessible using the `-r` option.  

```bash
vsnp_path_adder.py -d /full/path/to/vSNP_referenece_options
```

```bash
vSNP_step1.py -r1 *_R1*gz -r2 *_R2*gz -r Mycobacterium_AF2122
```

It is recommended the path to reference options be placed on shared storage available to both compute resources and subject matter expert.

For more run detail use `-h` option.

# Quick Setup

It is expected the setup user is familiar with the command-line and can installed conda.

Follow conda installation instructions:<br>
https://bioconda.github.io/user/install.html#set-up-channels<br>
Be sure to perform channel setup.<br>


Do not work in base.  If needed make new environmnet.

```bash
conda create --name myenv
```

Installation:

```bash
conda install vsnp
```

Run `vSNP_step1.py -h` to see usage details.

[macOS](./docs/macOS_special_instructions.md) users may need to follow these special instructions for Samtools.

Download and add reference options.  It is recommended options be placed on storage accessible to both compute resources and subject matter experts analyzing the output data.

```bash
git clone https://github.com/USDA-VS/vSNP_reference_options.git
```

Use `vsnp_path_adder.py` to add options.  See `vsnp_path_adder.py -h` for help.

Download test files:

```bash
git clone https://github.com/USDA-VS/fastq_data_set-tb_complex
```

Place each sample in its own directory, and on each directory run the following command:

```bash
vSNP_step1.py -r1 *_R1*gz -r2 *_R2*gz
```

See [help](./docs/run_guidance.md) running multiple samples at once.

Copy *zc.vcf output from step 1 into a directory for step 2. Only samples compared to the same reference can be analyzed together in step 2.<br>

Run `vSNP_step2.py` on this directory.

As with reference options, it is recommended to place output from step 1 and 2 on storage accessible to subject matter experts analyzing the data.  It may be necessary to use data from all three sources- dependencies, step 1 and step 2 to fully understand the relationships of the data.

Recommended directory structure for step 1 and 2 output:

![](./dependencies/step1_screenshot.png)

![](./dependencies/step2_screenshot.png)

# Archived version:
https://github.com/USDA-VS/vSNP_archive/tree/master/vSNP_version1