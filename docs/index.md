---
# You don't need to edit this file, it's empty on purpose.
# Edit theme's home layout instead if you wanna make some changes
# See: https://jekyllrb.com/docs/themes/#overriding-theme-defaults

title: Overview
layout: "default"
weight: 1
---

<h1><p style="text-align: center">Overview</p></h1>

-----
<br>
vSNP --> USDA APHIS Veterinary Services pipeline for Mycobacterium tuberculosis complex and Brucella sp.  Genotyping from high throughput sequence providing SNP tables and phylogentic trees with output to aid in SNP validation. 

vSNP can be installed on Linux and Mac platforms and is ran from the command-line.  It is written in Python 3 and relies on the Anaconda package manager.  It is used in a two step process.  Step 1 is called on a working directory containing paired FASTQ files.  BWA is used to align reads and SNPs are called using Freebayes, outputing VCF files.  In step 2, those VCF files are gathered to produce SNP alignments, tables and phylogenetic trees of grouped isolates.  vSNP is portable and can be ran with relatively little computer resources.  Because vSNP groups samples of similar isolates one is able to quickly validate SNP positions, and report sample comparisons as a high-quality validated SNP analysis.

Minimal computer requirements are 4 cores, and 8GB of memory, but more compute resources are advantageous when running multiple samples, FASTQ file sizes are excessively large or there are over 1,000 VCFs in a comparison.

# Step 1 - FASTQ to VCF
Step 1 is fairly straight forward.  Our main workflow includes Mycobacterium tuberculosis complex and Brucella sp. and therefore the script has been optimized for such.  The script begins by selecting the best reference (see `vSNP.py -h` for list of available references) and determining the spoliogtype for TB complex isolates and MLST for Brucella spieces.  The script then goes on to perform what is now a fairly standard bacterial SNP calling pipeline where BWA aligns and Frreebayes calls SNPs.  However, in addition to what is output by the Freebayes, zero coverage positions are added to the VCF files.  This is not a standard Freebayes output option and is done using Python.  Including zero coverage positions allows a more accurate SNP summary to be represented in step 2.

# Step 2 - VCF to SNP alignment
Step 2 is called on VCF files output from step 1.  References chosen in step 1 have been selected because they have been found to be relatively close to the isolate.  The closer the reference is to the isolate the less overall SNP calling error is seen.  VCF files analyzed in step 2 must all be output from the same reference.  Obviously VCF files analyzed by different references cannot be used in the same comparison.

In addition to using a closely related FASTA reference file, which minimizes SNP calling error, there are three additional external files, or dependencies, used to create high quality, informative SNP alignments.   The three dependent files are: filter file, defining SNPs, and gbk file. 
