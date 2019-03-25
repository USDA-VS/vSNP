---
layout: "default"
title: Setup
weight : 2
#permalink: /Setup/
---

<h1><p style="text-align: center">Setup</p></h1>

-----
<br>

vSNP INSTALLATION
=================

## Python environment setup

Linux or macOS required.  Minimum 4 cores, 8GB memory.

Special instructions provided for [Windows 10 Unbuntu app](https://usda-vs.github.io/vSNP/windows10.html).

Script written in Python 3 and must be ran using a conda build environment.  Currently tested with Python 3.6 and 3.7. 

Anaconda is a highly trusted Python package and scientific software distrubution platform.  

If Anaconda is not yet installed follow the Anaconda instructions to install on your platform.

    https://www.anaconda.com/download/     
    
Once Anaconda is installed close and reopen your terminal.

Clone script to home directory: 

    $ git clone https://github.com/USDA-VS/vSNP.git

    $ cd vSNP
    $ conda env create
    $ conda activate vsnp

 A dependency may still be written for legacy Python.  Check using...

    $ cat $(which vcffirstheader)
 
 If print statements do not contain parenthesis fix with the following command:

    $ sed -i 's/print line.strip()/print(line.strip())/' $(which vcffirstheader)

Put `vSNP` in your $PATH, or easier run lines below to put script in your anaconda PATH.

    $ ln -s {FULL PATH TO}/vSNP/vSNP.py ~/anaconda*/bin/

## Step 1 test

Test files can be downloaded at:

    ~$ git clone https://github.com/USDA-VS/fastq_data_set-tb_complex
    ~$ git clone https://github.com/USDA-VS/fastq_data_set-brucella
    
Files have been cut to 200,000 reads, which will give around 20X coverage.  This file size is convenent for downloading and testing.  They should not be added to any currated database or used in reporting.  The complete sequence files are available in SRA.

Test by making directory containing FASTQ files your working directory.

    ~$ cd ~/fastq_data_set-tb_complex

To aid in testing make a backup of files

    $ mkdir original test; cp *gz original; mv *gz test; cd test; ls

`vSNP.py` must only see `*fastq.gz` files in the working directory.  It will exit if any other file type is found.  `vSNP` will batch FASTQs based on available computer resources.  Run vSNP on the working directory containing FASTQ files.

    $ vSNP.py


## Step 2 test

In step 1 vSNP saw FASTQ files in the working directory and ran appropriately.  In step 2 it looks for VCF files.  `vSNP.py` must only see `*vcf` files in the working directory.  It will exit if any other file type is found.  

Test using VCF test files, or better yet use the VCFs you just produced from step 1 above.  From script 1 output use VCF files in the `alignment` directory ending in *_zc.vcf.  Make a working directory containing only those VCFs and call `vSNP.py`.  
    
If using VCF test files

    $ cd ~
    
    ~$ git clone https://github.com/USDA-VS/test_files.git
    
    ~$ cd test_files/bovis

    Unzip

    $ vSNP.py
    
For list of options:
    
    $ vSNP.py -h
    
<br>

---

<br>

To deactivate the conda environment, use:
    
    $ conda deactivate # will put you back into your base environment
    $ conda env list # to see all available environments
    $ conda activate vsnp # to jump back into vsnp enviroment
