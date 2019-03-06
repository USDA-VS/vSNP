---
layout: "default"
title: FAQ
weight: 7
permalink: /faq/
---

<h1><p style="text-align: center">Frequently Asked Questions</p></h1>

-----
<br>

## Samtools unable to access libraries

`samtools: error while loading shared libraries: libbz2.so.1.0: cannot open shared object file: No such file or directory`

The Samtools Anaconda install is unable to access libraries.  Fix: Install from different conda-forge channel

`~$ conda install -c conda-forge -c bioconda samtools bzip2`


The above should now print samtools list of commands.

## Picard fails

Error:  Picard fails to find Java, causing script to fail.

The Anaconda version of Picard can be tested using:

`~$ picard`

This should output the Picard command line usage.  If an error occurs check Java version.  

`~$ java -version`

Picard requires Java 1.8 (aka version 8).

## How do I start over with the Anaconda installation?

Simply remove your anaconda folder, `rm -rf ~/anaconda`, close and reopen your terminal, and restart with the setup instructions. 

## Panda's error occurs at around line 3400.

The error may represent itself as, `KeyError: 'reference_pos'` or `*** ValueError: all the input arrays must have same number of dimensions`

This may be caused by either an old version of pandas or the latest version of the sript is not being used.

History:  Pandas 0.18.1 contained a bug which required a work around.  This bug has been fixed and the script has been updated properly append tables.

Fix:  Update to the latest verion of vSNP and pandas.  To update vSNP preform a git pull or re-clone the repo.  To update pandas: `~$ conda update pandas`.  Version must be > 0.18.1.  Tested with 0.20.3.

-----
