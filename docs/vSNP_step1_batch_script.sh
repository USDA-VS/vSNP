#!/bin/bash

#SBATCH --job-name="vSNP_s1.py"

#Purpose: batches FASTQ.gz in sample directories and runs all at once

sleep 5

for i in *.fastq.gz; do 
    n=`echo $i | sed 's/_.*//' | sed 's/\..*//'`
    echo "n is : $n"
    mkdir -p $n; mv $i $n/
done

currentdir=`pwd`
for f in ./*/; do 
    echo "starting: $f: using $1"
    cd ./$f
    if [[ $1 ]]; then
        echo "Using reference:  $1"
        vSNP_step1.py -r1 *_R1*.fastq.gz -r2 *_R2*.fastq.gz -r $1 &
    else
        echo "Finding best reference"
        vSNP_step1.py -r1 *_R1*.fastq.gz -r2 *_R2*.fastq.gz &
    fi
    cd $currentdir
done
wait