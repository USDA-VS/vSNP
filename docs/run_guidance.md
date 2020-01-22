# vSNP_step1.py run recommendations

Contact your system administrator for optimizing compute resources.  If help is not available the below is meant to provide guidance.  Step 1 simply takes input for a single sample.  It is up to the user to run samples to maximize their compute resources.

See `vSNP_step1.py -h` for basic usage

## Multiple sample run guidance

In directory containing multiple FASTQ.gz files:<br>
```for fastq in *.fastq.gz; do name=$(echo $fastq | sed 's/[._].*//'); mkdir -p $name; mv -v $fastq $name/; done```

Check that samples grouped as expected:<br>
`ls *`

In Bash one way to get the number of computer cores available is:  `getconf _NPROCESSORS_ONLN`<br>
Dividing this number by 6 can be a good starting point for optimizing resources for vSNP.  For example if 24 cores are available, 24/6=4, set NUM_PER_CYCLE=4

```NUM_PER_CYCLE=4; starting_dir=$(pwd); for dir in ./*/; do (echo "starting: $dir"; cd ./$dir; vSNP_step1.py -r1 *_R1*.fastq.gz -r2 *_R2*.fastq.gz; cd $starting_dir) & let count+=1; [[ $((count%NUM_PER_CYCLE)) -eq 0 ]] && wait; done```

Provide `-r` option to the above if needed.  When no `-r` is used Mycobacterium TB complex and Brucella species can be called together.  Just make sure references get properly sorted for step 2.

## HPC

If HPC resources are available talk to your system administrator to utilize multiple nodes.  Some example batch scripts are provided.

`vSNP_step1_batch_script.sh` batch script to run a specified number of samples at once.<br>
`hpc_vSNP_step1_new.sh` script to sort samples into directories and run a batch script on each directory.

## Collection stats

It is convenient after running multiple samples to see all stats in a single file.  In working directory containing all subdirectories of the multiple samples.<br>
```mkdir stats; find . -name "*stat*xlsx" -exec cp -v {} stats \;; cd stats; excel_append_files.py```