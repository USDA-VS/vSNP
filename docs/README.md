# vSNP_step1.py run recommendations

See `vSNP_step1.py -h` for basic usage

Step 1 simply takes input for a single sample.  It is up to the user to run samples in the context of their compute environment to maximize resources.  Below are recommendations.

## Group sample into their own working directories
`for i in *.fastq.gz; do n=$(echo $i | sed 's/[._].*//'); mkdir -p $n; mv -v $i $n/; done`

If after grouping there are few enough samples to be handled at one time then the following could be used to run all samples at once.

`currentdir=$(pwd); for f in ./*/; do echo "starting: $f; cd ./$f; vSNP_step1.py -r1 *_R1*.fastq.gz -r2 *_R2*.fastq.gz -r <reference FASTA/option> & cd $currentdir; done`

Or if the reference option will be Mycobactium tuberculosis complex and/or Brucella species no `-r` option is required.  Also a mix of option types can be ran together.

`currentdir=$(pwd); for f in ./*/; do echo "starting: $f; cd ./$f; vSNP_step1.py -r1 *_R1*.fastq.gz -r2 *_R2*.fastq.gz & cd $currentdir; done`

The examples above will not work if only small amount of compute is available or many (>5) samples are ran.  In these cases bash scripts are recommended.  Example scripts are available above as examples.  Your system administrar will often want to help in situtations when maximizing compute resources, and will be able to answer system specific questions.

`vSNP_step1_batch_script.sh` 

After running multiple samples it is convenient to see all stats in a single file.  In working directory containing subdirectories of multiple samples.<br>
`mkdir stats; find . -name "*stat*xlsx" -exec cp -v {} stats \;; excel_append_files.py`