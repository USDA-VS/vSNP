#!/usr/bin/env bash

hflag=
rflag=
cflag=
while getopts ':hr:c::' OPTION; do
    case $OPTION in
    h) hflag=1
    ;;
    r) rflag=$OPTARG
    ;;
    c) cflag=$OPTARG
    ;;
    esac
done

help () {
    printf "\nPurpose: batches fastq.gz into dir[0-9]* and calls sbatch script on directories to run all at once\n\n"
    printf "### Working directory must contain .fastq.gz.  Other files are allowed, but no other directories, so no subdirectories in your working directory.  Just files.\n\n"
    echo "-h = help"
    echo "-r = reference, Optional"
    printf "\tIf reference not provided a best reference is selected if available\n"
    echo "-c = cpus given to each sample, Optional, Default: 8"
    printf "\tTherefore with 8 default 100 samples when 48 core per nodes 6 samples (12 files) will be ran per batch.  17 jobs will be sent.\n"
    printf "\nExamples:\n"
    printf '\thcp_vSNP_ste1_new.sh -r </path_to/reference.fasta|reference_option>\n\n'
    printf '\thcp_vSNP_ste1_new.sh # use best reference if tb or brucella\n\n'
    printf 'Get reference options with: vSNP_step1.py -t\n\n'
    printf 'Clean up from starting directory:\n'
    printf '\tmv ./dir*/* .; rm -r dir*; ls\n'
    printf 'mv ./dir*/* .; rm -r dir*; ls\n\n'
    printf '\tUsing zsh: mkdir stats; cp **/*.xlsx stats; cd stats; excel_append_files.py; mv combined_excelworksheets.xlsx ../; cd ..; rm -r stats\n\n'
    exit 1
}

if [ "$hflag" ]; then
    help
    exit 1
fi

if [ "$cflag" ]; then
    cpus_alloted_per_sample=$cflag
else
    cpus_alloted_per_sample=8
fi

reference=$rflag
printf "\nSpecified reference: $reference\n"

cpu_available=$(grep -c ^processor /proc/cpuinfo)
samples_per_node=$(awk -v x=$cpu_available -v y=$cpus_alloted_per_sample 'BEGIN { rounded = sprintf("%.0f", x/y); print rounded }')
files_per_node=$(( samples_per_node*2 )) #multiply by 2 for paired reads

printf "\nNumber of cpus:  $cpu_available\n"
printf "cpus begin used per sample:  $cpus_alloted_per_sample\n"
printf "Samples per node  $samples_per_node\n"

counter=0
file_count=$(ls *gz | wc -l)
while [ $file_count -gt 0 ]; do
    counter=$[$counter+1]
    counter_formated=$(printf "%02d" ${counter})
    mkdir dir${counter_formated}
    mv $(ls *gz | head -${files_per_node}) dir${counter_formated}
    echo "Moving $(ls dir${counter_formated}/*gz | wc -l) files to dir${counter_formated}"
    file_count=$(ls *gz | wc -l)
done
printf "\nNumber of nodes being used:  $counter\n"

root_dir=`pwd`; 
for each_dir in ./*/; do 
    cd $root_dir
    printf "$each_dir started: "
    cd ./$each_dir
    sbatch vSNP_step1_batch_script.sh ${reference}
done
cd $root_dir

printf '\nmv ./dir*/* .; rm -r dir*; ls\n'
printf 'Using zsh: mkdir stats; cp **/*.xlsx stats; cd stats; excel_append_files.py\n\n'