#!/bin/sh

# Use to make filter files from the text pasted from the Excel worksheet.
# working directory does not need to be set.
#################################################################################
#   Set variables:

tbNumberV='s/_.*//' #Remove all charaters at and beyond "_"
tbNumberW='s/\..*//' #Remove all charaters at and beyond "."

# Path to txt file containing paste from Excel worksheet.
filterFile="/home/shared/mycobacterium/bovis/scriptDependents/filterFile.txt"

# Number of columns in Excel worksheet
columns=`head $filterFile | awk 'BEGIN{ FS="\t"; OFS="\t" }  END {print NF}'`

# Location filter files are output to.
output="/home/shared/mycobacterium/bovis/scriptDependents/bovisGroups/"

# Number of Computer cores
NR_CPUS=24
#################################################################################


let columns=columns+1
rm ${output}*
echo "Number of columns: $columns"

count=1
while [ $count -lt ${columns} ]; do
    echo ${count}
    filename=`awk -v x=$count 'BEGIN{FS=OFS="\t"}{print $x}' $filterFile | head -n1`
    echo "Filename: $filename"
    awk -v x=$count 'BEGIN{FS=OFS="\t"} FNR>1 {print $x}' $filterFile | grep -v "^$" > ${output}${filename}.list
    let count=count+1
done

for i in ${output}*.list; do
    (base=`basename "$i"`
    readyfile=`echo $base | sed $tbNumberW`

    touch ${output}${readyfile}.txt

    mylist=`cat $i`

        for l in $mylist; do
            pos1=`echo $l | sed 's/-/ /g' | awk '{print $1}'`
            pos2=`echo $l | sed 's/-/ /g' | awk '{print $2}'`
            #echo $pos2
                if [[ -z "$pos2" ]]
                then
                let pos2=pos1+1
                    while [ $pos1 -lt $pos2 ]; do
                        #echo $pos1
                        echo $pos1 >> ${output}${readyfile}.txt
                        let pos1=pos1+1
                    done
                else
                let pos2=pos2+1
                    while [ $pos1 -lt $pos2 ]; do
                        #echo $pos1
                        echo $pos1 >> ${output}${readyfile}.txt
                        let pos1=pos1+1
                    done
                fi
     done) &
       let count+=1
       [[ $((count%NR_CPUS)) -eq 0 ]] && wait
 done
wait

rm ${output}*.list

#
#  Created by Stuber, Tod P - APHIS on 12/16/2013.
#

# echo {1..10} | sed 's/{//g' | sed 's/}//g' | tr ' ' '\n '
