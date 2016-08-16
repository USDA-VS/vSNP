#!/bin/sh

####################### USER DEFINING VARIABLES #######################
# Paths to set
picard='/usr/local/bin/picard-tools-1.141/picard.jar'
gatk='/usr/local/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar'
summary="./mlstreport"
sendreport="/scratch/report/mlstCheck.txt"
#######################################################################

tbNumberV='s/_.*//' #Remove all charaters at and beyond "_"
tbNumberW='s/\..*//' #Remove all charaters at and beyond "."
tbNumberOnly='s/.*\([0-9]\{2\}-[0-9]\{4,6\}\).*/\1/' #Only tb Number
dropEXT='s/\(.*\)\..*/\1/' #Just drop the extention from the file name.

#######################################################################

alias pause='read -p "$LINENO Enter"'

#Reference
echo ">ST1-MLST" > ST1-MLST.fasta
echo "CGTTTCCCCAAGGAAGTGGAGGTTGCAGGCGATACGATCGATGTTGGCTACGGCCCGATCAAGGTTCATGCCGTCCGCAACCCGGCCGAACTGCCGTGGAAGGAAGAAAACGTCGATATCGCCCTTGAATGCACCGGCATTTTCACCTCGCGCGACAAGGCAGCACTTCATCTTGAAGCTGGCGCCAAGCGCGTCATCGTCTCCGCTCCCGCAGACGGTGCCGATCTCACCGTCGTCTATGGTGTCAACAACGACAAGCTGACGAAGGACCATCTGGTCATCTCCAACGCTTCGTGTACCACCAACTGCCTTGCGCCGGTGGCTCAGGTTCTCAACGATACTATCGGTATCGAAAAGGGCTTCATGACCACGATCCACTCCTATACGGGCGACCAGCCGACGCTGGACACCATGCACAAGGATCTCTACCGCGCCCGCGCCGCTGCCCTTTCCATGATCCCGACCTCGACGGGTGCGGCCAAGGCCGTCGGTCTCGTTCTGCCGGAACTGAAAGGCAAGCTCGACGGCGTTGCCATTCGCGTCCCGACCCCAAATGTCTCGGTCGTTGATCTCACCTTCATCGCCAAGCGTGAAACCACCGTTGAAGAAGTCAACAATGCGATCCGCGAAGCCGCCAATGGCCGCCTCAAGGGCATTCTCGGCTATACCGATGAGAAGCTCGTCTCGCACGACTTCAACCACGATTCCCATTCCTCGGTCTTCCACACCGACCAGACCAAGGTTATGGACGGCACCATGGTGCGTATCCTGTCGTGGTACGACAATGAATGGGGCTTCTCCAGCCGCATGAGCGACACCGCCGTCGCTTTGGGCAAGCTGATCTGATAACGGCAACGCCTCTCCTTCACTGGCGAGGCGTTTTCATTTCTTGATAAGGACCGAGAGAAGAAACATGATGTTCCGCACCCTTGACGATGCCAATGTCCAATCCAAGCGCGTGCTGGTCCGTGTTGACCTCAACGTGCCGAAAATCCGCCGTTCTGCTCGCTGGTCTTAACACCCCGGGCGTCACCACCGTGATCGAGCCGGTCATGACGCGCGATCATACGGAAAAGATGCTGCAAGACTTTGGCGCAGACCTGACGGTTGAAACCGATAAGGATGGTGTGCGCCATATCCGTATTGTCGGCCAGGGCAAGCTTACCGGCCAGACCATCGACGTGCCGGGTGATCCCTCGTCAACGGCTTTTCCGCTGGTGGCCGCCCTTCTGGTCGAAGGTTCGGAGGTCACCATCCGCAATGTGCTGATGAACCCGACCCGCACCGGCCTGATCCTGACGTTGCAGGAAATGGGGGCGGATATCGAGATCATCGATCCACGCCTTGCCGGCGGCGAGGATGTCGCCGATCTGCGCGTCAAGGCCTCGAAGCTGAAAGGCGTTGTCGTTCCGCCGGAACGTGCGCCTTCGATGATCGATGAATATCCGGTTCTGGCCATTGCCGCGTCTTTTGCGGAAGGCGAAACCGTGATGGACGGTCTCGATGAACTGCGCGTCAAGGAATCGGATCGTCTGGCGGCCGTTGCGCGCGGCCTTGAAGCCAATGGTGTCGATTGTACCGAAGGCGAGATGTCGCTGACGGTTCGTGGCCGCCCCGGCGGCAAGGGGCTGGGCGGTGGCACGGTTGCAACCCACCTCGACCACCGCATCGCGATGAGTTTCCTCGTCATGGGCCTTGCATCGGAAAAGCCGGTTACGGTGGATGACAGCACCATGATCGCCACCTCTTTCCCGGAATTCATGGGCATGATGGCGGGGCTGGGGGCGAAGATTGCCGAAAGCGGTGCAGAATGAAATCGTTCGTCGTCGCCCCGTTCATTGTCGCCATTGACGGACCGGCCGCCTCGGGCAAGGGAACCCTTGCCCGGCGGATCGCGACACATTACGGGATGCCGCATCTCGATACGGGCCTGACCTATCGCGCGGTCGCCAAGAGCCGCGCTCTGTCATTCTGGCCGTGGCAGGCCCGGTGGACGGCGACGAGATCGACCTCACCAATTGCGACTGGGTCGTGCGTCCTAAAAAGATGATCGCTGATCTGGGCTTTGAAGACGTGACCGTCCTCAATGATTTCGAGGCGCAGGCCCTTGCCGTGGTTTCGCTGGAAGGCCACCATATGGAACAGATCGGCGGCAAACCGGAGGAGGCTGTTGCCACCCGCGTCGTGCTCGGCCCCGGCACGGGCCTTGGCGTGGCAGGTCTGTTTCGCACACGTCATGCATGGGTTCCGGTTCCCGGTGAAGGCGGTCATATCGATATCGGTCCACGCACCGAACGCGACTACCAGATTTTCCCGCATATCGAACGCATCGAAGGGCGTGTCACCGGCGAGCAAATTCTTAGCGGGCGGGGCCTGCGCAACCTCTATCTGGGCATCTGCGCCGCCGACAAGATCACGCCCACCCTTGAGACGCCAGTAGACATTACATCCGCCGGACTGGACGGCAGCAATCCACAAGCCGCAGAAACGCTTGACCTCTTCGCCACCTATCTGGGGCGGCTTGCGGGCGACCTTGCGCTCATTTTCATGGCGCATGGCGGCGTTTATCTTTCGGGTGGCATCCCGGTGCGCATCCTTTCCGCCCTCAAGGCCGGTTCGTTCCGCGCAACCTTCGAGGACAAGGCCCCGCACAAGGCCATCATGCGCGACATACCGGTCCGCGTTATCACATATCAACTGGCGGCCTTAACCGGGCTTTCCGCTTTCGCCCGCACCCCCTCGCGCTTTGAAGTTTCGACCGAGGGCCGCCGCTGGCGCATGCGCCGCTAGAGCATTTCCGAGCCAAAAGTGCGAAGCGGTTCCGTTTCCCAACGAGCCGACCGCGGCTGCGCTTGCCTATGGTCTCGACAAGAGCGAAGGCAAGACCATCGCTGTCTATGACCTTGGCGGCGGTACTTTCGACGTGTCGGTTCTGGAAATCGGCGACGGCGTTTTTGAAGTGAAGTCCACCAATGGCGACACGTTCCTTGGCGGTGAAGACTTCGATATTCGTCTGGTCGAATATCTGGTTGCCGAGTTCAAGAAGGAAAGTGGCATCGACCTGAAGAACGACAAGCTTGCCCTGCAGCGCCTCAAGGAAGCTGCCGAAAAGGCCAAGATCGAACTGTCGTCCTCGCAGCAGACCGAAATCAACCTGCCGTTCATCACGGCTGACCAGACTGGCCCGAAGCATCTGGCGATCAAGCTGTCGCGCGCCAAGTTTGAAAGCCTGGTCGATGATCTCGTGCAGCGCACGGTCGAGCCGTGCAAGGCGGCGCTCAAGGATGCCGGCCTCAAGGCTGGCGAAATTGACGAAGTGGTTCTGGTCGGCGGCATGACCCGCATGCCCAAGATTCAGGAAGTCGTGAAGGCCTTCTTCGGCAAGGAACCGCACAAGGGCGTGAACCCGGATGAAGTCGTGGCCATGGGCGCGGCGATCCAGGGCGGCGTTTTGCAGGGCGACGTCAAGGACGTGCTGCTGCTCGACGTGACCCCGCTTTCGCTCGGCATTGAAACGCTGGGCGGCGTGTTCACCCGCCTGATCGAACGCAACACCACTATCCCGACCAAGAAGTCGCAGACCTTCTCCACGGCTGAGGACAACCAGTCGGCCGTGACGATCCGCGTCTTCCAGGGCGAGCGTGAAATGGCAGCCGATAACAAGCTGCTTGGACAGTTCGACCTCGTTGGCATTCCGCCACGTCCCTGCCCGGAAAGCTTGCCGATTGCCAGGAGCGCGATCCGGCCAAGTCCGAAATCTTCATCGTCGAGGGCGATTCGGCAGGCGGTTCCGCCAAGAGCGGGCGCTCGCGCCAGAATCAGGCCATTCTGCCGCTGCGCGGCAAAATCCTGAACGTGGAACGCGTGCGTTTCGACCGGATGATTTCATCCGATCAGGTGGGCACCCTCATCACGGCGCTTGGCACCTCCATCGGCAAGGATGAAACGCACGGCTTCAACGCCGACAAGCTGCGTTATCACAAGATCATCATCATGACCGACGCCGACGTCGATGGCGCCCATATTCGTACGCTTCTGCTCACCTTCTTCTTCCGGCAGATGCCGGAACTGATCGAACGCGGGCATATCTATATCGCGCAGCCGCCGCTCTATAAGGTGACACGCGGCAAGTCTTCGCAATATATCAAGAACGAAGCCGCCTTTGAGGATTTCCTCATCGAAACCGGCCTTGAAGAAACGACACTGGAACTGGTGACTGGCGAAATGCGCGCCGGGCCGGATTTGCGCTCGGTGGTGGAGGATGCGCGCACGCTGCGTCAGCTTCTGCACGGCCTGCACACCCGCTATGACCGCAGCGTGGTGGAACAGGCGGCAATTGCCGGCCTGCTCAACCCCGATGCCTCAAGGGACAATGCAACGGCACAGCATTCCGCCGATACGGTTGCCAAGCGTCTCGACATGATTTCGGAAGAGACCGAGCGCGGCTGGAGCGGCCATGTGATGGAAGACGGCGGCTATCGCTTCGAGCGTATGGTGCGCGGTGTAAAGGATATCGCCATTCTCGACATGGCCCTGCTCGGCTCGGCCGATGCCCGCCAGGTCGACCGAGATCGAGATGTATTCCCGCCTGATCCATACGGTCGATCATATCGAAGGCCGCCTGCGTGACGGCATGGATGCGTTTGACGGCTTCCTCAGCCATGCATGGGCTGTGACGGTGACAGGCGCGCCGAAGCTGTGGGCAATGCGCTTTCTTGAGGAAAACGAACGCAGCCCGCGCGCATGGTATGGCGGCGCGATCGGCATGATGCATTTCAATGGCGATATGAATACAGGGCTGACGCTGCGCACCATCCGCATCAAGGATGGTGTGGCGGAAATCCGTGCAGGGGCGACGCTTCTGTTCGATTCCAACCCTGACGAGGAAGAAGCCGAGACCGAATTGAAGGCATCGGCCATGATTGCGGCTGTGCGGGACGCACAGAAGAGCAATCAGATCGCGGAAGAAAGTGTGGCGGCAAAGGTGGGTGAGGGGGTTTCGATCCTGCTGGTCGATCACGAGGATTCCTTCGTCCATACGCTTGCCAATTATTTCCGCCAGACGGGCGCCAAGGTTTCCACCGTGCGTTCACCGGTGGCAGAGGAGATATTCGACCGCGTCAATCCCGATCTGGTGGTGTTATCGCCGGGACCGGGCTCGCCGCAGGATTTCGATTGCAAGGCGACCATCGATAAGGCGCGCAAGCGCCAGCTTCCGATTTTTGGCGTCTGCCTCGGCCTTCAGGCACTGGCGGAAGCCTATGGCGGGGCGTTGCGCCAGCTTCGCGTTCCGGTGCATGGCAAGCCTTCACGCATCCGCGTATCAAAGCCGGAGCGCATTTTCTCCGGCTTGCCGGAGGAAGTGACGGTGGGGCGTTATCATTCGATCTTCGCCGATCCTGAACGCCTGCCGGATGATTTTCTCGTCACAGCCGAAACGGAAGACGGGATCATAGCCTGCGGTGGAGGTGGTGATGGTGCCGCCGGGCTCCAGCCTGCCTGCGGATGCGGGGCTTGTCGTGTTGCCCGGCACCAAATCCACGATTGCCGATCTGCTGGCGCTGCGTGAAAACGGCTGGGACCGCGAATTGGTCGCCCATGTGAAGCGGGGCGGGCATGTGCTTGGTATTTGCGGCGGGTTTCAAATGCTTGGACGGCGGATCAGTGACCCGGCGGGTATTGAAGGCAATGTGCGCGATATCGAGGGGCTGGGCCTTCTCGATATCGAGACGATGACGGAGCCGGAAAAAGTGGTTCGCAATGTTGAGGCGGTGTCGCTGCTGCATGATGAGCCGCTGGAGGGCTATGAAATCCACATCGGGCGCACCAGCGGGCCGGATATGGCGCGGCCATTTGCGCGTATCGGCGATCATGATGATGGGGCCGTCTCGCCCGATGGTCGTATCATGGGAACCTATCTCCACGGTATTTTCAGTGCGGATCGTTTCCGCCACCACTTTTTGCGCGCGCTGGGTGTGGAAGGCGGCCAGATGAATTATCGCGAGAGCGTCGAAGAGGCTCTGGGCGAACTGGCTGAAGGGCTGGAAGCCTCGCTGGATATTGATGGCCTGTTTGCGCTGGCATGATTGACGCCGCGAAGCCGAAAGCCTAGTGTCAAACCATGTGACAGGTTTTGCCGGAACGAATCCCCGGCAATACCAAAAGGGAATGCGACGGACGGACCCACGCCGGGCGTCTTTATCGCAGCCGACCCCGCGACTGTAGAGCGGAGAGGGAAGAGGCAAGCCGGGCAACCGGCAGCCACTGGAAATCAGATGCGATAATGCAACATCGCATTTTTGCCATCTTCTCGACAGATTATCTCCACACAATGGGGCATTTCGTGCCGCAATTACCCTCGATATGTCACCCCTGTCAGCGCGGCATGGGCGGTTTACTCCCGATGCTGCCCGCCCGATAAGGGACCGCGCAAAACGTAATTTGTGTAAGGAGAATGCCATGCGCACTCTTAAGTCTCTCGTAATCGTCTCGGCTGCGCTGCTGCCGTTCTCTGCGACCGCTTTTGCTGCCGACGCCATCCAGGAACAGCCTCCGGTTCCGGCTCCGGTTGAAGTAGCTCCCCAGTATAGCTGGGCTGGTGGCTATACCGGTCTTTACCTTGGCTATGGCTGGAACAAGGCCAAGACCAGCACCGTTGGCAGCATCAAGCCTGACGATTGGAAGGCTGGCGCCTTTGCTGGCTGGAACTTCCAGCAGGACCAGATCGTATACGGTGTTGAAGGTGATGCAGGTTATTCCTGGGCCAAGAAGTCCAAGGACGGCCTGGAAGTCAAGCAGGGCTTTGAAGGCTCGCTGCGTGCCCGCGTCGGCTACGACCTGAACCCGGTTATGCCGTACCTCACGGCTGGTATTGCCGGTTCGCAGATCAAGCTTAACAACGGCTTGGACGACGAAAGCAAGTTCCGCGTGGGTTGGACGGCTGGTGCCGGTCTCGAAGCCAAGCTGACGGACAACATCCTCGGCCGCGTTGAGTACCGTTACACCCAGTACGGCAACAAGAACTATGATCTGGCCGGTACGACTGTTCGCAACAAGCTGGACACGCAGGATATCCGCGTCGGCATCGGCTACAAGTTCTAATTATAGCATAATTGGACACGGAAAACCGGACAGGCAACTGTCCGGTTTTTTGTTGTCTGCAAAGGTCGAGAAAGCGCGGCAGAGCAACGGCGGCAGCCTGATTTTCAGGGGAAATGAAGTGGAGGCTTCTGTTGCCAGGTGCCTCCGAACCCCGCCTTAAGGGGCTAACCCTAAGGACTTTAGAGTGGGTTTCCCGCACCGCCATTAGGCAGCGAGAGCATAACCCTGAGCATTGTTGTCATTTGCAACTACTCTGTTGACCCGATAACGGTGGTATCATGCCGAGTAAAAGAGCGATCTTTACACCCTTGTCGATCCTGTTTCGCCCCCGCCACAACACAGCCTGATCGGCAAGCTGTGCTGTGGTGGAGGCGCCGGGTACCGCCCCCGGGTCCAATGGGTTTATTACACCGTCCGTTTATCACCATAGTCGGCTTGCGCCGACAGGACGTATATAGGCGTGGTTTTTACCGATTGGAAGGGGGCTTGTGCGTTTTCGCGCAAGACCGACAGAGGTGGTGCGGCCCTTCCGTTCATTTTCCATTGACAGCTTCCGCGTGCTGGTCAATCCTCACAATATATCGGGATCGGCCTTGAAGAGGCTTGGCGCAGCCGGGGCGGAAACCATGGCTGAAACGGGGACGATATGCCCCAATCGAAGGAGAGTGGATATATGAGTGAATATCTCGCGGATGTCCGTCGCTATGATGCTGCCGCCGATGAGGCCGTTGTCGAGAAAATCGTCAAGCATCTTGGCATTGCGCTTCGCAATCGCGATTCCTCGCTCGTTTCGGCAAGC" >> ST1-MLST.fasta

ref="ST1-MLST.fasta"

# Grab reads and reference and place them in variables
ref=`ls | grep .fasta`
echo "Reference Input:  $ref"

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"
revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

m=`basename ${forReads}`; n=`echo $m | sed 's/_.*//' | sed 's/\..*//'`
echo "***Sample naming convention:  ${n}"

m=`basename ${ref}`; nameref=`echo $m | sed 's/\..*//'`
bwa index $ref
samtools faidx $ref

java -Xmx2g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${nameref}.dict

if [ -s ${ref}.fai ] && [ -s ${nameref}.dict ]; then
    echo "Index and dict are present, continue script"
else
    sleep 5
    echo "Either index or dict for reference is missing, try making again"
    samtools faidx $ref

    java -Xmx2g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${nameref}.dict

    if [ -s ${ref}.fai ] && [ -s ${nameref}.dict ]; then
        read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
    fi
fi

echo "Read 2 file is present"
# Using -B 1 achieved a slightly higher reference coverage, M allows compatibility with Picard
bwa mem -M -t 5 -R @RG"\t"ID:"$n""\t"PL:ILLUMINA"\t"PU:"$n"_RG1_UNIT1"\t"LB:"$n"_LIB1"\t"SM:"$n" $ref $forReads $revReads > ${n}.sam

samtools view -bh -F4 -T $ref ${n}.sam > ${n}.raw.bam #Removed the ".fasta" from the $ref

echo "***Sorting Bam"
samtools sort ${n}.raw.bam -o ${n}.sorted.bam

echo "***Indexing Bam"
samtools index ${n}.sorted.bam

java -Xmx2g -jar ${gatk} -R $ref  -T UnifiedGenotyper -glm BOTH -out_mode EMIT_ALL_SITES -I ${n}.sorted.bam -o ${n}.vcf -nct 8

rm *.sam
rm *raw.bam
rm *.sorted.bam*
rm *.vcf.idx
rm *.dict
rm *fasta*

# After running Bruc_MLST.sh, gather vcfs in to single directory
# Position 1629 was too close to the end of glk sequence.  Reads would not assemble properly to call possilbe SNP, therefore 100 bases of the gene were added.  Because of this all positions beyond this point are 100 more.  Same with position 1645 and 2693.
echo "starting MLST positions"
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^231$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 0031.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^297$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 0097.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^363$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 0163.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^398$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 0198.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^429$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 0229.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^523$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 0323.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^631$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 0431.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^730$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 0530.txt
echo "Group1 collection complete"

for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^1247$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 0647.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^1296$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 0696.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^1342$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 0742.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^1381$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 0781.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^1648$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1048.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^1685$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1085.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^1741$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1141.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^1754$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1154.txt
echo "Group2 collection complete"

for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^2165$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1165.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^2224$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1224.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^2227$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1227.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^2297$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1297.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^2300$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1300.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^2344$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1344.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^2352$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1352.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^2403$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1403.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^2530$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1530.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^2557$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1557.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^2578$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1578.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^2629$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1629.txt
echo "Group3 collection complete"

for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^3045$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1645.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^3054$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1654.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^3118$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1718.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^3295$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1895.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^3328$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1928.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^3388$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 1988.txt
echo "Group4 collection complete"

for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^3966$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 2166.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^3969$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 2169.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^4167$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 2367.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^4271$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 2471.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^4296$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 2496.txt
echo "Group5 collection complete"

for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^4893$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 2693.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^4996$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 2796.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^4998$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 2798.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^5058$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 2858.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^5248$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3048.txt
echo "Group6 collection complete"

for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^5672$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3072.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^5737$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3137.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^5928$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3328.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^5963$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3363.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^5984$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3384.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^5987$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3387.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6025$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3425.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6045$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3445.txt
echo "Group7 collection complete"

for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6498$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3498.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6499$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3499.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6572$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3572.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6627$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3627.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6715$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3715.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6735$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3735.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6745$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3745.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6785$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3785.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6810$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3810.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6828$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3828.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6845$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3845.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6864$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3864.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^6875$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3875.txt
echo "Group8 collection complete"

for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^7382$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 3982.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^7432$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 4032.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^7464$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 4064.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^7594$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 4194.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^7660$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 4260.txt
for i in *.vcf; do awk 'BEGIN{OFS="\t"} $2 ~ /^7756$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > 4356.txt
echo "Group9, Last group  collection complete"

# Combine with:
#echo "" | awk 'BEGIN {OFS="\t"}{print "Isolate", "31", "97", "163", "198", "229", "323", "431", "530", "647", "696", "742", "781", "1048", "1085", "1141", "1154", "1165", "1224", "1227", "1297", "1300", "1344", "1352", "1403", "1530", "1557", "1578", "1629", "1645", "1654", "1718", "1895", "1928", "1988", "2166", "2169", "2367", "2471", "2496", "2693", "2796", "2798", "2858", "3048", "3072", "3137", "3328", "3363", "3384", "3387", "3425", "3445", "3498", "3499", "3572", "3627", "3715", "3735", "3745", "3785", "3810", "3828", "3845", "3864", "3875", "3982", "4032", "4064", "4194", "4260", "4356"}' > data.txt

#paste *.txt | awk 'BEGIN{OFS="\t"}{print $1, $2, $4, $6, $8, $10, $12, $14, $16, $18, $20, $22, $24, $26, $28, $30, $32, $34, $36, $38, $40, $42, $44, $46, $48, $50, $52, $54, $56, $58, $60, $62, $64, $66, $68, $70, $72, $74, $76, $78, $80, $82, $84, $86, $88, $90, $92, $94, $96, $98, $100, $102, $104, $106, $108, $110, $112, $114, $116, $118, $120, $122, $124, $126, $128, $130, $132, $134, $136, $138, $140, $142}' >> data.txt

snps=`paste *.txt | awk 'BEGIN{OFS="\t"}{print $2, $4, $6, $8, $10, $12, $14, $16, $18, $20, $22, $24, $26, $28, $30, $32, $34, $36, $38, $40, $42, $44, $46, $48, $50, $52, $54, $56, $58, $60, $62, $64, $66, $68, $70, $72, $74, $76, $78, $80, $82, $84, $86, $88, $90, $92, $94, $96, $98, $100, $102, $104, $106, $108, $110, $112, $114, $116, $118, $120, $122, $124, $126, $128, $130, $132, $134,  $136, $138, $140, $142}' | tr -d '\t'`

echo "This is the string of SNPs"
echo "$snps"

echo "This is the string of SNPs" > ${n}_mlst.txt
echo "$snps" >> ${n}_mlst.txt

rm [0-9]*.txt

for i in *.vcf; do
    grep -v "#" $i | awk '{ if ($2 >= 201 && $2 <= 789) print $0 }' > $i.section1.txt
    grep -v "#" $i | awk '{ if ($2 >= 1190 && $2 <= 1754) print $0 }' > $i.section2.txt
    grep -v "#" $i | awk '{ if ($2 >= 2155 && $2 <= 2629) print $0 }' > $i.section3.txt
    grep -v "#" $i | awk '{ if ($2 >= 3030 && $2 <= 3499) print $0 }' > $i.section4.txt
    grep -v "#" $i | awk '{ if ($2 >= 3900 && $2 <= 4368) print $0 }' > $i.section5.txt
    grep -v "#" $i | awk '{ if ($2 >= 4769 && $2 <= 5254) print $0 }' > $i.section6.txt
    grep -v "#" $i | awk '{ if ($2 >= 5655 && $2 <= 6076) print $0 }' > $i.section7.txt
    grep -v "#" $i | awk '{ if ($2 >= 6477 && $2 <= 6966) print $0 }' > $i.section8.txt
    grep -v "#" $i | awk '{ if ($2 >= 7367 && $2 <= 7796) print $0 }' > $i.section9.txt

    cat $i.section*.txt > $i.sectionall.txt

name=`echo ${i%.vcf}`
    
    rm $i.section1.txt
    rm $i.section2.txt
    rm $i.section3.txt
    rm $i.section4.txt
    rm $i.section5.txt
    rm $i.section6.txt
    rm $i.section7.txt
    rm $i.section8.txt
    rm $i.section9.txt

    awk '$2 !~ /^231$|^297$|^363$|^398$|^429$|^523$|^631$|^730$|^1247$|^1296$|^1342$|^1381$|^1648$|^1685$|^1741$|^1754$|^2165$|^2224$|^2227$|^2297$|^2300$|^2344$|^2352$|^2403$|^2530$|^2557$|^2578$|^2629$|^3045$|^3054$|^3118$|^3295$|^3328$|^3388$|^3966$|^3969$|^4167$|^4271$|^4296$|^4893$|^4996$|^4998$|^5058$|^5248$|^5672$|^5737$|^5928$|^5963$|^5984$|^5987$|^6025$|^6045$|^6498$|^6499$|^6572$|^6627$|^6715$|^6735$|^6745$|^6785$|^6810$|^6828$|^6845$|^6864$|^6875$|^7382$|^7432$|^7464$|^7594$|^7660$|^7756$/ {print $0}' $i.sectionall.txt > $i.target.vcf
done

    rm *.sectionall.txt

# Find possible snps with:
awk '$5 != "." {snps[$2]} END{for(i in snps) print i}' *.target.vcf | sort -n > ./possibleSNPs.txt

# Then:
for l in `cat ./possibleSNPs.txt`; do for i in *.target.vcf; do
    awk -v x=$l 'BEGIN{OFS="\t"} $2 ~ "^"x"$" {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' $i; done > ${l}.txt
done

if [ -s possibleSNPs.txt ]; then
	echo "There may be a new MLST SNP to be using."
	echo "See reference at position:"
	echo "$name may have an unused SNP:" >> $sendreport
	cat possibleSNPs.txt >> $sendreport
	cat possibleSNPs.txt
else
	rm possibleSNPs.txt
fi

rm *.vcf
###################################################

if [ $snps == CTCCCGGGGCGACCCGATCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA ]; then
	echo "$name --> MLST type 01"
    echo "$name --> MLST type 01" >> $summary
	echo "$name --> MLST type 01" >> $sendreport
elif [ $snps == CTCCCGGGGCGACCCGAGCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA ]; then
    echo "$name --> MLST type 02"
    echo "$name --> MLST type 02" >> $summary
	echo "$name --> MLST type 02" >> $sendreport
elif [ $snps == CTCCCGTGGCGACCCGAGCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA ]; then
    echo "$name --> MLST type 03"
    echo "$name --> MLST type 03" >> $summary
	echo "$name --> MLST type 03" >> $sendreport
elif [ $snps == CTCCCGGGGCGACCCGAGCGAAGCGGGAAGGCCAAGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA ]; then
    echo "$name --> MLST type 04"
    echo "$name --> MLST type 04" >> $summary
	echo "$name --> MLST type 04" >> $sendreport
 elif [ $snps == CTCCCGGGGCGACCCGATCGAAGCGGGAAGGCCACGGCGAGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA ]; then
    echo "$name --> MLST type 05"
    echo "$name --> MLST type 05" >> $summary
	echo "$name --> MLST type 05" >> $sendreport
 elif [ $snps == TTCCTGGGGCAACCCGAGCGAGGCAGGGAGGCCGCGGCTCGTGAGCGGTCGGGCATCTGTCCCGCGGGGTA ]; then
    echo "$name --> MLST type 06"
    echo "$name --> MLST type 06" >> $summary
	echo "$name --> MLST type 06" >> $sendreport
 elif [ $snps == CTTCCTGGCCGAGCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT ]; then
    echo "$name --> MLST type 07"
    echo "$name --> MLST type 07" >> $summary
	echo "$name --> MLST type 07" >> $sendreport
 elif [ $snps == CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCCGTCTCGCGGTGC ]; then
    echo "$name --> MLST type 08"
    echo "$name --> MLST type 08" >> $summary
	echo "$name --> MLST type 08" >> $sendreport
 elif [ $snps == CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCTGTCTCGTGGTGCT ]; then
    echo "$name --> MLST type 09"
    echo "$name --> MLST type 09" >> $summary
	echo "$name --> MLST type 09" >> $sendreport
 elif [ $snps == CTTCCTGGCCGACCCGAGTGAAGGGGGGGGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT ]; then
    echo "$name --> MLST type 10"
    echo "$name --> MLST type 10" >> $summary
	echo "$name --> MLST type 10" >> $sendreport
 elif [ $snps == CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCTGTCTCGCGGTGCT ]; then
    echo "$name --> MLST type 11"
    echo "$name --> MLST type 11" >> $summary
	echo "$name --> MLST type 11" >> $sendreport
elif [ $snps == CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT ]; then
    echo "$name --> MLST type 12"
    echo "$name --> MLST type 12" >> $summary
	echo "$name --> MLST type 12" >> $sendreport
 elif [ $snps == CCCCCGGGCCGACTCGAGCGAAGCGAAGAGGCCACGGCGCGTGAGTGACCAGGCACCTATCCCACGGGGTA ]; then
    echo "$name --> MLST type 13"
    echo "$name --> MLST type 13" >> $summary
	echo "$name --> MLST type 13" >> $sendreport
 elif [ $snps == CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAGTGCGTGAGTGGCCAGGCACCTGTCCCGCGGGGTA ]; then
    echo "$name --> MLST type 14"
    echo "$name --> MLST type 14" >> $summary
	echo "$name --> MLST type 14" >> $sendreport
 elif [ $snps == CCCCCGGGCCGACCCGGGCGAAGCGGGGAGGCTACGGTGCGTGAGTGGCCAGGCACCTGTCCCGCGAGGTA ]; then
    echo "$name --> MLST type 15"
    echo "$name --> MLST type 15" >> $summary
	echo "$name --> MLST type 15" >> $sendreport
 elif [ $snps == CCCCCGGCCCGACCCGGGCGAAGCGGGGAGGCTACGGTGCGTGAGTGGCCAGGCACCTGTCCCGCGAGGTA ]; then
    echo "$name --> MLST type 16"
    echo "$name --> MLST type 16" >> $summary
	echo "$name --> MLST type 16" >> $sendreport
 elif [ $snps == CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACCTGTCCCGCAGGGTA ]; then
    echo "$name --> MLST type 17"
    echo "$name --> MLST type 17" >> $summary
	echo "$name --> MLST type 17" >> $sendreport
 elif [ $snps == CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACCTGTCCCGCAGGCTA ]; then
    echo "$name --> MLST type 18"
    echo "$name --> MLST type 18" >> $summary
	echo "$name --> MLST type 18" >> $sendreport
 elif [ $snps == CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGACACGGCGCGTGAGTGGCCAGGCACCTGTCCCGCGGGGTA ]; then
    echo "$name --> MLST type 19"
    echo "$name --> MLST type 19" >> $summary
	echo "$name --> MLST type 19" >> $sendreport
 elif [ $snps == CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACATGTCCCGCAGGGTA ]; then
    echo "$name --> MLST type 20"
    echo "$name --> MLST type 20" >> $summary
	echo "$name --> MLST type 20" >> $sendreport
 elif [ $snps == CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACATGCCCCGCAGGGTA ]; then
    echo "$name --> MLST type 21"
    echo "$name --> MLST type 21" >> $summary
	echo "$name --> MLST type 21" >> $sendreport
 elif [ $snps == CCCCCGGGCCGACCCGAGCGAGGCGGGGAGGCCACGGCGCGGGAGTGGCCAGACACCTGTCCTGCGGGGTA ]; then
    echo "$name --> MLST type 22"
    echo "$name --> MLST type 22" >> $summary
	echo "$name --> MLST type 22" >> $sendreport
 elif [ $snps == CCCCCGGGCTGACCCGAGCGAAACGGGGAAGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA ]; then
    echo "$name --> MLST type 23"
    echo "$name --> MLST type 23" >> $summary
	echo "$name --> MLST type 23" >> $sendreport
 elif [ $snps == CCCCCGGGCTGACCCGAGCGGAACGGGGAAGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA ]; then
    echo "$name --> MLST type 23x"
    echo "$name --> MLST type 23x" >> $summary
	echo "$name --> MLST type 23x" >> $sendreport
 elif [ $snps == CCCCCGGGCCGACCCGAGCAAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA ]; then
    echo "$name --> MLST type 24"
    echo "$name --> MLST type 24" >> $summary
	echo "$name --> MLST type 24" >> $sendreport
 elif [ $snps == CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA ]; then
    echo "$name --> MLST type 25"
    echo "$name --> MLST type 25" >> $summary
	echo "$name --> MLST type 25" >> $sendreport
 elif [ $snps == CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAAGCACCTGTTCCGCGGGGTA ]; then
    echo "$name --> MLST type 26"
    echo "$name --> MLST type 26" >> $summary
	echo "$name --> MLST type 26" >> $sendreport
 elif [ $snps == CCCCCGGGCCGACCCGAGCGAAGCGGGGAGACCACGGCGCATAAGTGGCCAGGCACCTGTCCCGCGGGGTA ]; then
    echo "$name --> MLST type 27"
    echo "$name --> MLST type 27" >> $summary
	echo "$name --> MLST type 27" >> $sendreport
 elif [ $snps == CCCTCGGGCCGACCTGAGCGAAGCGGGGAGACCACGGCGCATAAGTGGCCAGGCTCCTGTCCCGCGGGGTA ]; then
    echo "$name --> MLST type 28"
    echo "$name --> MLST type 28" >> $summary
	echo "$name --> MLST type 28" >> $sendreport
else
    echo "$name --> No MLST type found"
    echo "$name --> No MLST type found" >> $summary
    echo "$name --> No MLST type found" >> $sendreport
fi

cat $summary >> ${n}_mlst.txt
rm $summary
#
#  Created by Stuber, Tod P - APHIS on 2014-09-02.
#

