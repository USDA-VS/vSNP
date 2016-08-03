#!/bin/sh

#################################################################################
#  Dependencies ---
#   Unix..ish setup :)
#################################################################################

echo "**********************START**********************"

# Starting working directory must be BWA-GATK folder with included /zips file containing 2 zipped fastq files.
echo "directory"
pwd
# Make fastq directory
mkdir ./../fastq
cp ./../zips/*.fastq.gz ./../fastq

# change working directory to /fastq
cd ./../fastq

echo "starting to unzip files"
# Unzip files
gunzip *.fastq.gz
echo "finished unzipping files"

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"

revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

n=`echo $revReads | sed 's/_.*//' | sed 's/\..*//'` #grab name minus the .vcf

##############

onemismatch () {
patt=($forward)
for ((i=0; i<${#patt[0]}; i++)); do
patt+=( "${patt[0]:0:i}.${patt[0]:i+1}" )
done
echo "${patt[*]}" | tr " " "\n" > searchpattern

patt=($reverse)
for ((i=0; i<${#patt[0]}; i++)); do
patt+=( "${patt[0]:0:i}.${patt[0]:i+1}" )
done
echo "${patt[*]}" | tr " " "\n" >> searchpattern
}

##############

#Spacer01
forward="TGATCCAGAGCCGGCGACCCTCTAT"
reverse="ATAGAGGGTCGCCGGCTCTGGATCA"
onemismatch
sp01=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer02
forward="CAAAAGCTGTCGCCCAAGCATGAGG"
reverse="CCTCATGCTTGGGCGACAGCTTTTG"
onemismatch
sp02=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer03
forward="CCGTGCTTCCAGTGATCGCCTTCTA"
reverse="TAGAAGGCGATCACTGGAAGCACGG"
onemismatch
sp03=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer05
forward="ACGTCATACGCCGACCAATCATCAG"
reverse="CTGATGATTGGTCGGCGTATGACGT"
onemismatch
sp04=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer05
forward="TTTTCTGACCACTTGTGCGGGATTA"
reverse="TAATCCCGCACAAGTGGTCAGAAAA"
onemismatch
sp05=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer06
forward="CGTCGTCATTTCCGGCTTCAATTTC"
reverse="GAAATTGAAGCCGGAAATGACGACG"
onemismatch
sp06=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer07
forward="GAGGAGAGCGAGTACTCGGGGCTGC"
reverse="GCAGCCCCGAGTACTCGCTCTCCTC"
onemismatch
sp07=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer08
forward="CGTGAAACCGCCCCCAGCCTCGCCG"
reverse="CGGCGAGGCTGGGGGCGGTTTCACG"
onemismatch
sp08=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer09
forward="ACTCGGAATCCCATGTGCTGACAGC"
reverse="GCTGTCAGCACATGGGATTCCGAGT"
onemismatch
sp09=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer10
forward="TCGACACCCGCTCTAGTTGACTTCC"
reverse="GGAAGTCAACTAGAGCGGGTGTCGA"
onemismatch
sp10=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer11
forward="GTGAGCAACGGCGGCGGCAACCTGG"
reverse="CCAGGTTGCCGCCGCCGTTGCTCAC"
onemismatch
sp11=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer12
forward="ATATCTGCTGCCCGCCCGGGGAGAT"
reverse="ATCTCCCCGGGCGGGCAGCAGATAT"
onemismatch
sp12=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer13
forward="GACCATCATTGCCATTCCCTCTCCC"
reverse="GGGAGAGGGAATGGCAATGATGGTC"
onemismatch
sp13=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer14
forward="GGTGTGATGCGGATGGTCGGCTCGG"
reverse="CCGAGCCGACCATCCGCATCACACC"
onemismatch
sp14=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer15
forward="CTTGAATAACGCGCAGTGAATTTCG"
reverse="CGAAATTCACTGCGCGTTATTCAAG"
onemismatch
sp15=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer16
forward="CGAGTTCCCGTCAGCGTCGTAAATC"
reverse="GATTTACGACGCTGACGGGAACTCG"
onemismatch
sp16=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer17
forward="GCGCCGGCCCGCGCGGATGACTCCG"
reverse="CGGAGTCATCCGCGCGGGCCGGCGC"
onemismatch
sp17=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer18
forward="CATGGACCCGGGCGAGCTGCAGATG"
reverse="CATCTGCAGCTCGCCCGGGTCCATG"
onemismatch
sp18=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer19
forward="TAACTGGCTTGGCGCTGATCCTGGT"
reverse="ACCAGGATCAGCGCCAAGCCAGTTA"
onemismatch
sp19=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer20
forward="TTGACCTCGCCAGGAGAGAAGATCA"
reverse="TGATCTTCTCTCCTGGCGAGGTCAA"
onemismatch
sp20=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer21
forward="TCGATGTCGATGTCCCAATCGTCGA"
reverse="TCGACGATTGGGACATCGACATCGA"
onemismatch
sp21=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer22
forward="ACCGCAGACGGCACGATTGAGACAA"
reverse="TTGTCTCAATCGTGCCGTCTGCGGT"
onemismatch
sp22=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer23
forward="AGCATCGCTGATGCGGTCCAGCTCG"
reverse="CGAGCTGGACCGCATCAGCGATGCT"
onemismatch
sp23=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer24
forward="CCGCCTGCTGGGTGAGACGTGCTCG"
reverse="CGAGCACGTCTCACCCAGCAGGCGG"
onemismatch
sp24=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer25
forward="GATCAGCGACCACCGCACCCTGTCA"
reverse="TGACAGGGTGCGGTGGTCGCTGATC"
onemismatch
sp25=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer26
forward="CTTCAGCACCACCATCATCCGGCGC"
reverse="GCGCCGGATGATGGTGGTGCTGAAG"
onemismatch
sp26=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer27
forward="GGATTCGTGATCTCTTCCCGCGGAT"
reverse="ATCCGCGGGAAGAGATCACGAATCC"
onemismatch
sp27=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer28
forward="TGCCCCGGCGTTTAGCGATCACAAC"
reverse="GTTGTGATCGCTAAACGCCGGGGCA"
onemismatch
sp28=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer29
forward="AAATACAGGCTCCACGACACGACCA"
reverse="TGGTCGTGTCGTGGAGCCTGTATTT"
onemismatch
sp29=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer30
forward="GGTTGCCCCGCGCCCTTTTCCAGCC"
reverse="GGCTGGAAAAGGGCGCGGGGCAACC"
onemismatch
sp30=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer31
forward="TCAGACAGGTTCGCGTCGATCAAGT"
reverse="ACTTGATCGACGCGAACCTGTCTGA"
onemismatch
sp31=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer32
forward="GACCAAATAGGTATCGGCGTGTTCA"
reverse="TGAACACGCCGATACCTATTTGGTC"
onemismatch
sp32=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer33
forward="GACATGACGGCGGTGCCGCACTTGA"
reverse="TCAAGTGCGGCACCGCCGTCATGTC"
onemismatch
sp33=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer34
forward="AAGTCACCTCGCCCACACCGTCGAA"
reverse="TTCGACGGTGTGGGCGAGGTGACTT"
onemismatch
sp34=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer35
forward="TCCGTACGCTCGAAACGCTTCCAAC"
reverse="GTTGGAAGCGTTTCGAGCGTACGGA"
onemismatch
sp35=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer36
forward="CGAAATCCAGCACCACATCCGCAGC"
reverse="GCTGCGGATGTGGTGCTGGATTTCG"
onemismatch
sp36=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer37
forward="CGCGAACTCGTCCACAGTCCCCCTT"
reverse="AAGGGGGACTGTGGACGAGTTCGCG"
onemismatch
sp37=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer38
forward="CGTGGATGGCGGATGCGTTGTGCGC"
reverse="GCGCACAACGCATCCGCCATCCACG"
onemismatch
sp38=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer39
forward="GACGATGGCCAGTAAATCGGCGTGG"
reverse="CCACGCCGATTTACTGGCCATCGTC"
onemismatch
sp39=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer40
forward="CGCCATCTGTGCCTCATACAGGTCC"
reverse="GGACCTGTATGAGGCACAGATGGCG"
onemismatch
sp40=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer41
forward="GGAGCTTTCCGGCTTCTATCAGGTA"
reverse="TACCTGATAGAAGCCGGAAAGCTCC"
onemismatch
sp41=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer42
forward="ATGGTGGGACATGGACGAGCGCGAC"
reverse="GTCGCGCTCGTCCATGTCCCACCAT"
onemismatch
sp42=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

#Spacer43
forward="CGCAGAATCGCACCGGGTGCGGGAG"
reverse="CTCCCGCACCCGGTGCGATTCTGCG"
onemismatch
sp43=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`

rm searchpattern

echo "spacer01	spacer02	spacer03	spacer04	spacer05	spacer06	spacer07	spacer08	spacer09	spacer10	spacer11	spacer12	spacer13	spacer14	spacer15	spacer16	spacer17	spacer18	spacer19	spacer20	spacer21	spacer22	spacer23	spacer24	spacer25	spacer26	spacer27	spacer28	spacer29	spacer30	spacer31	spacer32	spacer33	spacer34	spacer35	spacer36	spacer37	spacer38	spacer39	spacer40	spacer41	spacer42	spacer43" > $n.spacer.txt

echo "$sp01	$sp02	$sp03	$sp04	$sp05	$sp06	$sp07	$sp08	$sp09	$sp10	$sp11	$sp12	$sp13	$sp14	$sp15	$sp16	$sp17	$sp18	$sp19	$sp20	$sp21	$sp22	$sp23	$sp24	$sp25	$sp26	$sp27	$sp28	$sp29	$sp30	$sp31	$sp32	$sp33	$sp34	$sp35	$sp36	$sp37	$sp38	$sp39	$sp40	$sp41	$sp42	$sp43" >> $n.spacer.txt

cat $n.spacer.txt | awk 'NR==2 {for(i=1;i<=NF;i++) if ($i >= 5) print 1; else print 0}' | tr -cd "[:print:]" | fold -w3 > $n.myspacers

chmod 755 ./$n.myspacers

mybinaries=`cat ./$n.myspacers`

for i in $mybinaries; do 
if [ $i == 000 ]
then
	echo "0" >> $n.octalcode.txt
elif [ $i == 001 ]
then
	echo "1" >> $n.octalcode.txt
elif [ $i == 010 ]
	then
	echo "2" >> $n.octalcode.txt
elif [ $i == 011 ]
	then
	echo "3" >> $n.octalcode.txt
elif [ $i == 100 ]
	then
	echo "4" >> $n.octalcode.txt
elif [ $i == 101 ]
	then
	echo "5" >> $n.octalcode.txt
elif [ $i == 110 ]
	then
	echo "6" >> $n.octalcode.txt	
elif [ $i == 111 ]
	then	
	echo "7" >> $n.octalcode.txt
elif [ $i == 0 ]
	then	
	echo "0" >> $n.octalcode.txt
elif [ $i == 1 ]
	then
	echo "1" >> $n.octalcode.txt
else
	echo "***Error***" >> $n.octalcode.txt
fi
done

WGSpoligo=`cat $n.octalcode.txt | tr -cd "[:print:]"`

# Add infor to spoligoCheck.txt
echo "<----- $n ----->" >> /scratch/report/spoligoCheck.txt
echo "WGSpoligo:	$WGSpoligo" >> /scratch/report/spoligoCheck.txt

#Make FileMaker file
dateFile=`date "+%Y%m%d"`
printf "%s\t%s\n" "$n" "$WGSpoligo" >> "/bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/newFiles/${dateFile}_FileMakerSpoligoImport.txt"

# Add infor to spoligoCheck_all.txt
echo "<----- $n ----->" >> /scratch/report/spoligoCheck_all.txt
echo "WGSpoligo:	$WGSpoligo" >> /scratch/report/spoligoCheck_all.txt

# move back a directory to main sample folder
cd ..

#
#  Created by Stuber, Tod P - APHIS on 03/07/2013.
#
