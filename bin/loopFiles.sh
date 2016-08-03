#/bin/sh

#  bovis_loopFiles.sh
#  Working directory must have Brucella or Bovis zip files
#  Brucelle files must begin with a "B"
#  Bovis files must NOT begin with a "B"


date >> /scratch/report/coverageReport.txt
echo "" > /scratch/report/dailyReport.txt

# Reset spoligo and bruc mlst check file
echo "" > /scratch/report/spoligoCheck.txt
echo "" > /scratch/report/mlstCheck.txt
echo "WG Spoligo Check" >> /scratch/report/spoligoCheck.txt
echo "Brucella MLST Check" >> /scratch/report/mlstCheck.txt

#Reset file
dateFile=`date "+%Y%m%d"`
printf "%s\t%s\n" "TB Number" "Octal Code" > "/bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/newFiles/${dateFile}_FileMakerSpoligoImport.txt"

echo "Isolate Total_Bases AveDep %>Q15" | awk '{ printf("%-12s %-12s %-10s %-10s\n", $1, $2, $3, $4) }' >> /scratch/report/dailyReport.txt
echo "" >> /scratch/report/dailyReport.txt
echo ""  > /scratch/report/dailyStats.txt

currentdir=`pwd`

for i in *.fastq.gz; do

n=`echo $i | sed 's/_.*//' | sed 's/\..*//'`
echo "n is : $n"

mkdir -p $n
mv $i $n/
done

for f in ./*/; do
	echo "The cd is $currentdir"
	cd $currentdir
	echo "$f started"
	cd ./$f
	mkdir ./temp
	cp *R1*.fastq.gz ./temp
	`gunzip ./temp/*R1*.fastq.gz && oligo_identifier.sh ./temp/*R1*.fastq | tee tee_oligo_identifier_out1.txt` &

done

#
#  Created by Tod Stuber on 11/05/12, modified 1/22/15.
#
