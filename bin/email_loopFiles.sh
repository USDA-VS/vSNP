#!/bin/sh

#  email_loopFiles.sh

echo "Start Time: `date`" > /scratch/report/dailyTime
starttime=`date +%s`

echo "Please wait.  Searching for TB complex, Brucella and paratuberculosis oligos and then starting appropriate processZips.sh argument"
`loopFiles.sh` &&

echo "" >> /scratch/report/dailyTime
echo "End Time: `date`" >> /scratch/report/dailyTime
endtime=`date +%s`
runtime=$((endtime-starttime))
printf 'Runtime: %dh:%dm:%ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60)) >> /scratch/report/dailyTime

echo "e-mailing files"

cat /scratch/report/dailyTime > /scratch/report/email_processZips.txt
echo "" >> /scratch/report/email_processZips.txt
echo "ADD_MARKER" >> /scratch/report/email_processZips.txt
echo "" >> /scratch/report/dailyReport.txt
cat /scratch/report/dailyReport.txt >> /scratch/report/email_processZips.txt

cat /scratch/report/spoligoCheck.txt >> /scratch/report/email_processZips.txt
cat /scratch/report/mlstCheck.txt >> /scratch/report/email_processZips.txt

echo "ADD_MARKER" >> /scratch/report/email_processZips.txt

cat /scratch/report/dailyStats.txt >> /scratch/report/email_processZips.txt
echo "" >> /scratch/report/email_processZips.txt

grep -v '*' /scratch/report/email_processZips.txt | grep -v "Stats for BAM file" | sed 's/ADD_MARKER/******************************************/g' > /scratch/report/email_processZips2.txt

if [[ $1 == me ]]; then
	email_list="tod.p.stuber@aphis.usda.gov"
	else
	email_list="tod.p.stuber@aphis.usda.gov suelee.robbe-austerman@aphis.usda.gov patrick.m.camp@aphis.usda.gov David.T.Farrell@aphis.usda.gov Christine.R.Quance@aphis.usda.gov Robin.L.Swanson@aphis.usda.gov" 
fi

cat /scratch/report/email_processZips2.txt | mutt -s "WGS results" -- $email_list

date >> /scratch/report/mlstCheck_all.txt
cat /scratch/report/mlstCheck.txt >> /scratch/report/mlstCheck_all.txt


#
#  Created by Tod Stuber on 11/09/12.
#
