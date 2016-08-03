#!/bin/sh

# Usage: Place script in path and call from directory containing VCF files.
# VCFs: Created using GATK "Best practice" workflow
# VCFs: Zero mapped positions added.  See processZips.sh (optional)

# Dependancies: Unix --> Linux/OSX, raxmlHPC-SSE3, Newick Utilities

fastatable () {

# Make concatemer with the position and REF call.
# Factor in possible multiple chromosomes
# Get rid of duplicates in concatemer and list all the positions and REF calls
for i in *.vcf; do
	awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2, $4}' $i >> concatemer
done

# Get rid of duplicates in concatemer and list all the positions and REF calls
sort -k1,1 < concatemer | uniq > total_alt
awk '{print $1}' total_alt > total_pos

# Count the number of SNPs
totalSNPs=`wc -l  total_pos`
echo "Total SNPs: $totalSNPs" 

######################## FILTER FILE CREATOR ###########################
# filter poor QUAL and Map Quality

echo "`date` --> Finding positions to filter"
awk '{print $1}' total_pos > prepositionlist                                                   
        for n  in `cat prepositionlist`; do                                                    
        (front=`echo "$n" | sed 's/\(.*\)-\([0-9]*\)/\1/'`
        back=`echo "$n" | sed 's/\(.*\)-\([0-9]*\)/\2/'`                                       
        #echo "front: $front"
        #echo "back: $back"                                                                    

        positioncount=`awk -v f=$front -v b=$back ' $1 == f && $2 == b {count++} END {print count}' ./*vcf`
        #echo "position count: $positioncount"                                                 
        if [ $positioncount -gt 2 ]; then
                #printf "%s\t%s\n" "$front" "$back"
                echo "$n" >> positionlist
	else
		echo $n >> ${d}-DONOT_filtertheseposition.txt                                                      
        fi) &
        let count+=1                                                                           
        [[ $((count%NR_CPUS)) -eq 0 ]] && wait

        done
wait    
                
echo "`date` --> Filtering..."                                               
for p in `cat positionlist`; do                                                                
        (front=`echo "$p" | sed 's/\(.*\)-\([0-9]*\)/\1/'`                                     
        back=`echo "$p" | sed 's/\(.*\)-\([0-9]*\)/\2/'`                                       
        #echo "front: $front"                                                                  
        #echo "back: $back"                                                                    

        maxqual=`awk -v f=$front -v b=$back 'BEGIN{max=0} $1 == f && $2 == b {if ($6>max) max=$6} END {print max}' ./*vcf | sed 's/\..*//'`                                                   
	
	avequal=`awk -v f=$front -v b=$back '$6 != "." && $1 == f && $2 == b {print $6}' ./*vcf | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' | sed 's/\..*//'`

        maxmap=`awk -v f=$front -v b=$back ' $1 == f && $2 == b {print $8}' ./*vcf | sed 's/.*MQ=\(.....\).*/\1/' | awk 'BEGIN{max=0}{if ($1>max) max=$1} END {print max}' | sed 's/\..*//'`  
        
        avemap=`awk -v f=$front -v b=$back '$6 != "." && $1 == f && $2 == b {print $8}' ./*vcf | sed 's/.*MQ=\(.....\).*/\1/' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' | sed 's/\..*//'`
        
        #change maxmap from 52 to 56 2015-09-18                                                
        if [ $maxqual -lt 1300  ] || [ $avequal -lt 800 ] || [ $maxmap -lt 58  ] || [ $avemap -lt 57 ]; then
                echo "maxqual $maxqual" >> filterpositiondetail                                
                echo "avequal $avequal" >> filterpositiondetail
		echo "maxmap $maxmap" >> filterpositiondetail                                  
                echo "avemap $avemap" >> filterpositiondetail
                echo "position $p" >> filterpositiondetail
                echo ""  >> filterpositiondetail                                               
                echo "$p" >> ${d}-filtertheseposition.txt
        else
                echo "$p" >> ${d}-DONOT_filtertheseposition.txt                                
                #echo "maxqual $maxqual"                                                       
                #echo "maxmap $maxmap"
                #echo "avemap $avemap"
                #echo "position $p"
                #echo ""
        fi) &   
        let count+=1
        [[ $((count%NR_CPUS)) -eq 0 ]] && wait
        done    
wait    
sleep 10
rm positionlist
rm prepositionlist
               
# Filter VCF files
cat total_pos ${d}-DONOT_filtertheseposition.txt | sort -k1,1n | uniq -d > filtered_total_pos  
fgrep -f filtered_total_pos total_alt > filtered_total_alt   

rm ${d}-DONOT_filtertheseposition.txt

########################################################################

# Count the number of SNPs
totalSNPs=`wc -l total_pos`
echo "Total SNPs: $totalSNPs" 
filteredSNPs=`wc -l filtered_total_pos`
echo "Total SNPs after filtering $filteredSNPs"


for i in *.vcf; do
    	(n=${i%.vcf}
    	awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2, $5}' $i > $n.allsnps_alt
	#get SNPs of interest
	fgrep -f filtered_total_pos $n.allsnps_alt > $n.targetsnps_alt
	#if SNP not found in sample default call to reference, normalize.
	cat $n.targetsnps_alt filtered_total_alt | awk '{ if (a[$1]++ == 0) print $0; }' |  sort -nk1,1 > $n.filteredsnps_alt
    	
	# If position has zero map quality change alt call to -
    	# get positions being used
    	awk '{print $1}' $n.filteredsnps_alt > $n.filteredsnps_pos
    	# Get zero coverage positions.
   	 awk ' $0 !~ /^#/ && $10 ~ /\.\/\./ {print $1 "-" $2}' ${i} > ${n}.zeropositions

    	# if duplicate then zero mapped position found for sample
	cat $n.filteredsnps_pos ${n}.zeropositions | sort | uniq -d | awk '{print $1, "-"}' > ${n}.zerotomerge_alt #the - makes it and alt file

	#if zero positions found merge them to the SNPs found
	if [ -s ${n}.zerotomerge_alt ]; then
    		# merge zero updates to SNP file
    		cat ${n}.zerotomerge_alt $n.filteredsnps_alt | awk '{ if (a[$1]++ == 0) print $0; }' | sort -nk1,1 > ${n}.zerofilteredsnps_alt
    		#echo "***Found zero postions: $n"
		rm $n.filteredsnps_alt
    	else
		#echo "no zero positions found for $n"
    		mv $n.filteredsnps_alt ${n}.zerofilteredsnps_alt
    	fi
	
	rm $n.allsnps_alt
	rm $n.filteredsnps_pos
	rm $n.targetsnps_alt
    	rm ${n}.zeropositions
	rm ${n}.zerotomerge_alt
)  &
    	let count+=1
    	[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
wait

sleep 5
wait

echo "`date` --> Finding parsimony informative positions"
# Capture only positions that have more than one SNP type called at a position
cat *zerofilteredsnps_alt | sort -nk1,1 | uniq | awk '{print $1}' | uniq -d > parsimony_informative 
# This removes calls that are the same for all isolates being analyzed

# If many SNPs fgrep may not do much and be slow
fgrep -f parsimony_informative filtered_total_alt | sort -k1,1n > parsimony_filtered_total_alt
awk '{print $1}' parsimony_filtered_total_alt > parsimony_filtered_total_pos

# Create table and fasta
awk '{print $1}' parsimony_filtered_total_alt | awk 'BEGIN{print "reference_pos"}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}' >> ${d}.table.txt
awk '{print $2}' parsimony_filtered_total_alt | awk 'BEGIN{print "reference_call"}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}' >> ${d}.table.txt

for i in *zerofilteredsnps_alt; do
        (m=`basename "$i"`; n=`echo $m | sed 's/\..*//'`

        fgrep -f parsimony_filtered_total_pos $i | sort -k1,1n > $n.pretod

        ##############################################################
        # Change AC1s to IUPAC

        # get positions being used
        awk '{print $1}' ${n}.pretod > ${n}.usedpostions
        # get AC1 positions and iupac calls  that were changed to iupac
        awk -v Q="$QUAL" ' $0 !~ /#/ && $6 > Q && $8 ~ /^AC=1;/ {print $1 "-" $2, $5}' ${i%zerofilteredsnps_alt}vcf > ${n}.ac
        # get just positions of those AC1 grabbed above
        awk '{print $1}' ${n}.ac > ${n}.acpositions
        # AC duplicate positions will need to be kept
        cat ${n}.usedpostions ${n}.acpositions | sort | uniq -d > ${n}.actokeep
        # get AC1 position with iupac, these are only positions already in the pretod

        if [ -s ${n}.actokeep ]; then
                fgrep -f ${n}.actokeep ${n}.ac > ${n}.actomerge
                # merge iupac updates to filledcut
                cat ${n}.actomerge $n.pretod | awk '{ if (a[$1]++ == 0) print $0; }' | sort -nk1,1 > $n.tod
                rm ${n}.pretod
                rm ${n}.actomerge
        else
		#echo "else done"
                mv $n.pretod $n.tod
        fi
        rm ${n}.usedpostions
        rm ${n}.ac
        rm ${n}.acpositions
        rm ${n}.actokeep
        ##############################################################

        awk '{print $2}' $n.tod | tr -d [:space:] | sed "s/^/>$n;/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g' > $n.fas
        # Add each isolate to the table
        awk '{print $2}' $n.tod | awk -v number="$n" 'BEGIN{print number}1' | tr '\n' '\t' | sed 's/$//' | awk '{print $0}' >> ${d}.table.txt ) &
        let count+=1
        [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done

wait
sleep 5

#Create root sequence
awk '{print $2}' parsimony_filtered_total_alt > root
cat root | tr -cd "[:print:]" | sed "s/^/>root;/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g' > root.fas
echo "" >> root.fas

totalSNPs=`grep -c ".*" total_pos`
echo "Total informative SNPs: $totalSNPs" 

#Clean-up
rm concatemer
rm *.tod
mkdir fasta
mv *.fas ./fasta
rm total_pos
rm root
rm *vcf
rm filtered_total_alt
rm filtered_total_pos
rm $d-filtertheseposition.txt
rm parsimony_filtered_total_alt
rm parsimony_filtered_total_pos
rm parsimony_informative
rm total_alt
rm *zerofilteredsnps_alt

}

########################################################
########################################################

organizetable () {

cd fasta

echo "$d ********* RAxML started"

awk '{print $0}' *.fas | sed '/root/{N;d;}' >> fastaGroup.txt
awk '{print $0}' *.fas >> RAxMLfastaGroup.txt

raxmlHPC-SSE3 -s RAxMLfastaGroup.txt -n ${d} -m GTRCAT -p 12345 &>/dev/null && nw_reroot RAxML_bestTree.${d} root | nw_display -s -w 1000 -v 20 -b 'opacity:0' -i 'font-size:8' -l 'font-family:serif;font-style:italic' -d 'stroke-width:2;stroke:blue' - > ../${d}-tree.svg; nw_reroot RAxML_bestTree.${d} root > tableinput.${d}; nw_reroot RAxML_bestTree.${d} root > rooted_RAxML_bestTree.${d}; mv rooted_RAxML_bestTree.${d} RAxML_bestTree.${d} 
wait
rm RAxML_parsimonyTree*
for i in RAxML*Tree*; do mv $i ../${i}.tre; done

tr ":" "\n" < tableinput.${d} | tr "," "\n" | sed 's/(//g' | sed 's/)//g' | grep -v "\.[0-9]*" | grep -v "root" > cleanedAlignment.txt

awk 'NR==FNR{o[FNR]=$1; next} {t[$1]=$0} END{for(x=1; x<=FNR; x++){y=o[x]; print t[y]}}' cleanedAlignment.txt ../$d.table.txt > joined.txt
grep "reference" ../$d.table.txt > references.txt
cat references.txt joined.txt >> joined2.txt
mv joined2.txt ../$d.sortedTable.txt

rm joined.txt
rm references.txt

cd ..

#get the number of columns
columnCount=`awk '$0 !~ /^$/ {print $0}' *sortedTable.txt | awk '{print NF-1; exit}'`

#number=`jot - 1 $columnCount`
number=`seq $columnCount`
#echo "Numbers in list: $number"

#get the reference row
cat *sortedTable.txt | sed -n 2p | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' > referenceRow
cat referenceRow > out2.txt
cat *sortedTable.txt | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' > table
#remove first column from *sortedTable.txt

echo "countDif" > countOutput.txt
echo "countFirst" > firstOutput.txt

#iterate numbers up to number of columns
for n in $number; do
#use number to grab character in reference row i.e. C
letter=`awk -v x=$n '{print $x}' referenceRow`
#echo "Number: $n"
#echo "Letter: $letter"
awk -v var1="$n" '{print $var1}' table > column
grep "$letter" column | wc -l >> countOutput.txt
sed '1,2d' column > column2
cat column2 | awk -v var2=${letter} ' $0 !~ var2 {print NR; exit}' >> firstOutput.txt

done

#/var2/ {print NR; exit}' column
#Clear table2.txt
echo "" > table2.txt
#Add a \n to the end of file

cat *sortedTable.txt >> table2.txt
#sed '$a\' *sortedTable.txt >> table2.txt

#Prepare count line
cat countOutput.txt | sed 's/ //g' | tr "\n" "\t" > readyline.txt
#Create readytable.txt with count line.
cat table2.txt readyline.txt > orginizedTable2.txt
grep -v "^$" orginizedTable2.txt > orginizedTable3.txt

#Prepare firstOut line
cat firstOutput.txt | sed 's/ //g' | tr "\n" "\t" > readyFirstOut.txt
cat orginizedTable3.txt readyFirstOut.txt > orginizedTable4.txt

rm referenceRow
rm out2.txt
rm table
rm countOutput.txt
rm table2.txt
rm readyline.txt
rm orginizedTable2.txt

awk '{a[NR]=$0} END {print a[NR]; for (i=1;i<NR;i++) print a[i]}' orginizedTable4.txt > orginizedTable5.txt
awk '{a[NR]=$0} END {print a[NR]; for (i=1;i<NR;i++) print a[i]}' orginizedTable5.txt > orginizedTable6.txt

#Transpose
awk '{
for (i=1; i<=NF; i++)  {
a[NR,i] = $i
}
}
NF>p { p = NF }
END {
for(j=1; j<=p; j++) {
str=a[1,j]
for(i=2; i<=NR; i++){
str=str" "a[i,j];
}
print str
}
}' orginizedTable6.txt > orginizedTable7.txt

#Orgainize file based on 1st 2 columns
sort -n -k1 orginizedTable7.txt | sort -n -k2 > orginizedTable8.txt

#Convert spaces to tabs
awk -v OFS="\t" '$1=$1' orginizedTable8.txt > orginizedTable9.txt
awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' orginizedTable9.txt | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' > orginizedTable10.txt

#Transpose back
awk '{
for (i=1; i<=NF; i++)  {
a[NR,i] = $i
}
}
NF>p { p = NF }
END {
for(j=1; j<=p; j++) {
str=a[1,j]
for(i=2; i<=NR; i++){
str=str" "a[i,j];
}
print str
}
}' orginizedTable10.txt > orginizedTable11.txt

c=`basename $PWD`

#Convert spaces to tabs
awk -v OFS="\t" '$1=$1' orginizedTable11.txt > $c.organizedTable.txt

rm orginizedTable3.txt
rm orginizedTable4.txt
rm orginizedTable5.txt
rm orginizedTable6.txt
rm orginizedTable7.txt
rm orginizedTable8.txt
rm orginizedTable9.txt
rm orginizedTable10.txt
rm orginizedTable11.txt
rm column
rm column2
rm readyFirstOut.txt
rm firstOutput.txt
echo "Adding map qualities..."

# Add map qualities to sorted table

# Get just the position.  The chromosome must be removed
awk ' NR == 1 {print $0}' $d.sortedTable.txt | tr "\t" "\n" | sed "1d" | awk '{print NR, $0}' > $d-positions

echo "map-quality map-quality" > quality.txt
echo "`date` --> Sorted table map quality gathering for $c"
	while read p; do
		(rownumber=`echo $p | awk '{print $1}'`
		front=`echo "$p" | awk '{print $2}' | sed 's/\(.*\)-\([0-9]*\)/\1/'`
		back=`echo "$p" | awk '{print $2}' | sed 's/\(.*\)-\([0-9]*\)/\2/'`
		#echo "rownumber: $rownumber"
		#echo "front: $front"
		#echo "back: $back"
		avemap=`awk -v f=$front -v b=$back '$6 != "." && $1 == f && $2 == b {print $8}' ./starting_files/*vcf | sed 's/.*MQ=\(.....\).*/\1/' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' | sed 's/\..*//'`
		echo "$rownumber $avemap" >> quality.txt) &
	let count+=1
	[[ $((count%NR_CPUS)) -eq 0 ]] && wait
	done < $d-positions
wait

sort -nk1,1 < quality.txt | awk '{print $2}' | tr "\n" "\t" > qualitytransposed.txt

cat $d.sortedTable.txt qualitytransposed.txt | grep -v '^$' > $d-mapquality-sortedtable.txt
mv $d-mapquality-sortedtable.txt $d.sortedTable.txt

# Add map qualities to organized table

awk ' NR == 1 {print $0}' $c.organizedTable.txt | tr "\t" "\n" | sed "1d" | awk '{print NR, $0}' > $d-positions

echo "map-quality map-quality" > quality.txt
echo "`date` --> Organized table map quality gathering for $c"
	while read p; do
		(rownumber=`echo $p | awk '{print $1}'`
		front=`echo "$p" | awk '{print $2}' | sed 's/\(.*\)-\([0-9]*\)/\1/'`
		back=`echo "$p" | awk '{print $2}' | sed 's/\(.*\)-\([0-9]*\)/\2/'`
		#echo "rownumber: $rownumber"
		#echo "front: $front"
		#echo "back: $back"
		avemap=`awk -v f=$front -v b=$back '$6 != "." && $1 == f && $2 == b {print $8}' ./starting_files/*vcf | sed 's/.*MQ=\(.....\).*/\1/' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' | sed 's/\..*//'`
		echo "$rownumber $avemap" >> quality.txt) &
	let count+=1
	[[ $((count%NR_CPUS)) -eq 0 ]] && wait
	done < $d-positions
wait

sort -nk1,1 < quality.txt | awk '{print $2}' | tr "\n" "\t" > qualitytransposed.txt

cat $c.organizedTable.txt qualitytransposed.txt | grep -v '^$' > $d-mapquality-orgainizedtable.txt
mv $d-mapquality-orgainizedtable.txt $c.organizedTable.txt

rm quality.txt
rm qualitytransposed.txt
rm $d-positions

}

###############################################################
###############################################################
########################     Start     ########################
###############################################################
###############################################################
echo "***Start --> `date`"

######################## Set variables ########################
# Set the number of CPUs available
# Set to the numbe processors available
# If set to say 50, speed will improve dramatically
NR_CPUS=55
# Get working directory name
d=${PWD##*/}
# QUAL value cutoff
QUAL=300
echo "QUAL value cutoff: $QUAL"
echo ""

rax=`which raxmlHPC-SSE3`
if [ -z $rax ]; then
        echo "error: raxmlHPC-SSE3 is not in path"
        echo "see: https://github.com/stamatak/standard-RAxML"
        exit 1
fi

newutil=`which nw_reroot`
if [ -z $newutil ]; then
        echo "error: nw_reroot is not in path"
        echo "see: http://cegg.unige.ch/newick_utils"
        exit 1
fi

# Don't loose your original vcfs
mkdir starting_files
cp *vcf starting_files

for i in *vcf; do
	dos2unix $i >/dev/null 2>&1
	tr '\r' '\n' < $i | sed 's/"##/##/' > $i.temp
	mv $i.temp $i
done

# Remove VCFs with mean depth < 20X.
COUNTER=0
for i in *vcf; do 
	meandepth=`awk '$6 != "." {print $8}' $i | sed 's/.*;DP=\(...\).*/\1/' | sed 's/[;FS]//' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' | sed 's/\..*//'`
	if [ $meandepth -lt 20 ]; then
		echo "$i only ${meandepth}X mean depth and has been removed from comparision"
		rm $i
		let COUNTER=COUNTER+1
	fi	
done
echo ""
echo "$COUNTER VCFs deleted due to poor coverage."
echo ""

# Remove INDELs, because that's a different script
for i in *vcf; do 
	egrep "^#" $i > ${i%vcf}header
	egrep -v "^#" $i > ${i%vcf}body
	awk 'BEGIN{OFS="\t"} $4 !~ /..+/ {print $0}' ${i%vcf}body | awk 'BEGIN{OFS="\t"} $5 !~ /..+/ {print $0}' >> ${i%vcf}header
	mv ${i%vcf}header $i
	rm ${i%vcf}body
done

##################### Change AC1s to IUPAC ####################

echo "`date` --> Changing AC=1 to IUPAC..."

for i in *vcf; do

(awk '
BEGIN { OFS = "\t"}

{ if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /AG/ )
                print $1, $2, $3, $4, "R", $6, $7, $8, $9, $10
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /CT/ )
                print $1, $2, $3, $4, "Y", $6, $7, $8, $9, $10
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /GC/ )
                print $1, $2, $3, $4, "S", $6, $7, $8, $9, $10
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /AT/ )
                print $1, $2, $3, $4, "W", $6, $7, $8, $9, $10
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /GT/ )
                print $1, $2, $3, $4, "K", $6, $7, $8, $9, $10          
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /AC/ )
                print $1, $2, $3, $4, "M", $6, $7, $8, $9, $10          
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /GA/ )
                print $1, $2, $3, $4, "R", $6, $7, $8, $9, $10          
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /TC/ )
                print $1, $2, $3, $4, "Y", $6, $7, $8, $9, $10          
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /CG/ )
                print $1, $2, $3, $4, "S", $6, $7, $8, $9, $10          
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /TA/ )
                print $1, $2, $3, $4, "W", $6, $7, $8, $9, $10          
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /TG/ )
                print $1, $2, $3, $4, "K", $6, $7, $8, $9, $10          
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /CA/ )
                print $1, $2, $3, $4, "M", $6, $7, $8, $9, $10          
else
        print $0         
}' $i > ${i%vcf}temp
mv ${i%vcf}temp $i) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait

######################## Functions ########################
fastatable

organizetable

echo "***Done --> `date`"

# 2015-10-09 stuber
