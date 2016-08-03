#!/bin/sh

onemismatch () {
patt=($1)
for ((i=0; i<${#patt[0]}; i++)); do
    patt+=( "${patt[0]:0:i}.${patt[0]:i+1}" )
done
regex=$(IFS='|'; echo "${patt[*]}")

}

#################################################################################

echo "**********************ID FILES START**********************"

# updated oligos, 2014-12-11
pab1="AATTGTCGGATAGCCTGGCGATAACGACGC"
pab3="CACACGCGGGCCGGAACTGCCGCAAATGAC"
pab5="GCTGAAGCGGCAGACCGGCAGAACGAATAT"
pmel="TGTCGCGCGTCAAGCGGCGTGAAATCTCTG"
psuis1="TGCGTTGCCGTGAAGCTTAATTCGGCTGAT"
psuis2="GGCAATCATGCGCAGGGCTTTGCATTCGTC"
psuis3="CAAGGCAGATGCACATAATCCGGCGACCCG"
pceti1="GTGAATATAGGGTGAATTGATCTTCAGCCG"
pceti2="TTACAAGCAGGCCTATGAGCGCGGCGTGAA"
pcanis4="CTGCTACATAAAGCACCCGGCGACCGAGTT"
pcanis="ATCGTTTTGCGGCATATCGCTGACCACAGC"
povis="CACTCAATCTTCTCTACGGGCGTGGTATCC"

# updated oligos, 2015-02-02

primertb157="CTCTTCGTATACCGTTCCGTCGTCACCATGGTCCT"
primertb7="TCACGCAGCCAACGATATTCGTGTACCGCGACGGT"
primertbbov="CTGGGCGACCCGGCCGACCTGCACACCGCGCATCA"
primertb5="CCGTGGTGGCGTATCGGGCCCCTGGATCGCGCCCT"
primertb2="ATGTCTGCGTAAAGAAGTTCCATGTCCGGGAAGTA"
primertb3="GAAGACCTTGATGCCGATCTGGGTGTCGATCTTGA"
primertb4="CGGTGTTGAAGGGTCCCCCGTTCCAGAAGCCGGTG"
primertb6="ACGGTGATTCGGGTGGTCGACACCGATGGTTCAGA"

# paraTB specific
ppara="CCTTTCTTGAAGGGTGTTCG|CGAACACCCTTCAAGAAAGG"

base=`basename $1`
forReads=`echo $1 | grep _R1`
echo "Forward Reads:  $forReads"

name=`echo $base | grep _R1`
n=`echo $name | sed 's/_.*//' | sed 's/\..*//'`
echo "working on: $n"

ab1=`grep -c $pab1 $forReads`
ab3=`grep -c $pab3 $forReads`
ab5=`grep -c $pab5 $forReads`
mel=`grep -c $pmel $forReads`
suis1=`grep -c $psuis1 $forReads`
suis2=`grep -c $psuis2 $forReads`
suis3=`grep -c $psuis3 $forReads`
ceti1=`grep -c $pceti1 $forReads`
ceti2=`grep -c $pceti2 $forReads`
canis4=`grep -c $pcanis4 $forReads`
canis=`grep -c $pcanis $forReads`
ovis=`grep -c $povis $forReads`
para=`egrep -ic $ppara $forReads`

onemismatch $primertb157
tb157=`egrep -c $regex $forReads`
echo "tb157 $tb157"

onemismatch $primertb7
tb7=`egrep -c $regex $forReads`
echo "tb7 $tb7"

onemismatch $primertbbov
tbbov=`egrep -c $regex $forReads`
echo "tbbov $tbbov"

onemismatch $primertb5
tb5=`egrep -c $regex $forReads`
echo "tb5 $tb5"

onemismatch $primertb2
tb2=`egrep -c $regex $forReads`
echo "tb2 $tb2"

onemismatch $primertb3
tb3=`egrep -c $regex $forReads`
echo "tb3 $tb3"

onemismatch $primertb4
tb4=`egrep -c $regex $forReads`
echo "tb4 $tb4"

onemismatch $primertb6
tb6=`egrep -c $regex $forReads`
echo "tb6 $tb6"

bruccounts=`echo "$ab1 $ab3 $ab5 $mel $suis1 $suis2 $suis3 $ceti1 $ceti2 $canis4 $canis $ovis"`
tbcounts=`echo "$tb157 $tb7 $tbbov $tb5 $tb2 $tb3 $tb4 $tb6"`
paracounts=`echo "$para"`
echo $bruccounts
echo $tbcounts
echo $paracounts

brucbinary=`echo $bruccounts | awk '{for(i=1;i<=NF;i++) if ($i >= 1) print 1; else print 0}' | tr -cd "[:print:]"`
tbbinary=`echo $tbcounts | awk '{for(i=1;i<=NF;i++) if ($i >= 1) print 1; else print 0}' | tr -cd "[:print:]"`
parabinary=`echo $paracounts | awk '{for(i=1;i<=NF;i++) if ($i >= 10) print 1; else print 0}' | tr -cd "[:print:]"`

echo $brucbinary
echo $tbbinary
echo $parabinary

#count every occurance of 1 in binary.
check=`echo $brucbinary | grep -o "1" | wc -l`
echo "Brucella check= $check"

if [[ $check > 0 ]]; then
	echo "Brucella species found"
	tagname=`grep $n /bioinfo11/TStuber/Results/_Brucella/bruc_tags.txt`
	i=$brucbinary

	if [ $i == 111111111111 ]
	then
		catch=`echo "*** Odd Isolate, Unexpected findings ***"`
	    exit 1

	elif [ $i == 011111111111 ]
	then
		catch=`echo "Brucella abortus bv 1, 2 or 4"`

    	`processZips.sh ab1 $catch | tee tee_processZips_out.txt` &
    	echo "$catch" > tee_bruc_oligo_identifier_out2.txt

	elif [ $i == 101111111111 ]
		then
		catch=`echo "Brucella abortus bv 3"`
    	`processZips.sh ab1 $catch | tee tee_processZips_out.txt` &
    	echo "$catch" > tee_bruc_oligo_identifier_out2.txt

	elif [ $i == 110111111111 ]
    	then
    		catch=`echo "Brucella abortus bv 5, 6 or 9"`
    		`processZips.sh ab1 $catch | tee tee_processZips_out.txt` &
		echo "$catch" > tee_bruc_oligo_identifier_out2.txt

	elif [ $i == 111011111111 ]
	then
    		catch=`echo "Brucella melitensis"`
    		`processZips.sh mel $catch | tee tee_processZips_out.txt` &
		echo "$catch" > tee_bruc_oligo_identifier_out2.txt

	elif [ $i == 111101111111 ]
	then
		catch=`echo "Brucella suis bv1"`
    		`processZips.sh suis1 $catch | tee tee_processZips_out.txt` &
    		echo "$catch" > tee_bruc_oligo_identifier_out2.txt

	elif [ $i == 111110111111 ]
	then
		catch=`echo "Brucella suis bv2"`
    		`processZips.sh suis2 $catch | tee tee_processZips_out.txt` &
    		echo "$catch" > tee_bruc_oligo_identifier_out2.txt

	elif [ $i == 111111011111 ]
	then
		catch=`echo "Brucella suis bv3"`
    		`processZips.sh suis3 $catch | tee tee_processZips_out.txt` &
    		echo "$catch" > tee_bruc_oligo_identifier_out2.txt

	elif [ $i == 111111101111 ] || [ $i == 111111100111 ]
    	then
    		catch=`echo "Brucella ceti 1"`
    		`processZips.sh ceti1 $catch | tee tee_processZips_out.txt` &
    		echo "$catch" > tee_bruc_oligo_identifier_out2.txt

	elif [ $i == 111111110111 ]
    	then
    		catch=`echo "Brucella ceti 2"`
    		`processZips.sh ceti2 $catch | tee tee_processZips_out.txt` &
    		echo "$catch" > tee_bruc_oligo_identifier_out2.txt

	elif [ $i == 111111111011 ]
    	then
    		catch=`echo "Brucella suis bv4"`
    		`processZips.sh suis4 $catch | tee tee_processZips_out.txt` &
    		echo "$catch" > tee_bruc_oligo_identifier_out2.txt

	elif [ $i == 111111111001 ]
    	then
    		catch=`echo "Brucella canis"`
    		`processZips.sh canis $catch | tee tee_processZips_out.txt` &
    		echo "$catch" > tee_bruc_oligo_identifier_out2.txt

	elif [ $i == 111111111110 ]
    	then
    		catch=`echo "Brucella ovis"`
    		`processZips.sh ovis $catch | tee tee_processZips_out.txt` &
    		echo "$catch" > tee_bruc_oligo_identifier_out2.txt

	else
		catch=`echo "*** Odd Isolate, Unexpected findings, See /home/shared/brucella/bruc_oligo_identifier_output.txt ***"`
    		echo "***bruc_oligo_identifier cannot find a pattern for $n, see line $LINENO of script***"
		echo "${n} Unable to find a reference, oligo_identifier.sh stats: Oligo counts: ${bruccounts} ${tbcounts}, Binary: ${brucbinary} ${tbbinary}" >> /scratch/report/dailyReport.txt
	fi
fi

#count every occurance of 1 in binary.
check=`echo $tbbinary | grep -o "1" | wc -l`
echo "TB complex check= $check"

if [[ $check > 0 ]]; then
        echo "TB complex species found"
        tagname=`grep $n /bioinfo11/TStuber/Results/mycobacterium/Untitled.txt`
	i=$tbbinary

        if [ $i == 11101111 ] || [ $i == 11101101 ]
        then
        catch="TB1"
	`processZips.sh $catch | tee tee_processZips_out.txt` &
        echo "TB1" >> tee_tb_oligo_identifier_out2.txt
	echo "$tbbinary" >> tee_tb_oligo_identifier_out2.txt
	echo "$tbcounts" >> tee_tb_oligo_identifier_out2.txt
	
        elif [ $i == 01100111 ]
        then
        catch="TB2"
        `processZips.sh mungi $catch | tee tee_processZips_out.txt` &
        echo "TB2" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbbinary" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbcounts" >> tee_tb_oligo_identifier_out2.txt 

        elif [ $i == 01101011 ] || [ $i == 11101011 ]
        then
        catch="TB3"
        `processZips.sh $catch | tee tee_processZips_out.txt` &
        echo "TB3" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbbinary" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbcounts" >> tee_tb_oligo_identifier_out2.txt

        elif [ $i == 01101111 ]
        then
        catch="TB4a"
        `processZips.sh $catch | tee tee_processZips_out.txt` &
        echo "TB4a" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbbinary" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbcounts" >> tee_tb_oligo_identifier_out2.txt


        elif [ $i == 01101101 ] || [ $i == 11101101 ] || [ $i == 01101111 ]
        then
        catch="TB4b"
        `processZips.sh $catch | tee tee_processZips_out.txt` &
        echo "TB4b" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbbinary" >> tee_tb_oligo_identifier_out2.txt
        "$tbcounts" >> tee_tb_oligo_identifier_out2.txt

        elif [ $i == 11111111 ]
        then
        catch="TB5"
        `processZips.sh $catch | tee tee_processZips_out.txt` &
        echo "TB5" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbbinary" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbcounts" >> tee_tb_oligo_identifier_out2.txt

        elif [ $i == 11001111 ]
        then
        catch="TB6"
        `processZips.sh $catch | tee tee_processZips_out.txt` &
        echo "TB6" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbbinary" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbcounts" >> tee_tb_oligo_identifier_out2.txt

        elif [ $i == 10101110 ]
        then
        catch="TB7"
        `processZips.sh $catch | tee tee_processZips_out.txt` &
        echo "TB7 are as TBBOV" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbbinary" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbcounts" >> tee_tb_oligo_identifier_out2.txt

        elif [ $i == 11001110 ] || [ $i == 11011110 ] || [ $i == 11001100 ]
        then
        catch="TBBOV"
        `processZips.sh $catch | tee tee_processZips_out.txt` &
        echo "TBBOV" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbbinary" >> tee_tb_oligo_identifier_out2.txt
        echo "$tbcounts" >> tee_tb_oligo_identifier_out2.txt

        else
        catch="No match found"
        echo "***oligo_identifier cannot place $n with TB reference, see line $LINENO of script***"
		echo "Oligo counts:  ${tbcounts}, Binary:  $tbbinary"
		echo "${n} Unable to find a reference, oligo_identifier.sh stats: Oligo counts: ${bruccounts} ${tbcounts}, Binary: ${brucbinary} ${tbbinary}" >> /scratch/report/dailyReport.txt
        exit 1
	fi
fi

#count every occurance of 1 in binary.
check=`echo $parabinary | grep -o "1" | wc -l`
echo "M. paratb check= $check"

if [[ $check > 0 ]]; then
echo "MAC species found"
tagname=`grep $n /bioinfo11/TStuber/Results/mycobacterium/Untitled.txt`
i=$parabinary

    if [ $i == 1 ]; then
        catch="para"
        echo "M. paratuberculosis found"
        `processZips.sh para $catch | tee tee_processZips_out.txt` &
        echo "M. paratuberculosis" >> tee_tb_oligo_identifier_out2.txt
        echo "$parabinary" >> tee_tb_oligo_identifier_out2.txt
        echo "$paracounts" >> tee_tb_oligo_identifier_out2.txt
    else
        echo "oligo_identifier.sh could not find a match for $n"
        echo "oligo_identifier.sh could not find a match for $n" >> /scratch/report/dailyReport.txt
        echo "${n} Unable to find a reference, oligo_identifier.sh stats: Oligo counts: ${bruccounts} ${tbcounts} ${paracounts}, Binary: ${brucbinary} ${tbbinary} ${parabinary}" >> /scratch/report/dailyReport.txt
    fi
fi

allbinary=`echo ${brucbinary} ${tbbinary} ${parabinary} | grep -c "1"`
if [[ $allbinary == 0  ]]; then 
	echo "${n} NEEDS SPECIAL ATTENTION!!!" >> /scratch/report/dailyReport.txt
	echo "PLEASE GIVE TOD SPECIAL INSTRUCTIONS FOR ${n}.  This sample was NOT identified as Brucella, TB complex or avium complex.  If you know something about this isolate please send me an email.  CC all on email.  Thanks!" | mutt -s "${n} NEEDS SPECIAL ATTENTION!!!" -- "tod.p.stuber@usda.gov suelee.robbe-austerman@aphis.usda.gov" #Doris.M.Bravo@aphis.usda.gov john.b.fevold@aphis.usda.gov Patrick.M.Camp@aphis.usda.gov David.T.Farrell@aphis.usda.gov" 
fi

echo "Sample ${n}, ${tagname}, Oligo counts: Bruc ${bruccounts} TB ${tbcounts} MAC ${paracounts}, Binary: Bruc ${brucbinary} TB ${tbbinary} MAC ${parabinary}, ID:  ${catch}"

#Push to logfile
echo "Sample ${n}, ${tagname}, Oligo counts: Bruc ${bruccounts} TB ${tbcounts} MAC ${paracounts}, Binary: Bruc ${brucbinary} TB ${tbbinary} MAC ${parabinary}, ID:  ${catch}" >> /scratch/report/oligo_identifier_log.txt

#
#  Created by Stuber, Tod P - APHIS on 04/11/2014.
#
