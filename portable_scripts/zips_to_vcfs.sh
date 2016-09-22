#!/bin/sh

#  Clone to your home directory, place scripts folder in PATH
#  "~" as working directory call... git clone https://github.com/USDA-VS/public.git
#  Clone the repository "public" --> branch "portable" to your ~ directory
#  Place ~/public/SNP_analysis/bin into your PATH
#  Working directory must contain paired-end fastq reads
#  Reads must be included as _R1 and _R2
#  Usage: processZips.sh TBBOV

#################################################################################
#  Dependencies --- ALL DEPENDENCES MUST BE IN YOUR PATH
#   bwa, http://bio-bwa.sourceforge.net/bwa.shtml
#   samtools, http://samtools.sourceforge.net/samtools.shtml
#   picard, http://picard.sourceforge.net/command-line-overview.shtml
#   gatk, http://www.broadinstitute.org/gatk/
#   igvtools, http://www.broadinstitute.org/software/igv/igvtools_commandline
#   bamtools, https://github.com/pezmaster31/bamtools/wiki/Building-and-installing
#   abyss, http://www.bcgsc.ca/platform/bioinfo/software/abyss
#   RStudio IDE and R libraries ggplot2 and gsalib, Optional
#   File containing high quality SNPs, ~/public/SNP_analysis/script_dependents/Mycobacterium_bovis/HighestQualitySNPs.vcf
#   Reference, ~/public/SNP_analysis/script_dependents/Mycobacterium_bovis/NC_002945.fasta
#################################################################################

alias pause='read -p "$LINENO Enter"'

echo "**************************************************************************"
echo "**************************** START ${PWD##*/} ****************************"
echo "**************************************************************************"

# set path to script dependents folder
SCRIPT_DEPENDENTS="/home/tstuber/workspace/snp_analysis/portable_scripts/script_dependents"

if [ -e ! $SCRIPT_DEPENDENTS/Mycobacterium_bovis/DefiningSNPsGroupDesignations.txt ]; then
    echo "Set path to script dependents folder"
    echo "See line: $LINENO"
    exit 1
fi

picard=`which picard.jar`
if [[ -z $picard ]]; then
    echo "picard.jar not in PATH"
    echo "picard version >1.14"
    echo "Add picard.jar to PATH"
    echo "See line: $LINENO"
    exit 1
fi

BWA=`which bwa`
if [[ -z $BWA ]]; then
    echo "bwa is not in PATH"
    echo "Add bwa to PATH"
    echo "See line: $LINENO"
    exit 1
fi

SAMTOOLS=`which samtools`
    if [[ -z $SAMTOOLS ]]; then
    echo "samtools is not in PATH"
    echo "Add samtools to PATH"
    echo "See line: $LINENO"
    exit 1
fi

ABYSS=`which abyss-pe`
    if [[ -z $ABYSS ]]; then
    echo "abyss-pe is not in PATH"
    echo "Add abyss-pe to PATH"
    echo "See line: $LINENO"
    exit 1
fi

BAMTOOLS=`which bamtools`
if [[ -z $BAMTOOLS ]]; then
    echo "Bamtools is not in PATH"
    echo "Add Bamtools to PATH"
    echo "See line: $LINENO"
    exit 1
fi

GATK=`which GenomeAnalysisTK.jar`
if [[ -z $GATK ]]; then
    echo "GenomeAnalysisTK.jar is not in PATH"
    echo "Add GenomeAnalysisTK.jar to PATH"
    echo "See line: $LINENO"
    exit 1
fi

IGVTOOLS=`which igvtools.jar`
if [[ -z $IGVTOOLS ]]; then
    echo "igvtools.jar is not in PATH"
    echo "Add igvtools.jar to PATH"
    echo "See line: $LINENO"
    exit 1
fi

SPOLIGOSPACERFINDER=`which spoligoSpacerFinder.sh`
if [[ -z $SPOLIGOSPACERFINDER ]]; then
    echo "spoligoSpacerFinder.sh is not in PATH"
    echo "Add spoligoSpacerFinder.sh to PATH"
    echo "See line: $LINENO"
    exit 1
fi

echo "current directory"
pwd
startingdir=`pwd`

# Move zip files to their own directory
mkdir ./zips
mv *.fastq* ./zips
mkdir ./BWAmem-GATK
cd BWAmem-GATK/


# Make alias links in BWA-GATK directory to zip files
ls ../zips/*.fastq* | while read file; do ln -s $file; done

#################################################################################

# Lineage Bov-Afri, AF2122
if [ $1 == TBBOV ]; then
    cp $SCRIPT_DEPENDENTS/Mycobacterium_bovis/NC_002945.fasta ./
    if [[ ! -e NC_002945.fasta ]]; then
        echo "At line $LINENO check your path to reference"
        exit 1
    fi

    hqs=$SCRIPT_DEPENDENTS/Mycobacterium_bovis/HighestQualitySNPs.vcf
    if [[ -z $hqs ]]; then
        echo "Check your path to VCF containing high quality SNPs at line: $LINENO"
        exit 1
    fi

    # Run spoligoSpacerFinder.sh
    echo "Starting spoligoSpacerFinder.sh"
    ${SPOLIGOSPACERFINDER} &
    echo "Moving forward from spoligoSpacerFinder.sh"

else
    echo "Incorrect argument!  Must use argument: TBBOV"
    exit 1
fi

#################################################################################

# Grab reads and reference and place them in variables
ref=`ls | grep .fasta`
echo "Reference Input:  $ref"

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"

revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

#   retrieves reference name and name from sorted BAM file name
r=`echo $ref | sed 's/\..*//'`
n=`echo $revReads | sed 's/_.*//' | sed 's/\..*//'`

echo "***Reference naming convention:  $r"
echo "***Isolate naming convention:  $n"

samtools faidx $ref
java -Xmx4g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${r}.dict

if [ -s ${ref}.fai ] && [ -s ${r}.dict ]; then
    echo "Index and dict are present, continue script"
else
    sleep 5
    echo "Either index or dict for reference is missing, try making again"
    samtools faidx $ref
    java -Xmx4g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${r}.dict
        if [ -s ${ref}.fai ] && [ -s ${r}.dict ]; then
            read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
        fi
fi

# See echo comments
echo "***bwa index $r"
bwa index $ref

# -t sets the number of threads/cores
# -r ST	 Specify the read group in a format like ‘@RG\tID:foo\tSM:bar’ Needed for GATK
#adding -B 8 will require reads to have few mismatches to align to reference.  -B 1 will allow more mismatch per read.
echo "***Making Sam file"
bwa mem -M -t 16 -R @RG"\t"ID:"$n""\t"PL:ILLUMINA"\t"PU:"$n"_RG1_UNIT1"\t"LB:"$n"_LIB1"\t"SM:"$n" $ref $forReads $revReads > $n.sam
 
# -b	 Output in the BAM format.
# -h	 Include the header in the output.
#-F INT	 Skip alignments with bits present in INT [0]
echo "***Making Bam file"
samtools view -bh -F4 -T $ref $n.sam > $n.raw.bam

####### unmapped reads #######
#Bam with mapped and unmapped reads
samtools view -bh -T $ref $n.sam > $n.all.bam
#Strip off the unmapped reads
samtools view -h -f4 $n.all.bam > $n.unmappedReads.sam
#Create fastqs of unmapped reads to assemble
java -Xmx4g -jar ${picard} SamToFastq INPUT=$n.unmappedReads.sam FASTQ=${n}-unmapped_R1.fastq SECOND_END_FASTQ=${n}-unmapped_R2.fastq
rm $n.all.bam
rm $n.unmappedReads.sam
abyss-pe name=${n}_abyss k=64 in="${n}-unmapped_R1.fastq ${n}-unmapped_R2.fastq"

mkdir ../unmappedReads
mv ${n}-unmapped_R1.fastq ../unmappedReads
mv ${n}-unmapped_R2.fastq ../unmappedReads
mv ${n}_abyss-3.fa ../unmappedReads
mv ${n}_abyss-8.fa ../unmappedReads
mv ${n}_abyss-stats ../unmappedReads
mv *coverage* ../unmappedReads
rm *abyss*
######################

echo "***Sorting Bam"
samtools sort $n.raw.bam -o $n.sorted.bam
echo "***Indexing Bam"
samtools index $n.sorted.bam
# Remove duplicate molecules

echo "***Marking Duplicates"
java -Xmx4g -jar  ${picard} MarkDuplicates INPUT=$n.sorted.bam OUTPUT=$n.dup.bam METRICS_FILE=$n.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true

echo "***Index $n.dup.bam"
samtools index $n.dup.bam

# Creates file that is used in the next step
# locally realign reads such that the number of mismatching bases is minimized across all the reads
# http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_indels_RealignerTargetCreator.html
echo "***Realigner Target Creator"
java -Xmx4g -jar $GATK -T RealignerTargetCreator -I $n.dup.bam -R $ref -o $n.forIndelRealigner.intervals

if [ ! -e $n.forIndelRealigner.intervals ]; then
	java -Xmx4g -jar $GATK -T RealignerTargetCreator --fix_misencoded_quality_scores -I $n.dup.bam -R $ref -o $n.forIndelRealigner.intervals
fi

# Uses the RealignerTargetCreator output file to improve BAM alignment
# http://www.broadinstitute.org/gatk/guide/tagged?tag=indelrealigner
echo "***Target Intervals"
java -Xmx4g -jar $GATK -T IndelRealigner -I $n.dup.bam -R $ref -targetIntervals $n.forIndelRealigner.intervals -o $n.realignedBam.bam

if [ ! -e $n.realignedBam.bam ]; then
	echo "$n RealignedBam.bam failed to make.  Possible cause: Error in quality scores.  Try --fix_misencoded_quality_scores"
	echo "$n RealignedBam.bam failed to make.  Possible cause: Error in quality scores.  Try --fix_misencoded_quality_scores" > $n.errorReport
	#cat $n.errorReport | mutt -s "$n Alignment failure" -- tod.p.stuber@usda.gov
	java -Xmx4g -jar $GATK -T IndelRealigner --fix_misencoded_quality_scores -I $n.dup.bam -R $ref -targetIntervals $n.forIndelRealigner.intervals -o $n.realignedBam.bam
fi

# Uses a .vcf file which contains SNP calls of known high value to recalibrates base quality scores
# http://www.broadinstitute.org/gatk/guide/tagged?tag=baserecalibrator
echo "***Base Recalibrator"
java -Xmx4g -jar $GATK -T BaseRecalibrator -I $n.realignedBam.bam -R $ref -knownSites ${hqs} -o $n.recal_data.grp

if [ ! -e $n.recal_data.grp ]; then
	java -Xmx4g -jar $GATK -T BaseRecalibrator --fix_misencoded_quality_scores -I $n.realignedBam.bam -R $ref -knownSites ${hqs} -o $n.recal_data.grp
fi

# Make the finished "ready" .bam file
echo "***Print Reads"
java -Xmx4g -jar $GATK -T PrintReads -R $ref -I $n.realignedBam.bam -BQSR $n.recal_data.grp -o $n.ready-mem.bam

if [ ! -e $n.ready-mem.bam ]; then
	java -Xmx4g -jar $GATK -T PrintReads --fix_misencoded_quality_scores -R $ref -I $n.realignedBam.bam -BQSR $n.recal_data.grp -o $n.ready-mem.bam
fi

# Add zero positions to vcf
java -Xmx4g -jar $GATK -R $ref -T UnifiedGenotyper -out_mode EMIT_ALL_SITES -I ${n}.ready-mem.bam -o ${n}.allsites.vcf -nt 8
awk ' $0 ~ /#/ || $8 !~ /^AN=2;/ {print $0}' ${n}.allsites.vcf > $n.ready-mem.vcf
java -Xmx4g -jar $IGVTOOLS index $n.ready-mem.vcf

# SNP calling and .vcf making
# Threads used can be changed
# http://www.broadinstitute.org/gatk/guide/tagged?tag=unifiedgenotyper
echo "***HaplotypeCaller, aka calling SNPs"
java -Xmx4g -jar $GATK -R $ref -T HaplotypeCaller -I $n.ready-mem.bam -o $n.hapreadyAll.vcf -nct 8

echo "******Awk VCF leaving just SNPs******"
awk '/#/ || $4 ~ /^[ATGC]$/ && $5 ~ /^[ATGC]$/ {print $0}' $n.hapreadyAll.vcf > $n.hapreadyOnlySNPs.vcf

java -Xmx4g -jar $IGVTOOLS index $n.hapreadyOnlySNPs.vcf

echo "***Deleting Files"
rm $n.sam
rm $n.raw.bam
rm $n.dup.bam
rm $n.dup.bam.bai
rm $n.sorted.bam
rm $n.sorted.bam.bai
rm $n.realignedBam.bam
rm $n.realignedBam.bai
rm $forReads
rm $revReads
rm igv.log
rm ${n}.allsites.vcf
rm ${n}.allsites.vcf.idx

###################################
# The next 6 steps collect metrics
###################################

#Collect Depth of coverage info
echo "***Collect Depth of Coverage"
java -jar $GATK -T DepthOfCoverage -R $ref -I $n.ready-mem.bam --omitDepthOutputAtEachBase > $n.DepthofCoverage.xls

#Quality Score Distribution
echo "***Quality Score Distribution"
java -Xmx4g -jar ${picard} QualityScoreDistribution REFERENCE_SEQUENCE=$ref INPUT=$n.ready-mem.bam CHART_OUTPUT=$n.QualityScorceDistribution.pdf OUTPUT=$n.QualityScoreDistribution ASSUME_SORTED=true

#Mean Quality by Cycle
echo "***Mean Quality by Cycle"
java -Xmx4g -jar ${picard} CollectMultipleMetrics REFERENCE_SEQUENCE=$ref INPUT=$n.ready-mem.bam OUTPUT=$n.Quality_by_cycle PROGRAM=MeanQualityByCycle ASSUME_SORTED=true

#Collect Alignment Summary Metrics
echo "***Collect Alignment Summary Metrics"
java -Xmx4g -jar ${picard} CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=$ref INPUT=$n.ready-mem.bam OUTPUT=$n.AlignmentMetrics ASSUME_SORTED=true

#Collect GC Bias Error
echo "***Collect GC Bias Error"
java -Xmx4g -jar ${picard} CollectGcBiasMetrics REFERENCE_SEQUENCE=$ref INPUT=$n.ready-mem.bam OUTPUT=$n.CollectGcBiasMetrics CHART_OUTPUT=$n.GC.PDF ASSUME_SORTED=true

#Collect Insert Size Metrics
echo "***Collect Insert Size Metrics"
java -Xmx4g -jar ${picard} CollectInsertSizeMetrics REFERENCE_SEQUENCE=$ref INPUT=$n.ready-mem.bam HISTOGRAM_FILE=$n.InsertSize.pdf OUTPUT=$n.CollectInsertSizeMetrics ASSUME_SORTED=true

cat $n.DepthofCoverage.xls >> $n.Metrics_summary.xls
cat $n.AlignmentMetrics >> $n.Metrics_summary.xls
cat $n.CollectInsertSizeMetrics >> $n.Metrics_summary.xls

echo "***Organizing files"

#Move to QualityValues subfolder
mkdir QualityValues
mv $n.GC.PDF QualityValues/
mv $n.QualityScorceDistribution.pdf QualityValues/
mv $n.InsertSize.pdf QualityValues/
mv $n.Quality_by_cycle*.pdf QualityValues/

#Remove files
rm $n.DepthofCoverage.xls
rm $n.CollectInsertSizeMetrics
rm $n.forIndelRealigner.intervals
rm $n.recal_data.grp
rm $n.FilteredReads.xls
rm $n.Quality_by_cycle.quality_distribution_metrics
rm $n.Quality_by_cycle.quality_by_cycle_metrics
rm $n.Quality_by_cycle.alignment_summary_metrics
rm $n.CollectGcBiasMetrics
rm $n.QualityScoreDistribution

###########################
echo "***Getting stats for $n"

bamtools stats -in $n.ready-mem.bam > $n.stats2.txt

echo "fastq.gz file sizes:" >> $n.stats2.txt
ls -lh ../zips/ | awk '{print $5}' | egrep -v '^$' >> $n.stats2.txt

echo "Unmapped fastq file sizes:" >> $n.stats2.txt
ls -lh ../unmappedReads/*.fastq | awk '{print $5}' | egrep -v '^$' >> $n.stats2.txt

echo "Unmapped contig count:" >> $n.stats2.txt
grep -c ">" ../unmappedReads/${n}_abyss-3.fa >> $n.stats2.txt
echo "" >> $n.stats2.txt

bamtools coverage -in $n.ready-mem.bam | awk '{sum+=$3} END { print "Average coverage: ",sum/NR"X"}' >> $n.stats2.txt
bamtools coverage -in $n.ready-mem.bam | awk '{if ($3 < 1) ++b } END {print "Reference with coverage:  "((FNR-b)/FNR)*100 "%"}' >> $n.stats2.txt

cat $n.stats2.txt | grep -v "Failed" | grep -v "Duplicates" | grep -v "Proper-pairs" >> $n.stats.txt
rm $n.stats2.txt
echo "" >> $n.stats.txt
###########################

echo 'Mean_Read_Length:' >> $n.stats.txt
awk 'BEGIN {OFS="\t"} { print $16 }' $n.AlignmentMetrics | awk 'FNR == 10 {print $0}' >> $n.stats.txt

echo "" >> $n.stats.txt

#  Add SNP call numbers to stats.txt file
echo "Number of SNPs and Map-zero in ready-mem.vcf:" >> $n.stats.txt
egrep -v "#" $n.ready-mem.vcf | grep -c ".*" >> $n.stats.txt

echo "SNPs of AC2 and QUAL > 150:" >> $n.stats.txt
egrep -v "#" $n.ready-mem.vcf | egrep "AC=2" | awk '$6 > 150' | grep -c ".*" >> $n.stats.txt

#  Show Mean Coverage at Terminal and coverageReport
echo "Mean Coverage"
awk -v number="$n" 'BEGIN {OFS="\t"} $0 ~ number { print $1,$2,$3,$7 }' $n.Metrics_summary.xls | awk 'FNR == 2 {print $0}'

mv $n.Metrics_summary.xls QualityValues/
mv $n.stats.txt QualityValues/
rm $n.Quality_by_cycle.insert_size_metrics
rm $n.AlignmentMetrics
rm -r ${startingdir}/fastq

cp $0 ./

echo "**************************** END $n ****************************"

#
#  Created by Stuber, Tod P - APHIS on 3/21/15.
#
