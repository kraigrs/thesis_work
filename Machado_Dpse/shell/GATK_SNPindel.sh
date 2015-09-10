#!/bin/sh

#PBS -S /bin/sh
#PBS -N GATK_SNPindel
#PBS -A wittkopp_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l nodes=1:ppn=12,pmem=4000mb,walltime=10:00:00
#PBS -m abe
#PBS -j oe
#PBS -V

####  End PBS preamble
#  Include the next three lines always

if [ "x${PBS_NODEFILE}" != "x" ] ; then
   cat $PBS_NODEFILE   # contains a list of the CPUs you were using if run with PBS
fi

cd $PBS_O_WORKDIR

#Usage: var1=<ref>.[fa/fasta] var2=<file>.bam var3=<data directory> qsub samtools_SNPindel.sh

# merge and rmdup before these 3 steps

java -Xmx4g -jar $Picard_JAVA/AddOrReplaceReadGroups.jar INPUT=$var3\/$var2.bam OUTPUT=$var3\/$var2.group.bam RGLB=$var4 RGPL=Illumina RGPU=HiSeq RGSM=$var5

java -Xmx4g -jar $Picard_JAVA/SortSam.jar INPUT=$var3\/$var2.group.bam OUTPUT=$var3\/$var2.group.picard.bam SORT_ORDER=coordinate TMP_DIR=$scratch

java -Xmx4g -jar $Picard_JAVA/BuildBamIndex.jar INPUT=$var3\/$var2.group.picard.bam TMP_DIR=$scratch

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
     -T RealignerTargetCreator \
     -R $var1.fasta \
     -I $var3\/$var2.group.picard.bam \
     -o $var3\/$var2.intervals \

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
     -T IndelRealigner \
     -R $var1.fasta \
     -I $var3\/$var2.group.picard.bam \
     -targetIntervals $var3\/$var2.intervals \
     -o $var3\/$var2.group.picard.realigned.bam \

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
     -T ReduceReads \
     -R $var1.fasta \
     -I $var3\/$var2.group.picard.realigned.bam \
     -o $var3\/$var2.group.picard.reduced.bam \

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
     -T HaplotypeCaller \
     -R $var1.fasta \
     -I $var3\/$var2.group.picard.reduced.bam \
     -o $var3\/$var2.raw_variants.vcf \
     -nct 12
