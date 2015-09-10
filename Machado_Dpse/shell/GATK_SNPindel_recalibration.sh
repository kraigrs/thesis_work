#!/bin/sh

#PBS -S /bin/sh
#PBS -N GATK_SNPindel_recalibration
#PBS -A wittkopp_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l nodes=1:ppn=12,pmem=4000mb,walltime=10:00:00
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M kraigrs@umich.edu

####  End PBS preamble
#  Include the next three lines always

if [ "x${PBS_NODEFILE}" != "x" ] ; then
   cat $PBS_NODEFILE   # contains a list of the CPUs you were using if run with PBS
fi

cd $PBS_O_WORKDIR

#Usage: var1=<ref>.[fa/fasta] var2=<file>.bam var3=<data directory> qsub GATK_SNPindel_recalibration.sh

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
     -T BaseRecalibrator \
     -R $var1.fasta \
     -I $var3\/$var2.group.picard.realigned.bam \
     -knownSites $var3\/$var2.raw_SNPs_highQual.vcf \
     -o $var3\/$var2.recal.table \
     -nct 12

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
     -T AnalyzeCovariates \
     -R $var1.fasta \
     -BQSR $var3\/$var2.recal.table \
     -plots $var3\/$var2.BQSR_recal.pdf

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
     -T PrintReads \
     -R $var1.fasta \
     -I $var3\/$var2.group.picard.realigned.bam \
     -BQSR $var3\/$var2.recal.table \
     -o $var3\/$var2.group.picard.recal.bam \
     -nct 12

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
     -T ReduceReads \
     -R $var1.fasta \
     -I $var3\/$var2.group.picard.recal.bam \
     -o $var3\/$var2.group.picard.recal.reduced.bam \

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
     -T HaplotypeCaller \
     -R $var1.fasta \
     -I $var3\/$var2.group.picard.recal.reduced.bam \
     -o $var3\/$var2.BQSR_variants.vcf \
     -nct 12
