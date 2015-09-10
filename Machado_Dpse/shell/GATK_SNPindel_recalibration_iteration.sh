#!/bin/sh

#PBS -S /bin/sh
#PBS -N GATK_SNPindel_recalibration_iteration
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

#Usage: var1=<ref>.[fa/fasta] var2=<file>.bam var3=<data directory> qsub GATK_SNPindel_recalibration_iteration.sh

#java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
#     -T BaseRecalibrator \
#     -R $var1.fasta \
#     -I $var3\/$var2.group.picard.recal.bam \
#     -knownSites $var3\/$var2.BQSR_SNPs_highQual.vcf \
#     -o $var3\/$var2.recal_iter1.table \
#     -nct 12

#java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
#     -T AnalyzeCovariates \
#     -R $var1.fasta \
#     -before $var3\/$var2.recal.table \
#     -after $var3\/$var2.recal_iter1.table \
#     -plots $var3\/$var2.BQSR_recal_iter1.pdf

#java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
#     -T PrintReads \
#     -R $var1.fasta \
#     -I $var3\/$var2.group.picard.recal.bam \
#     -BQSR $var3\/$var2.recal_iter1.table \
#     -o $var3\/$var2.group.picard.recal_iter1.bam \
#     -nct 12

#java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
#     -T ReduceReads \
#     -R $var1.fasta \
#     -I $var3\/$var2.group.picard.recal_iter1.bam \
#     -o $var3\/$var2.group.picard.recal_iter1.reduced.bam \

#java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
#     -T HaplotypeCaller \
#     -R $var1.fasta \
#     -I $var3\/$var2.group.picard.recal_iter1.reduced.bam \
#     -o $var3\/$var2.BQSR_variants_iter1.vcf \
#     -nct 12

# run this after first iteration and add +1 to iter[n]

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
     -T BaseRecalibrator \
     -R $var1.fasta \
     -I $var3\/$var2.group.picard.recal_iter2.bam \
     -knownSites $var3\/$var2.BQSR_SNPs_iter2_highQual.vcf \
     -o $var3\/$var2.recal_iter3.table \
     -nct 12

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
     -T AnalyzeCovariates \
     -R $var1.fasta \
     -before $var3\/$var2.recal_iter2.table \
     -after $var3\/$var2.recal_iter3.table \
     -plots $var3\/$var2.BQSR_recal_iter3.pdf

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
     -T PrintReads \
     -R $var1.fasta \
     -I $var3\/$var2.group.picard.recal_iter2.bam \
     -BQSR $var3\/$var2.recal_iter3.table \
     -o $var3\/$var2.group.picard.recal_iter3.bam \
     -nct 12

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
     -T ReduceReads \
     -R $var1.fasta \
     -I $var3\/$var2.group.picard.recal_iter3.bam \
     -o $var3\/$var2.group.picard.recal_iter3.reduced.bam \

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
     -T HaplotypeCaller \
     -R $var1.fasta \
     -I $var3\/$var2.group.picard.recal_iter3.reduced.bam \
     -o $var3\/$var2.BQSR_variants_iter3.vcf \
     -nct 12
