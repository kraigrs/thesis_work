#!/bin/sh

#PBS -S /bin/sh
#PBS -N GATK_SNPindel_filtration
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=1,pmem=4000mb,walltime=10:00:00
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

#Usage: var1=<ref>.[fa/fasta] var2=<file>.bam var3=<data directory> qsub samtools_SNPindel.sh

#http://gatkforums.broadinstitute.org/discussion/2806/howto-apply-hard-filters-to-a-call-set
#http://gatkforums.broadinstitute.org/discussion/3225/how-can-i-filter-my-callset-if-i-cannot-use-vqsr-recalibrate-variants

#java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
#    -T SelectVariants \
#    -R $var1.fasta \
#    -V $var3\/$var2.raw_variants.vcf \
#    -selectType SNP \
#    -o $var3\/$var2.raw_SNPs.vcf

#java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
#    -T VariantFiltration \
#    -R $var1.fasta \
#    --variant $var3\/$var2.raw_SNPs.vcf \
#    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \
#    --filterName "FAIL" \
#    -o $var3\/$var2.raw_SNPs_filtered.vcf

#grep -v FAIL $var3\/$var2.raw_SNPs_filtered.vcf > $var3\/$var2.raw_SNPs_highQual.vcf

#java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
#    -T SelectVariants \
#    -R $var1.fasta \
#    -V $var3\/$var2.BQSR_variants.vcf \
#    -selectType SNP \
#    -o $var3\/$var2.BQSR_SNPs.vcf

#java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
#    -T VariantFiltration \
#    -R $var1.fasta \
#    --variant $var3\/$var2.BQSR_SNPs.vcf \
#    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \
#    --filterName "FAIL" \
#    -o $var3\/$var2.BQSR_SNPs_filtered.vcf

#grep -v FAIL $var3\/$var2.BQSR_SNPs_filtered.vcf > $var3\/$var2.BQSR_SNPs_highQual.vcf

# run this after first iteration and add +1 to iter[n]

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R $var1.fasta \
    -V $var3\/$var2.BQSR_variants_iter2.vcf \
    -selectType SNP \
    -o $var3\/$var2.BQSR_SNPs_iter2.vcf

java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R $var1.fasta \
    --variant $var3\/$var2.BQSR_SNPs_iter2.vcf \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "FAIL" \
    -o $var3\/$var2.BQSR_SNPs_iter2_filtered.vcf

grep -v FAIL $var3\/$var2.BQSR_SNPs_iter2_filtered.vcf > $var3\/$var2.BQSR_SNPs_iter2_highQual.vcf

# QD = variant confidence from the QUAL field divided by the unfiltered depth of non-reference samples
# FS = Phred-scaled p-value using Fisher’s Exact Test to detect strand bias (more bias --> false positives)
# MQ = Root Mean Square of the mapping quality of the reads across all samples
# HaplotypeScore = consistency of the site with two (and only two) segregating haplotypes
# MQRankSum = u-based z-approximation from the Mann-Whitney Rank Sum Test for mapping qualities 
#             (reads with ref bases vs. those with the alternate allele)
# ReadPosRankSum = u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele

#java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
#    -T SelectVariants \
#    -R $var1.fasta \
#    -V $var3\/$var2.raw_variants.vcf \
#    -L $var1.genes.bed \
#    -selectType INDEL \
#    -o $var3\/$var2.raw_indels.vcf

#java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
#    -T VariantFiltration \
#    -R $var1.fasta \
#    --variant $var3\/$var2.raw_indels.vcf \
#    --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
#    --filterName "FAIL" \
#    -o $var3\/$var2.GATK_indels.vcf

#java -Xmx4g -jar $GATK_JAVA/GenomeAnalysisTK.jar \
#    -T CombineVariants \
#    -R $var1.fasta \
#    --variant $var3\/$var2.GATK_SNPs.vcf \
#    --variant $var3\/$var2.GATK_indels.vcf \
#    -o $var3\/$var2.GATK_variants.vcf
