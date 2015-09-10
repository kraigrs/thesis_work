#!/bin/sh

# perform gDNA alignments using iterative 13b trimming

# make sure to specify --phred64

#var1=../../genomes/dpse_r3.1.fasta var2=../../data/Dpse_TL.fastq qsub bowtie2_SE_pipeline.sh
#var1=../../genomes/dpse_r3.1.fasta var2=../../data/Dbog_Toro1_R1.fastq var3=../../data/Dbog_Toro1_R2.fastq qsub bowtie2_PE_pipeline.sh

# treat everything as single-end

#var1=../../genomes/dpse_r3.1.fasta var2=../../data/Dpse_TL.fastq qsub bowtie2_SE_pipeline.sh
#var1=../../genomes/dpse_r3.1.fasta var2=../../data/Dbog_Toro1.fastq qsub bowtie2_SE_pipeline.sh

# convert SAM to BAM

#FILES=../../data/*.sam
#for sam in $FILES
#do
#  var1=$sam var2=../../genomes/dpse_r3.1.fasta qsub samtools_view.sh
#done

# merge BAM files

#for species in Dpse_TL Dbog_Toro1
#do
#  var1=../../data/$species.dpse_r3.1.bowtie2_very_sensitive.merged.bam var2=../../data/$species*.bam qsub samtools_cat.sh
#done

# sort BAM files

#for species in Dpse_TL Dbog_Toro1
#do
#  var1=../../data/$species.dpse_r3.1.bowtie2_very_sensitive.merged.bam qsub samtools_sort.sh
#done

# remove PCR duplicates

#for species in Dpse_TL Dbog_Toro1
#do
#  var1=../../data/$species.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.bam qsub samtools_rmdup.sh
#done

# used Picard to create a dictionary file:

#java -jar CreateSequenceDictionary.jar R= ../../../Wittkopp/Machado_Dpse/genomes/dpse_r3.1.fasta O= ../../../Wittkopp/Machado_Dpse/genomes/dpse_r3.1.dict

# create mpileups to determine avg depth of coverage

#var1=../../genomes/dpse_r3.1 var2=Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub samtools_mpileup.sh
#var1=../../genomes/dpse_r3.1 var2=Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub samtools_mpileup.sh

#var1=../../data/Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.mpileup.txt var2=../../data/Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.raw_variants.vcf var3=../../data/Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.raw_variants_highQual.vcf qsub filter_VCF.sh
#var1=../../data/Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.mpileup.txt var2=../../data/Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.raw_variants.vcf var3=../../data/Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.raw_variants_highQual.vcf qsub filter_VCF.sh

# call variants and indels

#var1=../../genomes/dpse_r3.1 var2=Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data var4=single var5=Dpse_TL qsub GATK_SNPindel.sh
#var1=../../genomes/dpse_r3.1 var2=Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data var4=single var5=Dbog_Toro1 qsub GATK_SNPindel.sh

# iterate these next steps
###############################

# isolate high-quality variants

#var1=../../genomes/dpse_r3.1 var2=Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub GATK_SNPindel_filtration.sh
#var1=../../genomes/dpse_r3.1 var2=Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub GATK_SNPindel_filtration.sh

# recalibrate based on high-quality variants

#var1=../../genomes/dpse_r3.1 var2=Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub GATK_SNPindel_recalibration.sh
#var1=../../genomes/dpse_r3.1 var2=Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub GATK_SNPindel_recalibration.sh

# iterate

#var1=../../genomes/dpse_r3.1 var2=Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub GATK_SNPindel_filtration.sh
#var1=../../genomes/dpse_r3.1 var2=Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub GATK_SNPindel_filtration.sh

#var1=../../genomes/dpse_r3.1 var2=Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub GATK_SNPindel_recalibration_iteration.sh
#var1=../../genomes/dpse_r3.1 var2=Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub GATK_SNPindel_recalibration_iteration.sh

#################################

# filter variants

#var1=../../genomes/dpse_r3.1 var2=Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub GATK_SNPindel_filtration.sh
#var1=../../genomes/dpse_r3.1 var2=Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub GATK_SNPindel_filtration.sh

# incorporate variants into annotated genes or exons

#perl convertVCF.pl ../../genomes/dpse_r3.1.fasta ../../genomes/dpse_r3.1.exons.bed ../../genomes/Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.BQSR_variants_iter3.vcf ../../genomes/Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.BQSR_variants_iter3.exons.vcf
#perl vcf2fasta.pl ../../genomes/dpse_r3.1.exons.fasta ../../genomes/Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.BQSR_variants_iter3.exons.vcf ../../genomes/dpse_r3.1.Dpse_TL.BQSR.exons.fasta

#perl convertVCF.pl ../../genomes/dpse_r3.1.fasta ../../genomes/dpse_r3.1.exons.bed ../../genomes/Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.BQSR_variants_iter3.vcf ../../genomes/Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.BQSR_variants_iter3.exons.vcf
#perl vcf2fasta.pl ../../genomes/dpse_r3.1.exons.fasta ../../genomes/Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.BQSR_variants_iter3.exons.vcf ../../genomes/dpse_r3.1.Dbog_Toro1.BQSR.exons.fasta

# pairwise alignments of Dpse_TL and Dbog_Toro1 exons

#var1=../../genomes/dpse_r3.1.Dpse_TL.BQSR.exons.fasta var2=../../genomes/dpse_r3.1.Dbog_Toro1.BQSR.exons.fasta qsub pairwise_aln_FSA.sh
#var1=../../genomes/dpse_r3.1.Dpse_TL.BQSR.exons_fsa.fasta var2=../../genomes/dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa.fasta var3=../../genomes/dpse_r3.1.Dpse_TL.Dbog_Toro1.BQSR.exons.SNPs.txt qsub compare_pairwise.sh