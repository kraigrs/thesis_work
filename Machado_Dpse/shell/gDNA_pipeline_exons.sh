#!/bin/sh

# perform gDNA alignments using iterative 13b trimming

#var1=../../genomes/dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.fasta var2=../../genomes/dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.fasta var3=../../data/Dpse_TL.fastq qsub bowtie_pipeline.sh
#var1=../../genomes/dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.fasta var2=../../genomes/dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.fasta var3=../../data/Dbog_Toro1.fastq qsub bowtie_pipeline.sh

#var1=../../genomes/dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.fasta var2=../../genomes/dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.fasta var3=../../data/Dbog_Toro1_R1.fastq qsub bowtie_pipeline.sh
#var1=../../genomes/dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.fasta var2=../../genomes/dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.fasta var3=../../data/Dbog_Toro1_R2.fastq qsub bowtie_pipeline.sh

# convert SAM to BAM

#FILES=../../data/*Dpse_TL.BQSR.exons_fsa_masked*.sam
#for sam in $FILES
#do
#  var1=$sam var2=../../genomes/dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.fasta qsub samtools_view.sh
#done

#FILES=../../data/*Dbog_Toro1.BQSR.exons_fsa_masked*.sam
#for sam in $FILES
#do
#  var1=$sam var2=../../genomes/dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.fasta qsub samtools_view.sh
#done

# merge BAM files

#for species in Dpse_TL Dbog_Toro1 
#do
#  var1=../../data/$species.dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bam var2=../../data/$species.dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.bowtie_v0_m1*.bam qsub samtools_cat.sh
#  var1=../../data/$species.dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bam var2=../../data/$species.dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.bowtie_v0_m1*.bam qsub samtools_cat.sh
#done

# sort BAM files

#FILES=../../data/*merged.bam
#for bam in $FILES
#do
#  var1=$bam qsub samtools_sort.sh
#done

# remove PCR duplicates

#FILES=../../data/*sorted.bam
#for bam in $FILES
#do
#  var1=$bam qsub samtools_rmdup.sh
#done

# convert BAM to BED

#FILES=../../data/*BQSR*merged.bam
#for bam in $FILES
#do
#  var1=$bam qsub bamToBed.sh
#done

# reformat BED

#FILES=../../data/*merged.bed
#for bed in $FILES
#do
#  var1=$bed qsub reformatBed.sh
#done

#var1=../../data/Dbog_Toro1_R1.dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bed.converted \
#var2=../../data/Dbog_Toro1_R1.dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bed.converted \
#var3=../../data/Dbog_Toro1_R2.dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bed.converted \
#var4=../../data/Dbog_Toro1_R2.dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bed.converted \
#var5=../../data/Dbog_Toro1 qsub Characterize_Each_Read.sh

#var1=../../genomes/dpse_r3.1.exons.bed \
#var2=Dpse_TL \
#var3=Dbog_Toro1 \
#var4=../../data/Dbog_Toro1.final_info_rmDup.txt \
#var5=../../data/Dbog_Toro1 qsub classify.sh



# treat all data as single-end

#cd ../single_end/shell

#var1=../../../data/Dpse_TL.dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bed.converted \
#var2=../../../data/Dpse_TL.dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bed.converted \
#var3=../../../data/Dpse_TL qsub Characterize_Each_Read.sh

#cd ../single_end/shell

#var1=../../../data/Dbog_Toro1.dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bed.converted \
#var2=../../../data/Dbog_Toro1.dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bed.converted \
#var3=../../../data/Dbog_Toro1 qsub Characterize_Each_Read.sh



var1=../../genomes/dpse_r3.1.exons.bed \
var2=Dpse_TL \
var3=Dbog_Toro1 \
var4=../../data/Dpse_TL.final_info_rmDup.txt \
var5=../../data/Dpse_TL qsub classify_SE.sh

var1=../../genomes/dpse_r3.1.exons.bed \
var2=Dpse_TL \
var3=Dbog_Toro1 \
var4=../../data/Dbog_Toro1.final_info_rmDup.txt \
var5=../../data/Dbog_Toro1 qsub classify_SE.sh




#var1=../../data/Dpse_TL.dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.sorted.rmdup.bed var2=../../genomes/dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.bed var3=../../data/Dpse_TL.site_coverage.txt qsub coverageBed.sh
#var1=../../data/Dbog_Toro1.dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.sorted.rmdup.bed var2=../../genomes/dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.bed var3=../../data/Dbog_Toro1.site_coverage.txt qsub coverageBed.sh

#var1=../../data/Dpse_TL.dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.sorted.rmdup.bed var2=../../genomes/dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.bed var3=../../data/Dpse_TL.exon_coverage.txt qsub coverageBed.sh
#var1=../../data/Dbog_Toro1.dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.sorted.rmdup.bed var2=../../genomes/dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.bed var3=../../data/Dbog_Toro1.exon_coverage.txt qsub coverageBed.sh