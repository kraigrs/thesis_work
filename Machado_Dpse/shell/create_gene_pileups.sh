#!/bin/sh

# perform gDNA alignments using iterative 13b trimming

# make sure to specify --phred64

#var1=../../genomes/dpse-all-gene-r3.1.Dpse_TL_samtools.fasta var2=../../data/Dpse_TL.fastq qsub bowtie2_SE_pipeline.sh
#var1=../../genomes/dpse-all-gene-r3.1.Dpse_TL_GATK.fasta var2=../../data/Dpse_TL.fastq qsub bowtie2_SE_pipeline.sh

#var1=../../genomes/dpse-all-gene-r3.1.Dbog_Toro1_samtools.fasta var2=../../data/Dbog_Toro1_R1.fastq var3=../../data/Dbog_Toro1_R2.fastq qsub bowtie2_PE_pipeline.sh
#var1=../../genomes/dpse-all-gene-r3.1.Dbog_Toro1_GATK.fasta var2=../../data/Dbog_Toro1_R1.fastq var3=../../data/Dbog_Toro1_R2.fastq qsub bowtie2_PE_pipeline.sh

# convert SAM to BAM

FILES=../../data/*.sam
for sam in $FILES
do
  var1=$sam var2=../../genomes/dpse-all-gene-r3.1.fasta qsub samtools_view.sh
done

# merge BAM files

#for species in Dpse_TL Dbog_Toro1
#do
#  samtools cat -o ../../data/$species.dpse-all-gene-r3.1.$species_samtools.bowtie2_very_sensitive.merged.bam ../../data/$species.dpse-all-gene-r3.1.$species_samtools*.bam
#  samtools cat -o ../../data/$species.dpse-all-gene-r3.1.$species_GATK.bowtie2_very_sensitive.merged.bam ../../data/$species.dpse-all-gene-r3.1.$species_GATK*.bam
#done

# sort BAM files

#for species in Dpse_TL Dbog_Toro1
#do
#  var1=../../data/$species.dpse-all-gene-r3.1.$species_samtools.bowtie2_very_sensitive.merged.bam var2=../../data/$species.dpse-all-gene-r3.1.$species_samtools.bowtie2_very_sensitive.merged.sorted qsub samtools_sort.sh
#  var1=../../data/$species.dpse-all-gene-r3.1.$species_GATK.bowtie2_very_sensitive.merged.bam var2=../../data/$species.dpse-all-gene-r3.1.$species_GATK.bowtie2_very_sensitive.merged.sorted qsub samtools_sort.sh
#done

# remove PCR duplicates

#qsub samtools_rmdup.sh

# create mpileups

#var1=../../genomes/dpse-all-gene-r3.1.Dpse_TL_samtools var2=Dpse_TL.dpse-all-gene-r3.1.Dpse_TL_samtools.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub samtools_mpileup.sh
#var1=../../genomes/dpse-all-gene-r3.1.Dpse_TL_GATK var2=Dpse_TL.dpse-all-gene-r3.1.Dpse_TL_GATK.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub samtools_mpileup.sh

#var1=../../genomes/dpse-all-gene-r3.1.Dbog_Toro1_samtools var2=Dbog_Toro1.dpse-all-gene-r3.1.Dbog_Toro1_samtools.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub samtools_mpileup.sh
#var1=../../genomes/dpse-all-gene-r3.1.Dbog_Toro1_GATK var2=Dbog_Toro1.dpse-all-gene-r3.1.Dbog_Toro1_GATK.bowtie2_very_sensitive.merged.sorted.rmdup var3=../../data qsub samtools_mpileup.sh