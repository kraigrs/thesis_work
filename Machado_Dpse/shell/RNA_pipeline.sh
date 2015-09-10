#!/bin/sh

# perform RNA alignments using iterative 13b trimming

#FILES=../../data/*carcass*.fastq
#for file in $FILES
#do
#  var1=../../genomes/dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.fasta var2=../../genomes/dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.fasta var3=$file qsub bowtie_pipeline.sh
#done

#FILES=../../data/*testes*.fastq
#for file in $FILES
#do
#  var1=../../genomes/dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.fasta var2=../../genomes/dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.fasta var3=$file qsub bowtie_pipeline.sh
#done

#FILES=../../data/*ovaries*.fastq
#for file in $FILES
#do
#  var1=../../genomes/dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.fasta var2=../../genomes/dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.fasta var3=$file qsub bowtie_pipeline.sh
#done

# convert SAM to BAM

#FILES=../../data/*Dpse_TL*.sam
#for sam in $FILES
#do
#  var1=$sam var2=../../genomes/dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.fasta qsub samtools_view.sh
#done

#FILES=../../data/*Dbog_Toro1*.sam
#for sam in $FILES
#do
#  var1=$sam var2=../../genomes/dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.fasta qsub samtools_view.sh
#done

# merge BAM files

#for sample in H5 H6 TL Toro1
#do
#  for tissue in F_carcass M_carcass ovaries testes
#  do
#    for mate in R1 R2 
#    do
#      for species in Dpse_TL Dbog_Toro1
#      do
#        var1=../../data/$sample\_$tissue\_$mate\.dpse_r3.1.$species.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bam var2=../../data/$sample\_$tissue\_$mate\.*$species*.bam qsub samtools_cat.sh
#      done
#    done
#  done
#done

# convert BAM to BED

#FILES=../../data/*merged.bam
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

#for sample in H5 H6 TL Toro1
#do
#  for tissue in F_carcass M_carcass ovaries testes
#  do
#    var1=../../data/$sample\_$tissue\_R1.dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bed.converted \
#    var2=../../data/$sample\_$tissue\_R1.dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bed.converted \
#    var3=../../data/$sample\_$tissue\_R2.dpse_r3.1.Dpse_TL.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bed.converted \
#    var4=../../data/$sample\_$tissue\_R2.dpse_r3.1.Dbog_Toro1.BQSR.exons_fsa_masked.bowtie_v0_m1.merged.bed.converted \
#    var5=../../data/$sample\_$tissue qsub Characterize_Each_Read.sh
#  done
#done

for sample in H5 H6 TL Toro1
do
  for tissue in F_carcass M_carcass ovaries testes
  do
    var1=../../genomes/dpse_r3.1.exons.bed \
    var2=Dpse_TL \
    var3=Dbog_Toro1 \
    var4=../../data/$sample\_$tissue.final_info_rmDup.txt \
    var5=../../data/$sample\_$tissue qsub classify.sh
  done
done