#!/bin/sh

# scripts used for head data

####################################################################################################
# build indices for reference genomes

#sh build_ebwt.sh

####################################################################################################
# align all reads (these files reference a script with "trim" in the title, but all trimming has been removed)

#sh align_reads_mel_mel.sh
#sh align_reads_sim_sec.sh
#sh align_reads_mel_sec.sh
#sh align_reads_mel_sim.sh

####################################################################################################
# convert SAM to BAM

#for genome in zhr_AMBIG_all z30_AMBIG_all droSim1_snp_indel sim_extra droSec1_snp_indel sec_extra
#do 
#  FILES=../../*/*$genome*.sam
#  for sam in $FILES
#  do
#    var1=$sam var2=../../genomes/$genome.fa qsub samtools_view.sh
#  done
#done

####################################################################################################
# convert BAM to BED

#FILES=../../*/*.bam
#for bam in $FILES
#do
#  var1=$bam qsub bamToBed.sh
#done

####################################################################################################
# lift reads to correct genomic space

#sh lift_reads.sh

####################################################################################################
# combine droSim1 with sim_extra and droSec1 with sec_extra

#sh combine.sh

####################################################################################################
# convert chromosomal coordinates to exon coordinates

#for genome in zhr_AMBIG_all z30_AMBIG_all tsim sec
#do
#  FILES=../../*/*.$genome.*.bed.lifted
#  for f in $FILES
#  do
#    var1=$f var2=../../genomes/Sim_sec_gaps_gt19.bed var3=../../genomes/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.sh
#  done
#done

####################################################################################################
# characterize each paired-end read's allele-specific information

#sh characterize.sh

####################################################################################################
# quantify gene- and exon-levels of ASE

#sh classify_reads.sh


