#!/bin/sh

# mel_mel

for sample in zhr_z30 zhrXz30 z30Xzhr
do
  for sex in F M
  do
    for rep in r1 r2 r3
    do
      for mate in m1 m2
      do
        for species in zhr_AMBIG_all z30_AMBIG_all
        do
          var1=../../mel_mel_data/$sample\_$sex\_$rep\_$mate\.$species.bowtie_v0_m1.merged.bam var2=../../mel_mel_data/$sample\_$sex\_$rep\_$mate\.$species.bowtie_v0_m1*.bam qsub samtools_cat.sh
        done
      done
    done
  done
done

# sim_sec

for sample in simXsec sim_sec
do
  for sex in F M
  do
    for rep in r1 r2 r3
    do
      for mate in m1 m2
      do
        for species in droSim1_snp_indel sim_extra droSec1_snp_indel sec_extra
        do
          var1=../../sim_sec_data/$sample\_$sex\_$rep\_$mate\.$species.bowtie_v0_m1.merged.bam var2=../../sim_sec_data/$sample\_$sex\_$rep\_$mate\.$species.bowtie_v0_m1*.bam qsub samtools_cat.sh
        done
      done
    done
  done
done

# mel_sec

for sample in zhr_sec zhrXsec
do
  for sex in F
  do
    for rep in r1 r2 r3
    do
      for mate in m1 m2
      do
        for species in zhr_AMBIG_all droSec1_snp_indel sec_extra
        do
          var1=../../mel_sec_data/$sample\_$sex\_$rep\_$mate\.$species.bowtie_v0_m1.merged.bam var2=../../mel_sec_data/$sample\_$sex\_$rep\_$mate\.$species.bowtie_v0_m1*.bam qsub samtools_cat.sh
        done
      done
    done
  done
done

# mel_sim

for sample in zhr_sim simXzhr
do
  for sex in F M
  do
    for rep in r1 r2 r3
    do
      for mate in m1 m2
      do
        for species in zhr_AMBIG_all droSim1_snp_indel sim_extra
        do
          var1=../../mel_sim_data/$sample\_$sex\_$rep\_$mate\.$species.bowtie_v0_m1.merged.bam var2=../../mel_sim_data/$sample\_$sex\_$rep\_$mate\.$species.bowtie_v0_m1*.bam qsub samtools_cat.sh
        done
      done
    done
  done
done