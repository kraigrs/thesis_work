#!/bin/sh

# mel_mel

for comp in zhrXz30 z30Xzhr zhr_z30
do
  for sex in F M
  do
    for (( rep = 1; rep <= 3; rep++ )) 
    do
      for (( mate = 1; mate <= 2; mate++ )) 
      do
        for ref in zhr_AMBIG_all z30_AMBIG_all
        do
          var1=../../mel_mel_data/$comp\_$sex\_r$rep\_m$mate.$ref qsub concatenate.sh
          #ls ../../mel_mel_data/$comp\_$sex\_r$rep\_m$mate.$ref.bowtie_v0_m1_all.bed.lifted.converted
        done
      done
    done
  done
done

# sim_sec

for comp in simXsec sim_sec
do
  for sex in F M
  do
    for (( rep = 1; rep <= 3; rep++ )) 
    do
      for (( mate = 1; mate <= 2; mate++ )) 
      do
        for ref in droSim1_snp_indel droSec1_snp_indel sim_extra sec_extra
        do
          var1=../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.$ref qsub concatenate.sh
          #ls ../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.$ref.bowtie_v0_m1_all.bed.lifted.converted
        done
      done
    done
  done
done

# mel_sim

for comp in simXzhr zhr_sim
do
  for sex in F M
  do
    for (( rep = 1; rep <= 3; rep++ )) 
    do
      for (( mate = 1; mate <= 2; mate++ )) 
      do
        for ref in droSim1_snp_indel sim_extra zhr_AMBIG_all
        do
          var1=../../mel_sim_data/$comp\_$sex\_r$rep\_m$mate.$ref qsub concatenate.sh
          #ls ../../mel_sim_data/$comp\_$sex\_r$rep\_m$mate.$ref.bowtie_v0_m1_all.bed.lifted.converted
        done
      done
    done
  done
done

# mel_sec

for comp in zhrXsec zhr_sec
do
  for sex in F
  do
    for (( rep = 1; rep <= 3; rep++ )) 
    do
      for (( mate = 1; mate <= 2; mate++ )) 
      do
        for ref in droSec1_snp_indel sec_extra zhr_AMBIG_all
        do
          var1=../../mel_sec_data/$comp\_$sex\_r$rep\_m$mate.$ref qsub concatenate.sh
          #ls ../../mel_sec_data/$comp\_$sex\_r$rep\_m$mate.$ref.bowtie_v0_m1_all.bed.lifted.converted
        done
      done
    done
  done
done