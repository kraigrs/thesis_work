#!/bin/sh

# sim_sec

for comp in simXsec sim_sec
do
  for sex in F M
  do
    for (( rep = 1; rep <= 3; rep++ )) 
    do
      for (( mate = 1; mate <= 2; mate++ )) 
      do
        var1=../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.droSim1_snp_indel.bowtie_v0_m1.merged.bed.lifted var2=../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.sim_extra.bowtie_v0_m1.merged.bed.lifted var3=../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.tsim.bowtie_v0_m1.merged.bed.lifted qsub concatenate.sh

        var1=../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.droSec1_snp_indel.bowtie_v0_m1.merged.bed.lifted var2=../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.sec_extra.bowtie_v0_m1.merged.bed.lifted var3=../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.sec.bowtie_v0_m1.merged.bed.lifted qsub concatenate.sh
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
        var1=../../mel_sim_data/$comp\_$sex\_r$rep\_m$mate.droSim1_snp_indel.bowtie_v0_m1.merged.bed.lifted var2=../../mel_sim_data/$comp\_$sex\_r$rep\_m$mate.sim_extra.bowtie_v0_m1.merged.bed.lifted var3=../../mel_sim_data/$comp\_$sex\_r$rep\_m$mate.tsim.bowtie_v0_m1.merged.bed.lifted qsub concatenate.sh
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
        var1=../../mel_sec_data/$comp\_$sex\_r$rep\_m$mate.droSec1_snp_indel.bowtie_v0_m1.merged.bed.lifted var2=../../mel_sec_data/$comp\_$sex\_r$rep\_m$mate.sec_extra.bowtie_v0_m1.merged.bed.lifted var3=../../mel_sec_data/$comp\_$sex\_r$rep\_m$mate.sec.bowtie_v0_m1.merged.bed.lifted qsub concatenate.sh
      done
    done
  done
done