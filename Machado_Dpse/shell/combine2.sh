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
        #cat ../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.droSim1_snp_indel.bowtie_v0_m1_all.bed.lifted.converted ../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.sim_extra.bowtie_v0_m1_all.bed.lifted.converted > ../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.tsim.bowtie_v0_m1_all.bed.lifted.converted

        #cat ../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.droSec1_snp_indel.bowtie_v0_m1_all.bed.lifted.converted ../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.sec_extra.bowtie_v0_m1_all.bed.lifted.converted > ../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.sec.bowtie_v0_m1_all.bed.lifted.converted

        cat ../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.droSim1_snp_indel.bowtie_v0_m1_all.bed.lifted ../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.sim_extra.bowtie_v0_m1_all.bed.lifted > ../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.tsim.bowtie_v0_m1_all.bed.lifted

        cat ../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.droSec1_snp_indel.bowtie_v0_m1_all.bed.lifted ../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.sec_extra.bowtie_v0_m1_all.bed.lifted > ../../sim_sec_data/$comp\_$sex\_r$rep\_m$mate.sec.bowtie_v0_m1_all.bed.lifted
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
        #cat ../../mel_sim_data/$comp\_$sex\_r$rep\_m$mate.droSim1_snp_indel.bowtie_v0_m1_all.bed.lifted.converted ../../mel_sim_data/$comp\_$sex\_r$rep\_m$mate.sim_extra.bowtie_v0_m1_all.bed.lifted.converted > ../../mel_sim_data/$comp\_$sex\_r$rep\_m$mate.tsim.bowtie_v0_m1_all.bed.lifted.converted

        cat ../../mel_sim_data/$comp\_$sex\_r$rep\_m$mate.droSim1_snp_indel.bowtie_v0_m1_all.bed.lifted ../../mel_sim_data/$comp\_$sex\_r$rep\_m$mate.sim_extra.bowtie_v0_m1_all.bed.lifted > ../../mel_sim_data/$comp\_$sex\_r$rep\_m$mate.tsim.bowtie_v0_m1_all.bed.lifted
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
        #cat ../../mel_sec_data/$comp\_$sex\_r$rep\_m$mate.droSec1_snp_indel.bowtie_v0_m1_all.bed.lifted.converted ../../mel_sec_data/$comp\_$sex\_r$rep\_m$mate.sec_extra.bowtie_v0_m1_all.bed.lifted.converted > ../../mel_sec_data/$comp\_$sex\_r$rep\_m$mate.sec.bowtie_v0_m1_all.bed.lifted.converted

        cat ../../mel_sec_data/$comp\_$sex\_r$rep\_m$mate.droSec1_snp_indel.bowtie_v0_m1_all.bed.lifted ../../mel_sec_data/$comp\_$sex\_r$rep\_m$mate.sec_extra.bowtie_v0_m1_all.bed.lifted > ../../mel_sec_data/$comp\_$sex\_r$rep\_m$mate.sec.bowtie_v0_m1_all.bed.lifted
      done
    done
  done
done