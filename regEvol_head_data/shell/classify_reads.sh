#!/bin/sh

# mel_mel

for comp in zhrXz30 z30Xzhr zhr_z30
do
  for sex in F M
  do
    for (( rep = 1; rep <= 3; rep++ )) 
    do
      var1=../../../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=z30 var4=../../mel_mel_data/$comp\_$sex\_r$rep.final_info_rmDup.txt var5=../../mel_mel_data/$comp\_$sex\_r$rep qsub classify.sh
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
      var1=../../../McManus/constitutive_regions_overlap_filtered.bed var2=tsim var3=sec var4=../../sim_sec_data/$comp\_$sex\_r$rep.final_info_rmDup.txt var5=../../sim_sec_data/$comp\_$sex\_r$rep qsub classify.sh
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
      var1=../../../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=tsim var4=../../mel_sim_data/$comp\_$sex\_r$rep.final_info_rmDup.txt var5=../../mel_sim_data/$comp\_$sex\_r$rep qsub classify.sh
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
      var1=../../../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=sec var4=../../mel_sec_data/$comp\_$sex\_r$rep.final_info_rmDup.txt var5=../../mel_sec_data/$comp\_$sex\_r$rep qsub classify.sh
    done
  done
done