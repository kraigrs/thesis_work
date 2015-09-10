#!/bin/sh

# mel_mel

#for comp in zhrXz30 z30Xzhr zhr_z30
#do
#  for sex in F M
#  do
#    for (( rep = 1; rep <= 3; rep++ )) 
#    do
#      var1=../../mel_mel_data/$comp\_$sex\_r$rep var2=zhr_AMBIG_all var3=z30_AMBIG_all qsub Characterize_Each_Read.sh
#    done
#  done
#done

# sim_sec

#for comp in simXsec sim_sec
#do
#  for sex in F M
#  do
#    for (( rep = 1; rep <= 3; rep++ )) 
#    do
#        var1=../../sim_sec_data/$comp\_$sex\_r$rep var2=tsim var3=sec qsub Characterize_Each_Read.sh
#    done
#  done
#done

#var1=../../sim_sec_data/simXsec_F_r1 var2=tsim var3=sec qsub Characterize_Each_Read.sh

# mel_sim

#for comp in simXzhr zhr_sim
#do
#  for sex in F M
#  do
#    for (( rep = 1; rep <= 3; rep++ )) 
#    do
#        var1=../../mel_sim_data/$comp\_$sex\_r$rep var2=zhr_AMBIG_all var3=tsim qsub Characterize_Each_Read.sh
#    done
#  done
#done

# mel_sec

for comp in zhrXsec zhr_sec
do
  for sex in F
  do
    for (( rep = 1; rep <= 3; rep++ )) 
    do
        var1=../../mel_sec_data/$comp\_$sex\_r$rep var2=zhr_AMBIG_all var3=sec qsub Characterize_Each_Read.sh
    done
  done
done