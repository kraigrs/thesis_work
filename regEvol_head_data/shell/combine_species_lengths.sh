#!/bin/sh

# mel_mel


for ref in zhr_AMBIG_all z30_AMBIG_all droSim1_snp_indel droSec1_snp_indel sim_extra sec_extra
do

  var1=../../*/*.$ref.bowtie_v0_m1*.bed var2=$ref qsub concatenate2.sh

#  for length in 76 63 50 37
#  do

#    var1=../../*/*.$ref.bowtie_v0_m1_$length.bed var2=$ref var3=$length qsub concatenate.sh

#  done
done