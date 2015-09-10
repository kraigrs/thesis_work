#!/bin/sh

# zhr

FILES=../../*/*zhr_AMBIG_all*merged.bed
for f in $FILES
do
  var1=$f var2=../../genomes/zhr/zhr_AMBIG_all_dm3.chain qsub liftOver.sh
done

# z30

FILES=../../*/*z30_AMBIG_all*merged.bed
for f in $FILES
do
  var1=$f var2=../../genomes/z30/z30_AMBIG_all_dm3.chain qsub liftOver.sh
done

# droSec1

FILES=../../*/*droSec1_snp_indel*merged.bed
for f in $FILES
do
  var1=$f var2=../../genomes/droSec1/droSec1_snp_indel_to_droSec1.over.chain var3=../../genomes/droSec1/droSec1ToDm3.over.chain qsub liftOver2.sh
done

# sec extra

FILES=../../*/*sec_extra*merged.bed
for f in $FILES
do
  var1=$f var2=../../genomes/droSec1/Sec_extra_To_dm3.over.chain qsub liftOver.sh
done

# tsim

FILES=../../*/*droSim1_snp_indel*merged.bed
for f in $FILES
do
  var1=$f var2=../../genomes/tsimbazaza/droSim1_snp_indel_to_droSim1.over.chain var3=../../genomes/tsimbazaza/droSim1ToDm3.over.chain qsub liftOver2.sh
done

# sim extra

FILES=../../*/*sim_extra*merged.bed
for f in $FILES
do
  var1=$f var2=../../genomes/tsimbazaza/Sim_extra_To_dm3.over.chain qsub liftOver.sh
done