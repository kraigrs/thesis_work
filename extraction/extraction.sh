#!/bin/sh

# for a set of dm3 regions, lift coordinates to 
# zhr, z30, tsim, and sec genomes and extract their relevant sequences

if [[ $1 =~ ([[:print:]]*).bed$ ]]
then
  file=${BASH_REMATCH[1]}
  #echo $1
else
  echo "***Error*** $1 appears not to have correct .bed extension"
  exit 1
fi

if [ $# != 2 ]
then
  prop=0.95
  #echo $prop
else
  prop=$2
  #echo $2
fi

# zhr

liftOver $file.bed dm3ToZhr_AMBIG_all.chain $file.zhr.bed $file.zhr.un -minMatch=$prop
perl get_sequences.pl $file.zhr.bed zhr_AMBIG_all.fa $file.zhr.fa

# z30

liftOver $file.bed dm3ToZ30_AMBIG_all.chain $file.z30.bed $file.z30.un -minMatch=$prop
perl get_sequences.pl $file.z30.bed z30_AMBIG_all.fa $file.z30.fa

# tsim - this needs 3 steps because there are 3 genomes: droSim1_snp_indel.fa and sim_extra.fa

liftOver $file.bed dm3ToDroSim1.over.chain $file.droSim1.bed $file.droSim1.un -minMatch=$prop
liftOver $file.droSim1.bed droSim1ToDroSim1_snp_indel.over.chain $file.droSim1_snp_indel.bed $file.droSim1_snp_indel.un -minMatch=$prop
perl get_sequences.pl $file.droSim1_snp_indel.bed droSim1_snp_indel.fa $file.droSim1_snp_indel.fa

liftOver $file.droSim1.un dm3ToSim_extra.over.chain $file.sim_extra.bed $file.sim_extra.un -minMatch=$prop
perl get_sequences.pl $file.sim_extra.bed sim_extra.fa $file.sim_extra.fa

cat $file.droSim1_snp_indel.fa $file.sim_extra.fa > $file.tsim.fa
rm $file.droSim1_snp_indel.fa $file.sim_extra.fa

# sec - this needs 3 steps because there are 3 genomes: droSec1_snp_indel.fa and sec_extra.fa

liftOver $file.bed dm3ToDroSec1.over.chain $file.droSec1.bed $file.droSec1.un -minMatch=$prop
liftOver $file.droSec1.bed droSec1ToDroSec1_snp_indel.over.chain $file.droSec1_snp_indel.bed $file.droSec1_snp_indel.un -minMatch=$prop
perl get_sequences.pl $file.droSec1_snp_indel.bed droSec1_snp_indel.fa $file.droSec1_snp_indel.fa

liftOver $file.droSec1.un dm3ToSec_extra.over.chain $file.sec_extra.bed $file.sec_extra.un -minMatch=$prop
perl get_sequences.pl $file.sec_extra.bed sec_extra.fa $file.sec_extra.fa

cat $file.droSec1_snp_indel.fa $file.sec_extra.fa > $file.sec.fa
rm $file.droSec1_snp_indel.fa $file.sec_extra.fa
