#!/bin/sh

#PBS -S /bin/sh
#PBS -N pairwise_aln_FSA
#PBS -A wittkopp_flux
#PBS -l qos=flux
#PBS -l nodes=1:ppn=12,pmem=4000mb,walltime=10:00:00
#PBS -q flux
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M kraigrs@umich.edu

####  End PBS preamble
#  Include the next three lines always

if [ "x${PBS_NODEFILE}" != "x" ] ; then
   cat $PBS_NODEFILE   # contains a list of the CPUs you were using if run with PBS
fi

cd $PBS_O_WORKDIR

#Usage: var1=<file1>... multiple_aln_FSA.sh

perl ../perl/pairwise_aln_FSA.pl $var1 $var2 $var3 $var4

#var1=zhr var2=z30 var3=../../references/zhr_z30/dm3_regulatory_regions.zhr.fa var4=../../references/zhr_z30/dm3_regulatory_regions.z30.fa qsub pairwise_aln_FSA.sh
#var1=tsim var2=droSec1 var3=../../references/tsim_droSec1/dm3_regulatory_regions.tsim.fa var4=../../references/tsim_droSec1/dm3_regulatory_regions.droSec1.fa qsub pairwise_aln_FSA.sh
#var1=zhr var2=tsim var3=../../references/zhr_tsim/dm3_regulatory_regions.zhr.fa var4=../../references/zhr_tsim/dm3_regulatory_regions.tsim.fa qsub pairwise_aln_FSA.sh






#var1=zhr var2=z30 var3=../../references/zhr_z30/dm3_regulatory_regions_merged.zhr.fa var4=../../references/zhr_z30/dm3_regulatory_regions_merged.z30.fa qsub pairwise_aln_FSA.sh

#var1=tsim var1=droSec1 var3=../../references/tsim_droSec1/dm3_regulatory_regions_merged.tsim.fa var4=../../references/tsim_droSec1/dm3_regulatory_regions_merged.droSec1.fa qsub pairwise_aln_FSA.sh

#var1=zhr var2=tsim var3=../../references/zhr_tsim/dm3_regulatory_regions_merged.zhr.fa var4=../../references/zhr_tsim/dm3_regulatory_regions_merged.tsim.fa qsub pairwise_aln_FSA.sh