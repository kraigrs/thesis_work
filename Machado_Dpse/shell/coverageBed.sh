#!/bin/sh

#PBS -S /bin/sh
#PBS -N coverageBed
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=1,walltime=1:00:00
#PBS -m abe
#PBS -j oe
#PBS -V

####  End PBS preamble
#  Include the next three lines always

if [ "x${PBS_NODEFILE}" != "x" ] ; then
   cat $PBS_NODEFILE   # contains a list of the CPUs you were using if run with PBS
fi

cd $PBS_O_WORKDIR

echo $var1

#coverageBed -a $var1 -b $var2 -d > $var3
coverageBed -a $var1 -b $var2 > $var3

#coverageBed -a ../../genomes/lifted/zhr_AMBIG_all.bowtie_v0_m1.bed -b ../../genomes/lifted/constitutive_regions_overlap_filtered2.zhr.bed > ../../genomes/lifted/zhr.coverage.txt