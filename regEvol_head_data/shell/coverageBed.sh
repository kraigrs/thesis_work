#!/bin/sh

#PBS -S /bin/sh
#PBS -N coverageBed
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=1,walltime=2:00:00
#PBS -l pmem=2000mb
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR
echo FILE = $var1

coverageBed -a $var1 -b $var2 > ../../genomes/lifted/$var3.coverage.txt

#coverageBed -a ../../genomes/lifted/zhr_AMBIG_all.bowtie_v0_m1.bed -b ../../genomes/lifted/constitutive_regions_overlap_filtered2.zhr.bed > ../../genomes/lifted/zhr.coverage.txt