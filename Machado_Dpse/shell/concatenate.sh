#!/bin/sh

#PBS -S /bin/sh
#PBS -N concatenate
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=1,walltime=1:00:00
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

#Usage: var1=<file1> qsub concatenate.sh

#cat $var1*.bed > $var1.bowtie_all.bed

#cat $var1*.converted > $var1.bowtie_v0_m1_all.bed.lifted.converted

cat $var1 > $var2.bowtie_v0_m1_$var3.bed