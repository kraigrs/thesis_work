#!/bin/sh

#PBS -S /bin/sh
#PBS -N samToBed
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

#Usage: var1=<file1> qsub samToBed.sh

perl ../perl/samToBed.pl $var1