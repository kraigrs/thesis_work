#!/bin/sh

#PBS -S /bin/sh
#PBS -N liftOver
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
echo FILE = $var1

#liftOver $var1 $var2 $var1.lifted $var1.un

liftOver $var1.bed $var2 $var1.$var3.bed $var1.$var3.un -minMatch=0.5
