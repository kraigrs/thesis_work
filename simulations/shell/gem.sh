#!/bin/sh

#PBS -S /bin/sh
#PBS -N GEM_mappability
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=1,walltime=10:00:00
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

#gem-do-index -i $var1 -o $var2
gem-mappability -I $var1 -l 36 -m $var2 -o $var3