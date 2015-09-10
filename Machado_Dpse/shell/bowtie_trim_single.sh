#!/bin/sh

#PBS -S /bin/sh
#PBS -N bowtie_trim_single
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

#Usage: var1=<file1> var2=<file2> qsub bowtie_trim_single.sh

perl ../perl/bowtie_trim_single.pl $var1 $var2