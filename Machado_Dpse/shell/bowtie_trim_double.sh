#!/bin/sh

#PBS -S /bin/sh
#PBS -N bowtie_trim_double
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

#Usage: var1=<ref1> var2=<ref2> var3=<readfile> qsub bowtie_trim_double.sh

perl ../perl/bowtie_trim_double.pl $var1 $var2 $var3