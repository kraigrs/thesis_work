#!/bin/sh

#PBS -S /bin/sh
#PBS -N convertChr2Gene
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=1,walltime=0:30:00
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR
echo FILE = $var1

perl ../perl/convertChr2Gene.pl $var1 $var2 $var3
