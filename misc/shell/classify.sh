#!/bin/sh

#PBS -S /bin/sh
#PBS -N classify.0.7
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -l nodes=1:ppn=1,walltime=1:00:00
#PBS -q flux
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo Host is $PBS_O_HOST
cd $PBS_O_WORKDIR

perl ../perlScripts/classify.0.7.pl $var1 $var2 $var3 $var4 $var5 $var6
