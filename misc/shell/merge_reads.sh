#!/bin/sh

#PBS -S /bin/sh
#PBS -N merge_reads
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -l nodes=1:ppn=1,walltime=3:00:00
#PBS -q flux
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

sh Characterize_Each_Read.sh $var1 $var2 $var3 $var4