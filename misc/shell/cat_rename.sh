#!/bin/sh

#PBS -S /bin/sh
#PBS -N cat_rename
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -l nodes=1:ppn=1,walltime=0:30:00
#PBS -q flux
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

cat $var1\/*R1*.fastq > $var2
cat $var1\/*R2*.fastq > $var3