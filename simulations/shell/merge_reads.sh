#!/bin/sh

#PBS -S /bin/sh
#PBS -N merge_reads
#PBS -A lsa_flux
#PBS -l qos=lsa_flux
#PBS -l nodes=1:ppn=1,walltime=10:00:00
#PBS -q flux
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

#Usage: qsub -v var1=<file1>,var2=<file2> merge_reads.sh

sh Characterize_Each_Read.sh $var1 $var2