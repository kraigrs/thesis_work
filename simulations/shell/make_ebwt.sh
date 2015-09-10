#!/bin/sh

#PBS -S /bin/sh
#PBS -N make_ebwt
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -l nodes=1:ppn=1,walltime=5:00:00
#PBS -q flux
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

#Usage: qsub -v var1=<fasta>,var2=<ebwt> make_ebwt.sh

bowtie-build $var1 $var2