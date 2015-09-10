#!/bin/sh

#PBS -S /bin/sh
#PBS -N pileup_SNP_ASE
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=1,walltime=2:00:00
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

perl ../perl/pileup_SNP_ASE.pl $var1 $var2 $var3 $var4