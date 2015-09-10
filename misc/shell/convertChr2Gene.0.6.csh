#!/bin/csh

#PBS -N convertChr2Gene.0.6 
#PBS -q flux
#PBS -l qos=flux
#PBS -A lsa_flux
#PBS -l nodes=1:ppn=1
#PBS -l walltime=2:00:00
#PBS -o $HOME/OUT.log
#PBS -e $HOME/ERR.err
#PBS -V
#PBS -m abe
#PBS -M $USER@umich.edu
setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

#Usage: qsub -v var1=<.bed file>,var2=<constExons> convertChr2Gene.0.6.csh

perl ../perlScripts/convertChr2Gene.0.6.pl $var1 $var2
