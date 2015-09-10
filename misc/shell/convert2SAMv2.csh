#!/bin/csh

#PBS -N convert2SAMv2
#PBS -q flux
#PBS -l qos=lsa_flux
#PBS -A lsa_flux
#PBS -l nodes=1:ppn=1
#PBS -l walltime=120:00:00
#PBS -o $HOME/OUT.log
#PBS -e $HOME/ERR.err
#PBS -V
##PBS -m abe
##PBS -M $USER@umich.edu
setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

#Usage: qsub -v var1=<dir>,var2=<prefix>,var3=<ref1>,var4=<ref2> convert2SAMv2.csh

perl ../perlScripts/convert2SAMv2.pl $var1 $var2 $var3 $var4

