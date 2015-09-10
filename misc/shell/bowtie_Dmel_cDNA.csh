#!/bin/csh

#PBS -N bowtie_Dmel_cDNA
#PBS -q flux
#PBS -l qos=lsa_flux
#PBS -A lsa_flux
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -o $HOME/OUT.log
#PBS -e $HOME/ERR.err
#PBS -V
##PBS -m abe
##PBS -M $USER@umich.edu
setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

perl ../perlScripts/bowtiePipe.pl ../Dmel/dm3.fa ../Dsec/droSec1.fa ../mel_sec_data/Dmel_r1.fa ../mel_sec_data/Dmel_r2.fa
