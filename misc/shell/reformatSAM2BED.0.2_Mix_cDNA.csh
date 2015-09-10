#!/bin/csh

#PBS -N reformatSAM2BED.0.2_Mix_cDNA
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

# Dmel+Dsec

# mate 1
perl ../perlScripts/reformatSAM2BED.0.2.pl ../mel_sec_data/Mix_r1.dm3.bowtie no NA
perl ../perlScripts/reformatSAM2BED.0.2.pl ../mel_sec_data/Mix_r1.droSec1.bowtie yes ../Dsec/droSec1ToDm3.over.chain 

# mate 2
perl ../perlScripts/reformatSAM2BED.0.2.pl ../mel_sec_data/Mix_r2.dm3.bowtie no NA
perl ../perlScripts/reformatSAM2BED.0.2.pl ../mel_sec_data/Mix_r2.droSec1.bowtie yes ../Dsec/droSec1ToDm3.over.chain


