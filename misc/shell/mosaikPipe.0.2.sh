#!/bin/csh

#PBS -N mosaikPipe.0.2
#PBS -q flux
#PBS -l qos=lsa_flux
#PBS -A lsa_flux
#PBS -l nodes=10:ppn=2
#PBS -l walltime=5:00:00
#PBS -V
#PBS -m abe
#PBS -M $USER@umich.edu

# qsub -v var1=../../../Dmel/dm3_ambig.fa,var2=../../tiled/constExons_single_bp50_error0_tiled.fa mosaikPipe.0.2.sh

perl ../perlScripts/mosaikPipe.0.2.pl $var1 $var2
