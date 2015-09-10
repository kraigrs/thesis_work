#!/bin/sh

#PBS -S /bin/sh
#PBS -N samtoolsPipe
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=4,walltime=6:00:00
#PBS -l pmem=8000mb
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

perl ../perl/samtoolsPipe.0.2.pl $var1 $var2

#var1=../../../Dmel/dm3_ref.fa var2=../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_n3_e161_l36_m1.sam qsub samtoolsPipe.sh