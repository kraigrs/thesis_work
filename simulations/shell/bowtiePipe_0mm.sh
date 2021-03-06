#!/bin/sh

#PBS -S /bin/sh
#PBS -N bowtiePipe_0mm
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=2,walltime=4:00:00
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

#Usage: qsub -v var1=<file1>,var2=<file2> bowtiePipe.0.2.sh

perl ../perl/bowtiePipe_0mm.pl $var1 $var2

#var1=../../../Dmel/dm3_ref.fa var2=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.fa qsub bowtiePipe_0mm.sh
#var1=../../../Dmel/dm3_alt_line_40.fa var2=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.fa qsub bowtiePipe.0.2.sh