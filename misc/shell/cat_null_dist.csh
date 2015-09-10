#!/bin/csh

#PBS -N imprinting.csh
#PBS -q flux
#PBS -l qos=lsa_flux
#PBS -A lsa_flux
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -o $HOME/OUT.log
#PBS -e $HOME/ERR.err
#PBS -V
#PBS -m abe
#PBS -M $USER@umich.edu
setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

cat null_dist*/chr2L/zhr_z30.chr2L_wind500000_step10000.imprinting.txt > null/chr2L/zhr_z30.chr2L_wind500000_step10000.imprinting.txt
cat null_dist*/chr2R/zhr_z30.chr2R_wind500000_step10000.imprinting.txt > null/chr2R/zhr_z30.chr2R_wind500000_step10000.imprinting.txt
cat null_dist*/chr3L/zhr_z30.chr3L_wind500000_step10000.imprinting.txt > null/chr3L/zhr_z30.chr3L_wind500000_step10000.imprinting.txt
cat null_dist*/chr3R/zhr_z30.chr3R_wind500000_step10000.imprinting.txt > null/chr3R/zhr_z30.chr3R_wind500000_step10000.imprinting.txt
cat null_dist*/chrX/zhr_z30.chrX_wind500000_step10000.imprinting.txt > null/chrX/zhr_z30.chrX_wind500000_step10000.imprinting.txt
cat null_dist*/chr4/zhr_z30.chr4_wind500000_step10000.imprinting.txt > null/chr4/zhr_z30.chr4_wind500000_step10000.imprinting.txt

cat null_dist*/chr2L/zhr_z30.chr2L_k50.imprinting.txt > null/chr2L/zhr_z30.chr2L_k50.imprinting.txt
cat null_dist*/chr2R/zhr_z30.chr2R_k50.imprinting.txt > null/chr2R/zhr_z30.chr2R_k50.imprinting.txt
cat null_dist*/chr3L/zhr_z30.chr3L_k50.imprinting.txt > null/chr3L/zhr_z30.chr3L_k50.imprinting.txt
cat null_dist*/chr3R/zhr_z30.chr3R_k50.imprinting.txt > null/chr3R/zhr_z30.chr3R_k50.imprinting.txt
cat null_dist*/chrX/zhr_z30.chrX_k50.imprinting.txt > null/chrX/zhr_z30.chrX_k50.imprinting.txt
cat null_dist*/chr4/zhr_z30.chr4_k50.imprinting.txt > null/chr4/zhr_z30.chr4_k50.imprinting.txt
