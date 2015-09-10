#!/bin/sh

#PBS -S /bin/sh
#PBS -N build_ebwt
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=1,walltime=1:00:00
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

bowtie-build ../../genomes/dm3/dm3.fa ../../genomes/dm3/dm3

#bowtie-build ../../genomes/zhr/zhr_AMBIG_all.fa ../../genomes/zhr/zhr_AMBIG_all

#bowtie-build ../../genomes/z30/z30_AMBIG_all.fa ../../genomes/z30/z30_AMBIG_all

#bowtie-build ../../genomes/tsimbazaza/droSim1_snp_indel.fa ../../genomes/tsimbazaza/droSim1_snp_indel
#bowtie-build ../../genomes/tsimbazaza/sim_extra.fa ../../genomes/tsimbazaza/sim_extra

#bowtie-build ../../genomes/droSec1/droSec1_snp_indel.fa ../../genomes/droSec1/droSec1_snp_indel
#bowtie-build ../../genomes/droSec1/sec_extra.fa ../../genomes/droSec1/sec_extra