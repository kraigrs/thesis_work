#!/bin/sh

#PBS -S /bin/sh
#PBS -N pairwise_aln_FSA
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -l procs=1,mem=3gb,walltime=8:00:00
#PBS -q flux
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $PBS_O_HOST
cd $PBS_O_WORKDIR

#Usage: var1=<file1>... multiple_aln_FSA.sh

perl pairwise_aln_FSA.pl $var1 $var2 $var3 $var4

#var1=zhr var2=z30 var3=constitutive_regions_overlap_filtered2.zhr.fa var4=constitutive_regions_overlap_filtered2.z30.fa qsub pairwise_aln_FSA.sh

#var1=sim var2=sec var3=constitutive_regions_overlap_filtered2.sim_snp_indel_extra.fa var4=constitutive_regions_overlap_filtered2.sec_snp_indel_extra.fa qsub pairwise_aln_FSA.sh

#var1=zhr var2=sim var3=constitutive_regions_overlap_filtered2.zhr.fa var4=constitutive_regions_overlap_filtered2.sim_snp_indel_extra.fa qsub pairwise_aln_FSA.sh