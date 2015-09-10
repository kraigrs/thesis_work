#!/bin/sh

#PBS -S /bin/sh
#PBS -N multiple_aln_FSA
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -l procs=1,mem=3gb,walltime=20:00:00
#PBS -q flux
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $PBS_O_HOST
cd $PBS_O_WORKDIR

#Usage: var1=<file1>... qsub multiple_aln_FSA.sh

#var1=dm3,var2=zhr,var3=z30,var4=sim,var5=sec,var6=constitutive_regions_overlap_filtered2.dm3.fa,var7=constitutive_regions_overlap_filtered2.zhr.fa,var8=constitutive_regions_overlap_filtered2.z30.fa,var9=constitutive_regions_overlap_filtered2.sim_snp_indel.fa,varA=constitutive_regions_overlap_filtered2.sec_snp_indel.fa multiple_aln_FSA.sh

#perl multiple_aln_FSA.pl $var1 $var2 $var3 $var4 $var5 $var6 $var7 $var8 $var9 $varA

perl multiple_aln_FSA.pl dm3 zhr z30 sim sec constitutive_regions_overlap_filtered2.dm3.fa constitutive_regions_overlap_filtered2.zhr.fa constitutive_regions_overlap_filtered2.z30.fa constitutive_regions_overlap_filtered2.sim_snp_indel.fa constitutive_regions_overlap_filtered2.sec_snp_indel.fa
