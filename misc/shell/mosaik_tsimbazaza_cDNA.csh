#!/bin/csh

#PBS -N mosaik_tsimbazaza_cDNA
#PBS -q flux
#PBS -l qos=lsa_flux
#PBS -A lsa_flux
#PBS -l nodes=1:ppn=8
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

perl ../perlScripts/mosaikPipe.pl ../mel_mel_data/Resequencing/resequencing-assembly/zhr/zhr_sanger_all.fa ../mel_sim_data/Resequencing/resequencing-assembly/droSim1_snp_indel_GATC.fa ../mel_sim_data/mRNA-Seq/tsimbazaza_cDNA/tsimbazaza.mate1.fastq ../mel_sim_data/mRNA-Seq/tsimbazaza_cDNA/tsimbazaza.mate2.fastq
