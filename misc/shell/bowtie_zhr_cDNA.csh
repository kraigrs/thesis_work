#!/bin/csh

#PBS -N bowtie_zhr_cDNA
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

perl ../perlScripts/bowtiePipe.pl ../mel_mel_data/Resequencing/resequencing-assembly/zhr/zhr_sanger_all.fa ../mel_sim_data/Resequencing/resequencing-assembly/droSim1_snp_indel_GATC.fa ../mel_sim_data/mRNA-Seq/zhr_cDNA/zhr.mate1.fastq ../mel_sim_data/mRNA-Seq/zhr_cDNA/zhr.mate2.fastq

