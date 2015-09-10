#!/bin/csh

#PBS -N bowtie_z30_cDNA
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

perl ../perlScripts/bowtiePipe.pl ../mel_mel_data/Resequencing/resequencing-assembly/zhr/zhr_sanger_all.fa ../mel_mel_data/Resequencing/resequencing-assembly/z30/z30_sanger_all.fa ../mel_mel_data/mRNA-Seq/z30_cDNA/z30.mate1.fastq ../mel_mel_data/mRNA-Seq/z30_cDNA/z30.mate2.fastq

