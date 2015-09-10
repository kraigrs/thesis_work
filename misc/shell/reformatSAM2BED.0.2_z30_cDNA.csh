#!/bin/csh

#PBS -N reformatSAM2BED.0.2_z30_cDNA
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

# z30

# mate 1
perl ../perlScripts/reformatSAM2BED.0.2.pl ../mel_mel_data/mRNA-Seq/z30_cDNA/z30.mate1.zhr_sanger_all.bowtie ../mel_mel_data/Resequencing/resequencing-assembly/zhr/zhr_sanger_rsq.chain yes
perl ../perlScripts/reformatSAM2BED.0.2.pl ../mel_mel_data/mRNA-Seq/z30_cDNA/z30.mate1.zhr_sanger_all.mosaik ../mel_mel_data/Resequencing/resequencing-assembly/zhr/zhr_sanger_rsq.chain yes
perl ../perlScripts/reformatSAM2BED.0.2.pl ../mel_mel_data/mRNA-Seq/z30_cDNA/z30.mate1.z30_sanger_all.bowtie ../mel_mel_data/Resequencing/resequencing-assembly/z30/z30_sanger_rsq.chain yes
perl ../perlScripts/reformatSAM2BED.0.2.pl ../mel_mel_data/mRNA-Seq/z30_cDNA/z30.mate1.z30_sanger_all.mosaik ../mel_mel_data/Resequencing/resequencing-assembly/z30/z30_sanger_rsq.chain yes

# mate 2
perl ../perlScripts/reformatSAM2BED.0.2.pl ../mel_mel_data/mRNA-Seq/z30_cDNA/z30.mate2.zhr_sanger_all.bowtie ../mel_mel_data/Resequencing/resequencing-assembly/zhr/zhr_sanger_rsq.chain yes
perl ../perlScripts/reformatSAM2BED.0.2.pl ../mel_mel_data/mRNA-Seq/z30_cDNA/z30.mate2.zhr_sanger_all.mosaik ../mel_mel_data/Resequencing/resequencing-assembly/zhr/zhr_sanger_rsq.chain yes
perl ../perlScripts/reformatSAM2BED.0.2.pl ../mel_mel_data/mRNA-Seq/z30_cDNA/z30.mate2.z30_sanger_all.bowtie ../mel_mel_data/Resequencing/resequencing-assembly/z30/z30_sanger_rsq.chain yes
perl ../perlScripts/reformatSAM2BED.0.2.pl ../mel_mel_data/mRNA-Seq/z30_cDNA/z30.mate2.z30_sanger_all.mosaik ../mel_mel_data/Resequencing/resequencing-assembly/z30/z30_sanger_rsq.chain yes
