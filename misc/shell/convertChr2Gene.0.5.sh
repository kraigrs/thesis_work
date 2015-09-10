#!/bin/sh

#PBS -S /bin/sh
#PBS -N convertChr2Gene.0.5
#PBS -l nodes=1:ppn=1,walltime=1:00:00
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -M $USER@umich.edu
#PBS -m abe
#PBS -j oe
#PBS -V
echo "I ran on:"
cat $PBS_NODEFILE
cd $PBS_O_WORKDIR

#Usage: var1=<.bed file> var2=<gap file> var3=<constExons> qsub convertChr2Gene.0.5.sh

perl ../perlScripts/convertChr2Gene.0.5.pl $var1 $var2 $var3

#perl ../perlScripts/convertChr2Gene.0.5.pl ../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.zhr.mosaik.bed.lifted ../sim_sec_data/Sim_sec_gaps_gt19.bed ../McManus/constitutive_regions_overlap_filtered.bed
