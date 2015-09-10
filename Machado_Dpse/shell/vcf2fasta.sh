#!/bin/sh

#PBS -S /bin/sh
#PBS -N vcf2fasta
#PBS -A wittkopp_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=10:00:00
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu

####  End PBS preamble
#  Include the next three lines always

if [ "x${PBS_NODEFILE}" != "x" ] ; then
   cat $PBS_NODEFILE   # contains a list of the CPUs you were using if run with PBS
fi

cd $PBS_O_WORKDIR

#Usage: var1=<pileup> var2=<VCF> var3=<output> qsub filter_VCF.sh

#var1=../../data/Dpse_TL.dpse-all-gene-r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.mpileup.txt var2=../../data/Dpse_TL.dpse-all-gene-r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.raw_variants.vcf var3=../../data/Dpse_TL.dpse-all-gene-r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.GATK_highQual.vcf qsub filter_VCF.sh

#var1=../../data/Dbog_Toro1.dpse-all-gene-r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.mpileup.txt var2=../../data/Dbog_Toro1.dpse-all-gene-r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.raw_variants.vcf var3=../../data/Dbog_Toro1.dpse-all-gene-r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.GATK_highQual.vcf qsub filter_VCF.sh

perl ../perl/vcf2fasta.pl $var1 $var2 $var3