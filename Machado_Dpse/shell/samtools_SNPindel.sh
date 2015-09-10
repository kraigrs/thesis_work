#!/bin/sh

#PBS -S /bin/sh
#PBS -N samtools_SNPindel
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=4:00:00
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

#Usage: var1=<ref>.[fa/fasta] var2=<file>.bam var3=<data directory> qsub samtools_SNPindel.sh

samtools mpileup -uf $var1.fasta $var3\/$var2.bam | bcftools view -bvcg - > $var3\/$var2.raw.bcf
bcftools view $var3\/$var2.raw.bcf | vcfutils.pl varFilter -d20 -Q30 > $var3\/$var2.flt.vcf
#vcfutils.pl vcf2fq $var3\/$var2.flt.vcf > $var2.SAMtools_alt.fq
