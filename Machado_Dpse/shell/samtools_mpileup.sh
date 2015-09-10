#!/bin/sh

#PBS -S /bin/sh
#PBS -N samtools_mpileup
#PBS -A wittkopp_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=15:00:00
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

samtools mpileup -f $var1.fasta $var3\/$var2.bam > $var3\/$var2.mpileup.txt
