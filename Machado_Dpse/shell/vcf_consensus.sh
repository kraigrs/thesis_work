#!/bin/sh

#PBS -S /bin/sh
#PBS -N vcf_consensus
#PBS -A wittkopp_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=1:00:00
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

bgzip -c $var2 > $var2.gz
tabix -p vcf $var2\.gz
cat $var1.fasta | vcf-consensus $var2\.gz > $var1.$var3.fasta
