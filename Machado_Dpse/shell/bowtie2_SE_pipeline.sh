#!/bin/sh

#PBS -S /bin/sh
#PBS -N bowtie2_SE_pipeline
#PBS -A wittkopp_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l nodes=1:ppn=12,pmem=4000mb,walltime=10:00:00
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M kraigrs@umich.edu

####  End PBS preamble
#  Include the next three lines always

if [ "x${PBS_NODEFILE}" != "x" ] ; then
   cat $PBS_NODEFILE   # contains a list of the CPUs you were using if run with PBS
fi

cd $PBS_O_WORKDIR

#Usage: var1=<file1> var2=<file2> qsub bowtie2_SE_pipeline.sh

perl ../perl/bowtie2_SE_pipeline.pl $var1 $var2