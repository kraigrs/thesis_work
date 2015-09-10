#!/bin/sh

#PBS -S /bin/sh
#PBS -N samtools_sort
#PBS -A wittkopp_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=1,pmem=4000mb,walltime=10:00:00
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

tmp=${var1%*.bam}
samtools sort $tmp.bam $tmp.sorted
