#!/bin/sh

#PBS -S /bin/sh
#PBS -N samtools_cat
#PBS -A wittkopp_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=1,pmem=4000mb,walltime=24:00:00
#PBS -j oe
#PBS -V

####  End PBS preamble
#  Include the next three lines always

if [ "x${PBS_NODEFILE}" != "x" ] ; then
   cat $PBS_NODEFILE   # contains a list of the CPUs you were using if run with PBS
fi

cd $PBS_O_WORKDIR

#echo $var1

samtools cat -o $var1 $var2
#echo samtools cat -o $var1 $var2
