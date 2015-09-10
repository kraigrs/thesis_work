#!/bin/sh

#PBS -S /bin/sh
#PBS -N bowtie_trim_single
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=1,pmem=4000mb,walltime=10:00:00
#PBS -j oe
#PBS -V

####  End PBS preamble
#  Include the next three lines always

if [ "x${PBS_NODEFILE}" != "x" ] ; then
   cat $PBS_NODEFILE   # contains a list of the CPUs you were using if run with PBS
fi

cd $PBS_O_WORKDIR

#Usage: var1=<file1> var2=<file2> qsub bowtie_trim_single.sh

perl ../perl/bowtie_trim_single.pl $var1 $var2 $var3