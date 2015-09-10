#!/bin/sh

#PBS -S /bin/sh
#PBS -N bowtie_trim_single_double
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

#Usage: var1=<ref1> var2=<ref2> var3=<readfile> qsub bowtie_trim_double.sh

perl ../perl/bowtie_trim_single_double.pl $var1 $var2 $var3 $var4