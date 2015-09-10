#!/bin/sh

#PBS -S /bin/sh
#PBS -N liftOver
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

liftOver $var1 $var2 $var1.lifted $var1.un

#tmp=${var1%*.bed}
#liftOver $var1.bed $var2 $var1.$var3.bed $var1.$var3.un -minMatch=0.5
