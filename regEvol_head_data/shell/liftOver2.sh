#!/bin/sh

#PBS -S /bin/sh
#PBS -N liftOver2
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

liftOver $var1 $var2 $var1.temp.lifted $var1.un
liftOver $var1.temp.lifted $var3 $var1.lifted $var1.temp.un

#tmp=${var1%*.bed}

#liftOver $tmp.bed $var2 $tmp.bed.temp $tmp.un.temp -minMatch=0.5
#liftOver $tmp.bed.temp $var3 $tmp.$var5.bed $tmp.$var5.un -minMatch=0.5

#liftOver $tmp.un.temp $var4 $tmp.$var5\_extra.bed $tmp.$var5\_extra.un -minMatch=0.5
