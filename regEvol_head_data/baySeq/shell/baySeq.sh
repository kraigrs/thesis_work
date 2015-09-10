#!/bin/sh

#PBS -S /bin/sh
#PBS -N baySeq
#PBS -A wittkopp_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l nodes=1:ppn=12,pmem=4000mb,walltime=72:00:00
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

R CMD BATCH --vanilla $var1\.R $var1\.out