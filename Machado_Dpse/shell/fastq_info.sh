#!/bin/sh

#PBS -S /bin/sh
#PBS -N fastq_info
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=5:00:00
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

fastx_quality_stats -i $var1.fastq -o $var1\_stats.txt
fastq_quality_boxplot_graph.sh -i $var1\_stats.txt -o $var1\_quality.png -t $var1
fastx_nucleotide_distribution_graph.sh -i $var1\_stats.txt -o $var1\_nuc.png -t $var1
