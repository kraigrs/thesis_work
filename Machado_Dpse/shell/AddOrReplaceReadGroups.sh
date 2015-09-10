#!/bin/sh

#PBS -S /bin/sh
#PBS -N AddOrReplaceReadGroups
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=5:00:00
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

#java -Xmx4g -jar $Picard_JAVA/AddOrReplaceReadGroups.jar INPUT=../../data/Dpse_TL.dpse-all-gene-r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.bam OUTPUT=../../data/Dpse_TL.dpse-all-gene-r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.group.bam RGLB=single RGPL=Illumina RGPU=HiSeq RGSM=Dpse_TL

#java -Xmx4g -jar $Picard_JAVA/AddOrReplaceReadGroups.jar INPUT=../../data/Dbog_Toro1.dpse-all-gene-r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.bam OUTPUT=../../data/Dbog_Toro1.dpse-all-gene-r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.group.bam RGLB=paired RGPL=Illumina RGPU=HiSeq RGSM=Dbog_Toro1

java -Xmx4g -jar $Picard_JAVA/AddOrReplaceReadGroups.jar INPUT=../../data/Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.bam OUTPUT=../../data/Dpse_TL.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.group.bam RGLB=single RGPL=Illumina RGPU=HiSeq RGSM=Dpse_TL

java -Xmx4g -jar $Picard_JAVA/AddOrReplaceReadGroups.jar INPUT=../../data/Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.bam OUTPUT=../../data/Dbog_Toro1.dpse_r3.1.bowtie2_very_sensitive.merged.sorted.rmdup.group.bam RGLB=paired RGPL=Illumina RGPU=HiSeq RGSM=Dbog_Toro1