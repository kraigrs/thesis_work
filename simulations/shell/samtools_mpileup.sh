#!/bin/sh

#PBS -S /bin/sh
#PBS -N samtools_mpileup
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -l nodes=1:ppn=1,walltime=3:00:00
#PBS -q flux
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -M $USER@umich.edu
#setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

#Usage: qsub -v var1=<sam1>,var2=<sam2>,var3=<ref>,var4=<bam2> samtools_mpileup.sh

#samtools view -S -b -T $var3 $var1 > temp1.bam
#samtools view -S -b -T $var3 $var2 > temp2.bam

#samtools sort temp1.bam temp1.bam.sorted
#samtools sort temp2.bam temp2.bam.sorted

#samtools index temp1.bam.sorted
#samtools index temp1.bam.sorted

samtools mpileup -f $var1 $var2 $var3 $var4 > $var5

#qsub -v var1=../../../Graze/berlin-updated-exonic-regions.fasta,var2=../../../Graze/SRR389095.mate1_2.berlin-updated-exonic-regions.bowtie_v1_m1.sorted.bam,var3=../../../Graze/SRR389095.mate1_7.berlin-updated-exonic-regions.bowtie_v1_m1.sorted.bam,var4=../../../Graze/SRR389095.mate1_8.berlin-updated-exonic-regions.bowtie_v1_m1.sorted.bam,var5=../../../Graze/SRR389095.berlin-updated-exonic-regions.bowtie_v1_m1.pileup.txt samtools_mpileup.sh

#qsub -v var1=../../../Graze/berlin-updated-exonic-regions.fasta,var2=../../../Graze/SRR389095.mate1_2.berlin-updated-exonic-regions.bowtie_v2_m1.sorted.bam,var3=../../../Graze/SRR389095.mate1_7.berlin-updated-exonic-regions.bowtie_v2_m1.sorted.bam,var4=../../../Graze/SRR389095.mate1_8.berlin-updated-exonic-regions.bowtie_v2_m1.sorted.bam,var5=../../../Graze/SRR389095.berlin-updated-exonic-regions.bowtie_v2_m1.pileup.txt samtools_mpileup.sh

#qsub -v var1=../../../Graze/berlin-updated-exonic-regions.fasta,var2=../../../Graze/SRR389095.mate1_2.berlin-updated-exonic-regions.bowtie_v3_m1.sorted.bam,var3=../../../Graze/SRR389095.mate1_7.berlin-updated-exonic-regions.bowtie_v3_m1.sorted.bam,var4=../../../Graze/SRR389095.mate1_8.berlin-updated-exonic-regions.bowtie_v3_m1.sorted.bam,var5=../../../Graze/SRR389095.berlin-updated-exonic-regions.bowtie_v3_m1.pileup.txt samtools_mpileup.sh


#qsub -v var1=../../../Graze/c1674-updated-exonic-regions.fasta,var2=../../../Graze/SRR389095.mate1_2.c1674-updated-exonic-regions.bowtie_v1_m1.sorted.bam,var3=../../../Graze/SRR389095.mate1_7.c1674-updated-exonic-regions.bowtie_v1_m1.sorted.bam,var4=../../../Graze/SRR389095.mate1_8.c1674-updated-exonic-regions.bowtie_v1_m1.sorted.bam,var5=../../../Graze/SRR389095.c1674-updated-exonic-regions.bowtie_v1_m1.pileup.txt samtools_mpileup.sh

#qsub -v var1=../../../Graze/c1674-updated-exonic-regions.fasta,var2=../../../Graze/SRR389095.mate1_2.c1674-updated-exonic-regions.bowtie_v2_m1.sorted.bam,var3=../../../Graze/SRR389095.mate1_7.c1674-updated-exonic-regions.bowtie_v2_m1.sorted.bam,var4=../../../Graze/SRR389095.mate1_8.c1674-updated-exonic-regions.bowtie_v2_m1.sorted.bam,var5=../../../Graze/SRR389095.c1674-updated-exonic-regions.bowtie_v2_m1.pileup.txt samtools_mpileup.sh

#qsub -v var1=../../../Graze/c1674-updated-exonic-regions.fasta,var2=../../../Graze/SRR389095.mate1_2.c1674-updated-exonic-regions.bowtie_v3_m1.sorted.bam,var3=../../../Graze/SRR389095.mate1_7.c1674-updated-exonic-regions.bowtie_v3_m1.sorted.bam,var4=../../../Graze/SRR389095.mate1_8.c1674-updated-exonic-regions.bowtie_v3_m1.sorted.bam,var5=../../../Graze/SRR389095.c1674-updated-exonic-regions.bowtie_v3_m1.pileup.txt samtools_mpileup.sh



#qsub -v var1=../../../Graze/berlin-updated-exonic-regions_fsa_masked.fasta,var2=../../../Graze/SRR389095.mate1_2.berlin-updated-exonic-regions_fsa_masked.bowtie_v0_m1.sorted.bam,var3=../../../Graze/SRR389095.mate1_7.berlin-updated-exonic-regions_fsa_masked.bowtie_v0_m1.sorted.bam,var4=../../../Graze/SRR389095.mate1_8.berlin-updated-exonic-regions_fsa_masked.bowtie_v0_m1.sorted.bam,var5=../../../Graze/SRR389095.mate1.berlin-updated-exonic-regions_fsa_masked.bowtie_v0_m1.pileup.txt samtools_mpileup.sh

#qsub -v var1=../../../Graze/c1674-updated-exonic-regions_fsa_masked.fasta,var2=../../../Graze/SRR389095.mate1_2.c1674-updated-exonic-regions_fsa_masked.bowtie_v0_m1.sorted.bam,var3=../../../Graze/SRR389095.mate1_7.c1674-updated-exonic-regions_fsa_masked.bowtie_v0_m1.sorted.bam,var4=../../../Graze/SRR389095.mate1_8.c1674-updated-exonic-regions_fsa_masked.bowtie_v0_m1.sorted.bam,var5=../../../Graze/SRR389095.mate1.c1674-updated-exonic-regions_fsa_masked.bowtie_v0_m1.pileup.txt samtools_mpileup.sh