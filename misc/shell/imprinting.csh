#!/bin/csh

#PBS -N imprinting.csh
#PBS -q flux
#PBS -l qos=lsa_flux
#PBS -A lsa_flux
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -o $HOME/OUT.log
#PBS -e $HOME/ERR.err
#PBS -V
#PBS -m abe
#PBS -M $USER@umich.edu
setenv NPROCS `wc -l < $PBS_NODEFILE`
echo This job is allocated $NPROCS cpus
echo host is $HOST
cd $PBS_O_WORKDIR

# qsub -v var1=1250,var2=134,var3=../McManus/constitutive_genes.txt,var4=../mel_mel_data/zhr_z30_MB_cut20_genes.txt,var5=../mel_mel_data/null_dist/all_genes_intersect.bed,var6=../mel_mel_data/null_dist,var7=500000,var8=10000,var9=50 imprinting.csh

# qsub -v var1=1250,var2=134,var3=../McManus/constitutive_genes.txt,var4=../mel_mel_data/zhr_z30_MB_cut20_genes.txt,var5=../mel_mel_data/null_dist2/all_genes_intersect.bed,var6=../mel_mel_data/null_dist2,var7=500000,var8=10000,var9=50 imprinting.csh

# qsub -v var1=1250,var2=134,var3=../McManus/constitutive_genes.txt,var4=../mel_mel_data/zhr_z30_MB_cut20_genes.txt,var5=../mel_mel_data/null_dist3/all_genes_intersect.bed,var6=../mel_mel_data/null_dist3,var7=500000,var8=10000,var9=50 imprinting.csh

# qsub -v var1=1250,var2=134,var3=../McManus/constitutive_genes.txt,var4=../mel_mel_data/zhr_z30_MB_cut20_genes.txt,var5=../mel_mel_data/null_dist4/all_genes_intersect.bed,var6=../mel_mel_data/null_dist4,var7=500000,var8=10000,var9=50 imprinting.csh

# qsub -v var1=1250,var2=134,var3=../McManus/constitutive_genes.txt,var4=../mel_mel_data/zhr_z30_MB_cut20_genes.txt,var5=../mel_mel_data/null_dist5/all_genes_intersect.bed,var6=../mel_mel_data/null_dist5,var7=500000,var8=10000,var9=50 imprinting.csh

# qsub -v var1=1250,var2=134,var3=../McManus/constitutive_genes.txt,var4=../mel_mel_data/zhr_z30_MB_cut20_genes.txt,var5=../mel_mel_data/null_dist6/all_genes_intersect.bed,var6=../mel_mel_data/null_dist6,var7=500000,var8=10000,var9=50 imprinting.csh

# qsub -v var1=1250,var2=134,var3=../McManus/constitutive_genes.txt,var4=../mel_mel_data/zhr_z30_MB_cut20_genes.txt,var5=../mel_mel_data/null_dist7/all_genes_intersect.bed,var6=../mel_mel_data/null_dist7,var7=500000,var8=10000,var9=50 imprinting.csh

# qsub -v var1=1250,var2=134,var3=../McManus/constitutive_genes.txt,var4=../mel_mel_data/zhr_z30_MB_cut20_genes.txt,var5=../mel_mel_data/null_dist8/all_genes_intersect.bed,var6=../mel_mel_data/null_dist8,var7=500000,var8=10000,var9=50 imprinting.csh

../perlScripts/imprinting_nullDist.pl $var1 $var2 $var3 $var4 $var5 $var6 $var7 $var8 $var9
