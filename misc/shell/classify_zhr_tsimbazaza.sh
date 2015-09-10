#!/bin/sh

#zhr

var1=../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=tsimbazaza var4=../mel_sim_data/mRNA-Seq/zhr/final_info_rmDup.txt var5=../mel_sim_data/mRNA-Seq/zhr var6=zhr_mRNA qsub classify.sh

#tsimbazaza

var1=../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=tsimbazaza var4=../mel_sim_data/mRNA-Seq/tsimbazaza/final_info_rmDup.txt var5=../mel_sim_data/mRNA-Seq/tsimbazaza var6=tsimbazaza_mRNA qsub classify.sh

#zhr+tsimbazaza

var1=../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=tsimbazaza var4=../mel_sim_data/mRNA-Seq/zhr+tsimbazaza/final_info_rmDup.txt var5=../mel_sim_data/mRNA-Seq/zhr+tsimbazaza var6=zhr+tsimbazaza_mRNA qsub classify.sh

#zhrXtsimbazaza

var1=../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=tsimbazaza var4=../mel_sim_data/mRNA-Seq/zhrXtsimbazaza/final_info_rmDup.txt var5=../mel_sim_data/mRNA-Seq/zhrXtsimbazaza var6=zhrXtsimbazaza_mRNA qsub classify.sh

#tsimbazazaXzhr

var1=../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=tsimbazaza var4=../mel_sim_data/mRNA-Seq/tsimbazazaXzhr/final_info_rmDup.txt var5=../mel_sim_data/mRNA-Seq/tsimbazazaXzhr var6=tsimbazazaXzhr_mRNA qsub classify.sh
