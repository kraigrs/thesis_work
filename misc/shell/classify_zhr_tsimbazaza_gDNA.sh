#!/bin/sh

#zhr

var1=../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=tsimbazaza var4=../mel_sim_data/Resequencing/zhr/final_info_rmDup.txt var5=../mel_sim_data/Resequencing/zhr var6=zhr_gDNA qsub classify.sh

#tsimbazaza

var1=../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=tsimbazaza var4=../mel_sim_data/Resequencing/tsimbazaza/final_info_rmDup.txt var5=../mel_sim_data/Resequencing/tsimbazaza var6=tsimbazaza_gDNA qsub classify.sh
