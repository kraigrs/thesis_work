#!/bin/sh

#zhr

var1=../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=z30 var4=../mel_mel_data/Resequencing/zhr/final_info_rmDup.txt var5=../mel_mel_data/Resequencing/zhr var6=zhr_gDNA qsub classify.sh

#z30

var1=../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=z30 var4=../mel_mel_data/Resequencing/z30/final_info_rmDup.txt var5=../mel_mel_data/Resequencing/z30 var6=z30_gDNA qsub classify.sh
