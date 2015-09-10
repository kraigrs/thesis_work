#!/bin/sh

#zhr

var1=../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=z30 var4=../mel_mel_data/mRNA-Seq/zhr/final_info_rmDup.txt var5=../mel_mel_data/mRNA-Seq/zhr var6=zhr_mRNA qsub classify.sh

#z30

var1=../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=z30 var4=../mel_mel_data/mRNA-Seq/z30/final_info_rmDup.txt var5=../mel_mel_data/mRNA-Seq/z30 var6=z30_mRNA qsub classify.sh

#zhrXz30

var1=../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=z30 var4=../mel_mel_data/mRNA-Seq/zhrXz30/final_info_rmDup.txt var5=../mel_mel_data/mRNA-Seq/zhrXz30 var6=zhrXz30_mRNA qsub classify.sh

#z30Xzhr

var1=../McManus/constitutive_regions_overlap_filtered.bed var2=zhr var3=z30 var4=../mel_mel_data/mRNA-Seq/z30Xzhr/final_info_rmDup.txt var5=../mel_mel_data/mRNA-Seq/z30Xzhr var6=z30Xzhr_mRNA qsub classify.sh
