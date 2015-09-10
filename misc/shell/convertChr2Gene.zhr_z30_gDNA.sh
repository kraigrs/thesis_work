#!/bin/sh

#zhr
var1=../mel_mel_data/Resequencing/zhr/zhr.mate1.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/Resequencing/zhr/zhr.mate2.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/Resequencing/zhr/zhr.mate1.z30.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/Resequencing/zhr/zhr.mate2.z30.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh

#z30
var1=../mel_mel_data/Resequencing/z30/z30.mate1.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/Resequencing/z30/z30.mate2.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/Resequencing/z30/z30.mate1.z30.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/Resequencing/z30/z30.mate2.z30.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
