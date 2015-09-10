#!/bin/sh

#zhr
var1=../mel_sim_data/Resequencing/zhr/zhr_gDNA.mate1.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/Resequencing/zhr/zhr_gDNA.mate2.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/Resequencing/zhr/zhr_gDNA.mate1.tsimbazaza.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/Resequencing/zhr/zhr_gDNA.mate2.tsimbazaza.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh

#tsimbazaza
var1=../mel_sim_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate1.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate2.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate1.tsimbazaza.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate2.tsimbazaza.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
