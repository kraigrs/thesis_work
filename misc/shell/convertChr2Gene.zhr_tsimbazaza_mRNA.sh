#!/bin/sh

#zhr
var1=../mel_sim_data/mRNA-Seq/zhr/zhr.mate1.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/zhr/zhr.mate2.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/zhr/zhr.mate1.tsimbazaza.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/zhr/zhr.mate2.tsimbazaza.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh

#tsimbazaza
var1=../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate1.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate2.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate1.tsimbazaza.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate2.tsimbazaza.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh

#zhr+tsimbazaza
var1=../mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza.mate1.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza.mate2.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza.mate1.tsimbazaza.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza.mate2.tsimbazaza.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh

#zhrXtsimbazaza
var1=../mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza.mate1.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza.mate2.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza.mate1.tsimbazaza.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza.mate2.tsimbazaza.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh

#tsimbazazaXzhr
var1=../mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr.mate1.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr.mate2.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr.mate1.tsimbazaza.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr.mate2.tsimbazaza.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
