#!/bin/sh

#zhr
var1=../mel_sim_data/Resequencing/zhr/zhr_gDNA.mate1.zhr.mosaik.bed.lifted.converted var2=../mel_sim_data/Resequencing/zhr/zhr_gDNA.mate1.tsimbazaza.mosaik.bed.lifted.converted var3=../mel_sim_data/Resequencing/zhr/zhr_gDNA.mate2.zhr.mosaik.bed.lifted.converted var4=../mel_sim_data/Resequencing/zhr/zhr_gDNA.mate2.tsimbazaza.mosaik.bed.lifted.converted qsub merge_reads.sh

#tsimbazaza
var1=../mel_sim_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate1.zhr.mosaik.bed.lifted.converted var2=../mel_sim_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate1.tsimbazaza.mosaik.bed.lifted.converted var3=../mel_sim_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate2.zhr.mosaik.bed.lifted.converted var4=../mel_sim_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate2.tsimbazaza.mosaik.bed.lifted.converted qsub merge_reads.sh
