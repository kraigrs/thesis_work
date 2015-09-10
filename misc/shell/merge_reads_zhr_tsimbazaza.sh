#!/bin/sh

#zhr
var1=../mel_sim_data/mRNA-Seq/zhr/zhr.mate1.zhr.mosaik.bed.lifted.converted var2=../mel_sim_data/mRNA-Seq/zhr/zhr.mate1.tsimbazaza.mosaik.bed.lifted.converted var3=../mel_sim_data/mRNA-Seq/zhr/zhr.mate2.zhr.mosaik.bed.lifted.converted var4=../mel_sim_data/mRNA-Seq/zhr/zhr.mate2.tsimbazaza.mosaik.bed.lifted.converted qsub merge_reads.sh

#tsimbazaza
var1=../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate1.zhr.mosaik.bed.lifted.converted var2=../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate1.tsimbazaza.mosaik.bed.lifted.converted var3=../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate2.zhr.mosaik.bed.lifted.converted var4=../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate2.tsimbazaza.mosaik.bed.lifted.converted qsub merge_reads.sh

#zhrXtsimbazaza
var1=../mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza.mate1.zhr.mosaik.bed.lifted.converted var2=../mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza.mate1.tsimbazaza.mosaik.bed.lifted.converted var3=../mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza.mate2.zhr.mosaik.bed.lifted.converted var4=../mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza.mate2.tsimbazaza.mosaik.bed.lifted.converted qsub merge_reads.sh

#tsimbazazaXzhr
var1=../mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr.mate1.zhr.mosaik.bed.lifted.converted var2=../mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr.mate1.tsimbazaza.mosaik.bed.lifted.converted var3=../mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr.mate2.zhr.mosaik.bed.lifted.converted var4=../mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr.mate2.tsimbazaza.mosaik.bed.lifted.converted qsub merge_reads.sh

#zhr+tsimbazaza
var1=../mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza.mate1.zhr.mosaik.bed.lifted.converted var2=../mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza.mate1.tsimbazaza.mosaik.bed.lifted.converted var3=../mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza.mate2.zhr.mosaik.bed.lifted.converted var4=../mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza.mate2.tsimbazaza.mosaik.bed.lifted.converted qsub merge_reads.sh