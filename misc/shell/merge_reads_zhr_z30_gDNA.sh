#!/bin/sh

#zhr
var1=../mel_mel_data/Resequencing/zhr/zhr.mate1.zhr.mosaik.bed.lifted.converted var2=../mel_mel_data/Resequencing/zhr/zhr.mate1.z30.mosaik.bed.lifted.converted var3=../mel_mel_data/Resequencing/zhr/zhr.mate2.zhr.mosaik.bed.lifted.converted var4=../mel_mel_data/Resequencing/zhr/zhr.mate2.z30.mosaik.bed.lifted.converted qsub merge_reads.sh

#z30
var1=../mel_mel_data/Resequencing/z30/z30.mate1.zhr.mosaik.bed.lifted.converted var2=../mel_mel_data/Resequencing/z30/z30.mate1.z30.mosaik.bed.lifted.converted var3=../mel_mel_data/Resequencing/z30/z30.mate2.zhr.mosaik.bed.lifted.converted var4=../mel_mel_data/Resequencing/z30/z30.mate2.z30.mosaik.bed.lifted.converted qsub merge_reads.sh