#!/bin/sh

#zhr
var1=../mel_mel_data/mRNA-Seq/zhr/zhr.mate1.zhr.mosaik.bed.lifted.converted var2=../mel_mel_data/mRNA-Seq/zhr/zhr.mate1.z30.mosaik.bed.lifted.converted var3=../mel_mel_data/mRNA-Seq/zhr/zhr.mate2.zhr.mosaik.bed.lifted.converted var4=../mel_mel_data/mRNA-Seq/zhr/zhr.mate2.z30.mosaik.bed.lifted.converted qsub merge_reads.sh

#z30
var1=../mel_mel_data/mRNA-Seq/z30/z30.mate1.zhr.mosaik.bed.lifted.converted var2=../mel_mel_data/mRNA-Seq/z30/z30.mate1.z30.mosaik.bed.lifted.converted var3=../mel_mel_data/mRNA-Seq/z30/z30.mate2.zhr.mosaik.bed.lifted.converted var4=../mel_mel_data/mRNA-Seq/z30/z30.mate2.z30.mosaik.bed.lifted.converted qsub merge_reads.sh

#zhrXz30
var1=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.zhr.mosaik.bed.lifted.converted var2=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.z30.mosaik.bed.lifted.converted var3=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.zhr.mosaik.bed.lifted.converted var4=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.z30.mosaik.bed.lifted.converted qsub merge_reads.sh

#z30Xzhr
var1=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.zhr.mosaik.bed.lifted.converted var2=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.z30.mosaik.bed.lifted.converted var3=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.zhr.mosaik.bed.lifted.converted var4=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.z30.mosaik.bed.lifted.converted qsub merge_reads.sh

# gap filtered files

#zhr
#qsub -v var1=../mel_mel_data/mRNA-Seq/zhr/zhr.mate1.zhr.mosaik.bed.lifted.rmGaps.converted,var2=../mel_mel_data/mRNA-Seq/zhr/zhr.mate1.z30.mosaik.bed.lifted.rmGaps.converted,var3=../mel_mel_data/mRNA-Seq/zhr/zhr.mate2.zhr.mosaik.bed.lifted.rmGaps.converted,var4=../mel_mel_data/mRNA-Seq/zhr/zhr.mate2.z30.mosaik.bed.lifted.rmGaps.converted merge_reads.sh

#z30
#qsub -v var1=../mel_mel_data/mRNA-Seq/z30/z30.mate1.zhr.mosaik.bed.lifted.rmGaps.converted,var2=../mel_mel_data/mRNA-Seq/z30/z30.mate1.z30.mosaik.bed.lifted.rmGaps.converted,var3=../mel_mel_data/mRNA-Seq/z30/z30.mate2.zhr.mosaik.bed.lifted.rmGaps.converted,var4=../mel_mel_data/mRNA-Seq/z30/z30.mate2.z30.mosaik.bed.lifted.rmGaps.converted merge_reads.sh

#zhrXz30
#qsub -v var1=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.zhr.mosaik.bed.lifted.rmGaps.converted,var2=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.z30.mosaik.bed.lifted.rmGaps.converted,var3=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.zhr.mosaik.bed.lifted.rmGaps.converted,var4=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.z30.mosaik.bed.lifted.rmGaps.converted merge_reads.sh

#z30Xzhr
#qsub -v var1=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.zhr.mosaik.bed.lifted.rmGaps.converted,var2=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.z30.mosaik.bed.lifted.rmGaps.converted,var3=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.zhr.mosaik.bed.lifted.rmGaps.converted,var4=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.z30.mosaik.bed.lifted.rmGaps.converted merge_reads.sh