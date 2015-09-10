#!/bin/csh

#zhr

qsub -v var1=../mel_mel_data/mRNA-Seq/zhr/zhr.mate1.zhr.mosaik.bed,var2=../mel_mel_data/Resequencing/zhr/zhr_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/zhr/zhr.mate1.zhr.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/zhr/zhr.mate1.zhr.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_mel_data/mRNA-Seq/zhr/zhr.mate2.zhr.mosaik.bed,var2=../mel_mel_data/Resequencing/zhr/zhr_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/zhr/zhr.mate2.zhr.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/zhr/zhr.mate2.zhr.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_mel_data/mRNA-Seq/zhr/zhr.mate1.z30.mosaik.bed,var2=../mel_mel_data/Resequencing/z30/z30_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/zhr/zhr.mate1.z30.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/zhr/zhr.mate1.z30.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_mel_data/mRNA-Seq/zhr/zhr.mate2.z30.mosaik.bed,var2=../mel_mel_data/Resequencing/z30/z30_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/zhr/zhr.mate2.z30.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/zhr/zhr.mate2.z30.mosaik.bed.un liftOver.csh

#z30

qsub -v var1=../mel_mel_data/mRNA-Seq/z30/z30.mate1.zhr.mosaik.bed,var2=../mel_mel_data/Resequencing/zhr/zhr_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/z30/z30.mate1.zhr.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/z30/z30.mate1.zhr.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_mel_data/mRNA-Seq/z30/z30.mate2.zhr.mosaik.bed,var2=../mel_mel_data/Resequencing/zhr/zhr_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/z30/z30.mate2.zhr.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/z30/z30.mate2.zhr.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_mel_data/mRNA-Seq/z30/z30.mate1.z30.mosaik.bed,var2=../mel_mel_data/Resequencing/z30/z30_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/z30/z30.mate1.z30.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/z30/z30.mate1.z30.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_mel_data/mRNA-Seq/z30/z30.mate2.z30.mosaik.bed,var2=../mel_mel_data/Resequencing/z30/z30_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/z30/z30.mate2.z30.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/z30/z30.mate2.z30.mosaik.bed.un liftOver.csh

#zhrXz30

qsub -v var1=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.zhr.mosaik.bed,var2=../mel_mel_data/Resequencing/zhr/zhr_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.zhr.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.zhr.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.zhr.mosaik.bed,var2=../mel_mel_data/Resequencing/zhr/zhr_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.zhr.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.zhr.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.z30.mosaik.bed,var2=../mel_mel_data/Resequencing/z30/z30_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.z30.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.z30.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.z30.mosaik.bed,var2=../mel_mel_data/Resequencing/z30/z30_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.z30.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.z30.mosaik.bed.un liftOver.csh

#z30Xzhr

qsub -v var1=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.zhr.mosaik.bed,var2=../mel_mel_data/Resequencing/zhr/zhr_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.zhr.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.zhr.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.zhr.mosaik.bed,var2=../mel_mel_data/Resequencing/zhr/zhr_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.zhr.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.zhr.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.z30.mosaik.bed,var2=../mel_mel_data/Resequencing/z30/z30_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.z30.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.z30.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.z30.mosaik.bed,var2=../mel_mel_data/Resequencing/z30/z30_AMBIG_all_dm3.chain,var3=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.z30.mosaik.bed.lifted,var4=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.z30.mosaik.bed.un liftOver.csh
