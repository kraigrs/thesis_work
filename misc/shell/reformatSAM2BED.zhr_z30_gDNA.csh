#!/bin/csh

#zhr
qsub -v var1=../mel_mel_data/Resequencing/zhr/zhr.mate1.zhr.mosaik.sam,var2=yes,var3=../mel_mel_data/Resequencing/resequencing-assembly/zhr/zhr_AMBIG_Sep.chain reformatSAM2BED.0.3.csh
qsub -v var1=../mel_mel_data/Resequencing/zhr/zhr.mate2.zhr.mosaik.sam,var2=yes,var3=../mel_mel_data/Resequencing/resequencing-assembly/zhr/zhr_AMBIG_Sep.chain reformatSAM2BED.0.3.csh
qsub -v var1=../mel_mel_data/Resequencing/zhr/zhr.mate1.z30.mosaik.sam,var2=yes,var3=../mel_mel_data/Resequencing/resequencing-assembly/z30/z30_AMBIG_Sep.chain reformatSAM2BED.0.3.csh
qsub -v var1=../mel_mel_data/Resequencing/zhr/zhr.mate2.z30.mosaik.sam,var2=yes,var3=../mel_mel_data/Resequencing/resequencing-assembly/z30/z30_AMBIG_Sep.chain reformatSAM2BED.0.3.csh

#z30
qsub -v var1=../mel_mel_data/Resequencing/z30/z30.mate1.zhr.mosaik.sam,var2=yes,var3=../mel_mel_data/Resequencing/resequencing-assembly/zhr/zhr_AMBIG_Sep.chain reformatSAM2BED.0.3.csh
qsub -v var1=../mel_mel_data/Resequencing/z30/z30.mate2.zhr.mosaik.sam,var2=yes,var3=../mel_mel_data/Resequencing/resequencing-assembly/zhr/zhr_AMBIG_Sep.chain reformatSAM2BED.0.3.csh
qsub -v var1=../mel_mel_data/Resequencing/z30/z30.mate1.z30.mosaik.sam,var2=yes,var3=../mel_mel_data/Resequencing/resequencing-assembly/z30/z30_AMBIG_Sep.chain reformatSAM2BED.0.3.csh
qsub -v var1=../mel_mel_data/Resequencing/z30/z30.mate2.z30.mosaik.sam,var2=yes,var3=../mel_mel_data/Resequencing/resequencing-assembly/z30/z30_AMBIG_Sep.chain reformatSAM2BED.0.3.csh
