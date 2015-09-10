#!/bin/csh

#qsub -v var1=../mel_sec_data/Hyb_r1.droSec1.mosaik.bed,var2=../Dsec/droSec1ToDm3.over.chain,var3=../mel_sec_data/Hyb_r1.droSec1.mosaik.bed.lifted,var4=../mel_sec_data/Hyb_r1.droSec1.mosaik.bed.un liftOver.csh

#qsub -v var1=../mel_sec_data/Hyb_r2.droSec1.mosaik.bed,var2=../Dsec/droSec1ToDm3.over.chain,var3=../mel_sec_data/Hyb_r2.droSec1.mosaik.bed.lifted,var4=../mel_sec_data/Hyb_r2.droSec1.mosaik.bed.un liftOver.csh

#qsub -v var1=../mel_sec_data/Mix_r1.droSec1.mosaik.bed,var2=../Dsec/droSec1ToDm3.over.chain,var3=../mel_sec_data/Mix_r1.droSec1.mosaik.bed.lifted,var4=../mel_sec_data/Mix_r1.droSec1.mosaik.bed.un liftOver.csh

#qsub -v var1=../mel_sec_data/Mix_r2.droSec1.mosaik.bed,var2=../Dsec/droSec1ToDm3.over.chain,var3=../mel_sec_data/Mix_r2.droSec1.mosaik.bed.lifted,var4=../mel_sec_data/Mix_r2.droSec1.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_sec_data/mRNA-Seq/dm3/dm3.mate1.droSec1.mosaik.bed,var2=../Dsec/droSec1ToDm3.over.chain,var3=../mel_sec_data/mRNA-Seq/dm3/dm3.mate1.droSec1.mosaik.bed.lifted,var4=../mel_sec_data/mRNA-Seq/dm3/dm3.mate1.droSec1.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_sec_data/mRNA-Seq/dm3/dm3.mate2.droSec1.mosaik.bed,var2=../Dsec/droSec1ToDm3.over.chain,var3=../mel_sec_data/mRNA-Seq/dm3/dm3.mate2.droSec1.mosaik.bed.lifted,var4=../mel_sec_data/mRNA-Seq/dm3/dm3.mate2.droSec1.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_sec_data/mRNA-Seq/droSec1/droSec1.mate1.droSec1.mosaik.bed,var2=../Dsec/droSec1ToDm3.over.chain,var3=../mel_sec_data/mRNA-Seq/droSec1/droSec1.mate1.droSec1.mosaik.bed.lifted,var4=../mel_sec_data/mRNA-Seq/droSec1/droSec1.mate1.droSec1.mosaik.bed.un liftOver.csh

qsub -v var1=../mel_sec_data/mRNA-Seq/droSec1/droSec1.mate2.droSec1.mosaik.bed,var2=../Dsec/droSec1ToDm3.over.chain,var3=../mel_sec_data/mRNA-Seq/droSec1/droSec1.mate2.droSec1.mosaik.bed.lifted,var4=../mel_sec_data/mRNA-Seq/droSec1/droSec1.mate2.droSec1.mosaik.bed.un liftOver.csh
