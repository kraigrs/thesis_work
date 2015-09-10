#!/bin/sh

#tsimbazaza

qsub -v var1=../McManus/constitutive_regions_overlap_filtered.bed,var2=tsimbazaza,var3=droSec1,var4=../sim_sec_data/Resequencing/tsimbazaza/tsimbazaza_gDNA_rmDup.txt,var5=../sim_sec_data/Resequencing/tsimbazaza,var6=tsimbazaza_gDNA_v2 classify.sh

#droSec1

qsub -v var1=../McManus/constitutive_regions_overlap_filtered.bed,var2=tsimbazaza,var3=droSec1,var4=../sim_sec_data/Resequencing/droSec1/droSec1_gDNA_rmDup.txt,var5=../sim_sec_data/Resequencing/droSec1,var6=droSec1_gDNA_v2 classify.sh
