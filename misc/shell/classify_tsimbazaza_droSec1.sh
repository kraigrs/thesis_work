#!/bin/sh

#tsimbazaza

qsub -v var1=../McManus/constitutive_regions_overlap_filtered.bed,var2=tsimbazaza,var3=droSec1,var4=../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza_info_rmDup.txt,var5=../sim_sec_data/mRNA-Seq/tsimbazaza,var6=tsimbazaza_v2 classify.sh

#droSec1

qsub -v var1=../McManus/constitutive_regions_overlap_filtered.bed,var2=tsimbazaza,var3=droSec1,var4=../sim_sec_data/mRNA-Seq/droSec1/droSec1_info_rmDup.txt,var5=../sim_sec_data/mRNA-Seq/droSec1,var6=droSec1_v2 classify.sh

#tsimbazaza+droSec1

qsub -v var1=../McManus/constitutive_regions_overlap_filtered.bed,var2=tsimbazaza,var3=droSec1,var4=../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1_info_rmDup.txt,var5=../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1,var6=tsimbazaza+droSec1_v2 classify.sh

#tsimbazazaXdroSec1

qsub -v var1=../McManus/constitutive_regions_overlap_filtered.bed,var2=tsimbazaza,var3=droSec1,var4=../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1_info_rmDup.txt,var5=../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1,var6=tsimbazazaXdroSec1_v2 classify.sh

#droSec1Xtsimbazaza

qsub -v var1=../McManus/constitutive_regions_overlap_filtered.bed,var2=tsimbazaza,var3=droSec1,var4=../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza_info_rmDup.txt,var5=../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza,var6=droSec1Xtsimbazaza_v2 classify.sh
