#!/bin/sh

#tsimbazaza
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate1.tsimbazaza.mosaik.bed.lifted.converted,var2=../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate1.droSec1.mosaik.bed.lifted.converted,var3=../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate2.tsimbazaza.mosaik.bed.lifted.converted,var4=../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate2.droSec1.mosaik.bed.lifted.converted merge_reads.sh

#droSec1
qsub -v var1=../sim_sec_data/mRNA-Seq/droSec1/droSec1.mate1.tsimbazaza.mosaik.bed.lifted.converted,var2=../sim_sec_data/mRNA-Seq/droSec1/droSec1.mate1.droSec1.mosaik.bed.lifted.converted,var3=../sim_sec_data/mRNA-Seq/droSec1/droSec1.mate2.tsimbazaza.mosaik.bed.lifted.converted,var4=../sim_sec_data/mRNA-Seq/droSec1/droSec1.mate2.droSec1.mosaik.bed.lifted.converted merge_reads.sh

#tsimbazazaXdroSec1
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1.mate1.tsimbazaza.mosaik.bed.lifted.converted,var2=../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1.mate1.droSec1.mosaik.bed.lifted.converted,var3=../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1.mate2.tsimbazaza.mosaik.bed.lifted.converted,var4=../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1.mate2.droSec1.mosaik.bed.lifted.converted merge_reads.sh

#droSec1Xtsimbazaza
qsub -v var1=../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza.mate1.tsimbazaza.mosaik.bed.lifted.converted,var2=../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza.mate1.droSec1.mosaik.bed.lifted.converted,var3=../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza.mate2.tsimbazaza.mosaik.bed.lifted.converted,var4=../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza.mate2.droSec1.mosaik.bed.lifted.converted merge_reads.sh

#tsimbazaza+droSec1
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1.mate1.tsimbazaza.mosaik.bed.lifted.converted,var2=../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1.mate1.droSec1.mosaik.bed.lifted.converted,var3=../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1.mate2.tsimbazaza.mosaik.bed.lifted.converted,var4=../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1.mate2.droSec1.mosaik.bed.lifted.converted merge_reads.sh