#!/bin/sh

#tsimbazaza
qsub -v var1=../sim_sec_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate1.tsimbazaza.mosaik.bed.lifted.converted,var2=../sim_sec_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate1.droSec1.mosaik.bed.lifted.converted,var3=../sim_sec_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate2.tsimbazaza.mosaik.bed.lifted.converted,var4=../sim_sec_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate2.droSec1.mosaik.bed.lifted.converted merge_reads.sh

#droSec1
qsub -v var1=../sim_sec_data/Resequencing/droSec1/droSec1_gDNA.mate1.tsimbazaza.mosaik.bed.lifted.converted,var2=../sim_sec_data/Resequencing/droSec1/droSec1_gDNA.mate1.droSec1.mosaik.bed.lifted.converted,var3=../sim_sec_data/Resequencing/droSec1/droSec1_gDNA.mate2.tsimbazaza.mosaik.bed.lifted.converted,var4=../sim_sec_data/Resequencing/droSec1/droSec1_gDNA.mate2.droSec1.mosaik.bed.lifted.converted merge_reads.sh
