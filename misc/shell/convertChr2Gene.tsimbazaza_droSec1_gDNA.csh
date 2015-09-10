#!/bin/csh

#tsimbazaza
qsub -v var1=../sim_sec_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate1.tsimbazaza.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate2.tsimbazaza.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate1.droSec1.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mate2.droSec1.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh

#droSec1
qsub -v var1=../sim_sec_data/Resequencing/droSec1/droSec1_gDNA.mate1.tsimbazaza.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/Resequencing/droSec1/droSec1_gDNA.mate2.tsimbazaza.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/Resequencing/droSec1/droSec1_gDNA.mate1.droSec1.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/Resequencing/droSec1/droSec1_gDNA.mate2.droSec1.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
