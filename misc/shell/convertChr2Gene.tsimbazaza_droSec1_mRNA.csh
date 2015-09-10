#!/bin/csh

#tsimbazaza
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate1.tsimbazaza.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate2.tsimbazaza.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate1.droSec1.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate2.droSec1.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh

#droSec1
qsub -v var1=../sim_sec_data/mRNA-Seq/droSec1/droSec1.mate1.tsimbazaza.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/droSec1/droSec1.mate2.tsimbazaza.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/droSec1/droSec1.mate1.droSec1.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/droSec1/droSec1.mate2.droSec1.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh

#tsimbazaza+droSec1
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1.mate1.tsimbazaza.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1.mate2.tsimbazaza.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1.mate1.droSec1.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1.mate2.droSec1.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh

#tsimbazazaXdroSec1
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1.mate1.droSec1.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1.mate2.droSec1.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1.mate1.tsimbazaza.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1.mate2.tsimbazaza.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh

#droSec1Xtsimbazaza
qsub -v var1=../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza.mate1.droSec1.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza.mate2.droSec1.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza.mate1.tsimbazaza.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
qsub -v var1=../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza.mate2.tsimbazaza.mosaik.bed.lifted,var2=../sim_sec_data/Sim_sec_gaps_gt19.bed,var3=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.5.csh
