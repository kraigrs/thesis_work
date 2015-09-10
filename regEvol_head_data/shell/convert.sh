#!/bin/sh

FILES=../../*/*.bed.lifted
for f in $FILES
do
  var1=$f var2=../../../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../../../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.sh
done