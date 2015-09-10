#!/bin/sh

#zhr
#qsub -v var1=../mel_mel_data/mRNA-Seq/zhr/zhr.mate1.zhr.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh
#qsub -v var1=../mel_mel_data/mRNA-Seq/zhr/zhr.mate2.zhr.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh
#qsub -v var1=../mel_mel_data/mRNA-Seq/zhr/zhr.mate1.z30.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh
#qsub -v var1=../mel_mel_data/mRNA-Seq/zhr/zhr.mate2.z30.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh

#z30
#qsub -v var1=../mel_mel_data/mRNA-Seq/z30/z30.mate1.zhr.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh
#qsub -v var1=../mel_mel_data/mRNA-Seq/z30/z30.mate2.zhr.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh
#qsub -v var1=../mel_mel_data/mRNA-Seq/z30/z30.mate1.z30.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh
#qsub -v var1=../mel_mel_data/mRNA-Seq/z30/z30.mate2.z30.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh

#zhrXz30
#qsub -v var1=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.zhr.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh
#qsub -v var1=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.zhr.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh
#qsub -v var1=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.z30.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh
#qsub -v var1=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.z30.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh

#z30Xzhr
#qsub -v var1=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.zhr.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh
#qsub -v var1=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.zhr.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh
#qsub -v var1=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.z30.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh
#qsub -v var1=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.z30.mosaik.bed.lifted,var2=../McManus/constitutive_regions_overlap_filtered.bed convertChr2Gene.0.6.sh

# gap filtering

#zhr
var1=../mel_mel_data/mRNA-Seq/zhr/zhr.mate1.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/mRNA-Seq/zhr/zhr.mate2.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/mRNA-Seq/zhr/zhr.mate1.z30.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/mRNA-Seq/zhr/zhr.mate2.z30.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh

#z30
var1=../mel_mel_data/mRNA-Seq/z30/z30.mate1.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/mRNA-Seq/z30/z30.mate2.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/mRNA-Seq/z30/z30.mate1.z30.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/mRNA-Seq/z30/z30.mate2.z30.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh

#zhrXz30
var1=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.z30.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.z30.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh

#z30Xzhr
var1=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.zhr.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.z30.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
var1=../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.z30.mosaik.bed.lifted var2=../sim_sec_data/Sim_sec_gaps_gt19.bed var3=../McManus/constitutive_regions_overlap_filtered.bed qsub convertChr2Gene.0.5.sh
