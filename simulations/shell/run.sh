#!/bin/sh

# mel_mel

#var1=../../../Dmel/dm3_ref.fa var2=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.fa qsub bowtiePipe_3mm.sh
#var1=../../../Dmel/dm3_ref.fa var2=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.sam qsub samtoolsPipe.sh

#var1=../../tiled/line_40/DGRP_line_40_SNPs_const.txt var2=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_ref.bowtie_v0_m1.pileup.txt var3=36 var4=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_ref.bowtie_v0_m1.SNPs.txt qsub pileup_SNP_ASE.sh
#var1=../../tiled/line_40/DGRP_line_40_SNPs_const.txt var2=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.pileup.txt var3=36 var4=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.SNPs.txt qsub pileup_SNP_ASE.sh
#var1=../../tiled/line_40/DGRP_line_40_SNPs_const.txt var2=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.pileup.txt var3=36 var4=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.SNPs.txt qsub pileup_SNP_ASE.sh
#var1=../../tiled/line_40/DGRP_line_40_SNPs_const.txt var2=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.pileup.txt var3=36 var4=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.SNPs.txt qsub pileup_SNP_ASE.sh

#var1=../../tiled/line_40/DGRP_line_40_SNPs_const.txt var2=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_alt_line_40.bowtie_v0_m1.pileup.txt var3=36 var4=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_alt_line_40.bowtie_v0_m1.SNPs.txt qsub pileup_SNP_ASE.sh
#var1=../../tiled/line_40/DGRP_line_40_SNPs_const.txt var2=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_alt_line_40.bowtie_v1_m1.pileup.txt var3=36 var4=../../tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_alt_line_40.bowtie_v1_m1.SNPs.txt qsub pileup_SNP_ASE.sh

# mel_sim

#var1=../../../Graze/berlin-updated-exonic-regions.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.fa qsub bowtiePipe_0mm.sh
#var1=../../../Graze/berlin-updated-exonic-regions_fsa_masked.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.fa qsub bowtiePipe_0mm.sh
#var1=../../../Graze/berlin-updated-exonic-regions.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.fa qsub bowtiePipe_1mm.sh
#var1=../../../Graze/berlin-updated-exonic-regions.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.fa qsub bowtiePipe_2mm.sh
#var1=../../../Graze/berlin-updated-exonic-regions.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.fa qsub bowtiePipe_3mm.sh

#var1=../../../Graze/c1674-updated-exonic-regions.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.fa qsub bowtiePipe_0mm.sh
#var1=../../../Graze/c1674-updated-exonic-regions_fsa_masked.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.fa qsub bowtiePipe_0mm.sh
#var1=../../../Graze/c1674-updated-exonic-regions.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.fa qsub bowtiePipe_1mm.sh

#var1=../../../Graze/berlin-updated-exonic-regions.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v0_m1.sam qsub samtoolsPipe.sh
#var1=../../../Graze/berlin-updated-exonic-regions_fsa_masked.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions_fsa_masked.bowtie_v0_m1.sam qsub samtoolsPipe.sh
#var1=../../../Graze/berlin-updated-exonic-regions.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v1_m1.sam qsub samtoolsPipe.sh
#var1=../../../Graze/berlin-updated-exonic-regions.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v2_m1.sam qsub samtoolsPipe.sh

#var1=../../../Graze/c1674-updated-exonic-regions.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.c1674-updated-exonic-regions.bowtie_v0_m1.sam qsub samtoolsPipe.sh
#var1=../../../Graze/c1674-updated-exonic-regions_fsa_masked.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.c1674-updated-exonic-regions_fsa_masked.bowtie_v0_m1.sam qsub samtoolsPipe.sh
#var1=../../../Graze/c1674-updated-exonic-regions.fasta var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.c1674-updated-exonic-regions.bowtie_v1_m1.sam qsub samtoolsPipe.sh

#var1=../../../Graze/simulation/berlin_c1674_fsa.SNPs.txt var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions_fsa_masked.bowtie_v0_m1.pileup.txt var3=36 var4=../../../Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions_fsa_masked.bowtie_v0_m1.SNPs.txt qsub pileup_SNP_ASE.sh

#var1=../../../Graze/simulation/berlin-updated-exonic-regions_fsa.SNPs.txt var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v1_m1.pileup.txt var3=36 var4=../../../Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v1_m1.SNPs.txt qsub pileup_SNP_ASE.sh
#var1=../../../Graze/simulation/berlin-updated-exonic-regions_fsa.SNPs.txt var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v2_m1.pileup.txt var3=36 var4=../../../Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v2_m1.SNPs.txt qsub pileup_SNP_ASE.sh
#var1=../../../Graze/simulation/berlin-updated-exonic-regions_fsa.SNPs.txt var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v3_m1.pileup.txt var3=36 var4=../../../Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v3_m1.SNPs.txt qsub pileup_SNP_ASE.sh

#var1=../../../Graze/simulation/berlin_c1674_fsa.SNPs.txt var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.c1674-updated-exonic-regions_fsa_masked.bowtie_v0_m1.pileup.txt var3=36 var4=../../../Graze/simulation/berlin_c1674.tiled_36bp.c1674-updated-exonic-regions_fsa_masked.bowtie_v0_m1.SNPs.txt qsub pileup_SNP_ASE.sh

#var1=../../../Graze/simulation/c1674-updated-exonic-regions_fsa.SNPs.txt var2=../../../Graze/simulation/berlin_c1674.tiled_36bp.c1674-updated-exonic-regions.bowtie_v1_m1.pileup.txt var3=36 var4=../../../Graze/simulation/berlin_c1674.tiled_36bp.c1674-updated-exonic-regions.bowtie_v1_m1.SNPs.txt qsub pileup_SNP_ASE.sh



# mappability

#var1=../../../Dmel/dm3_ref.fa var2=../../../Dmel/dm3_ref qsub gem.sh
#var1=../../../Dmel/dm3_alt_line_40.fa var2=../../../Dmel/dm3_alt_line_40 qsub gem.sh

#var1=../../../Dmel/dm3_ref var2=0 var3=../../../Dmel/dm3_ref.l36_m0.mappability.txt qsub gem.sh
#var1=../../../Dmel/dm3_ref var2=1 var3=../../../Dmel/dm3_ref.l36_m1.mappability.txt qsub gem.sh
#var1=../../../Dmel/dm3_ref var2=2 var3=../../../Dmel/dm3_ref.l36_m2.mappability.txt qsub gem.sh
#var1=../../../Dmel/dm3_ref var2=3 var3=../../../Dmel/dm3_ref.l36_m3.mappability.txt qsub gem.sh

#var1=../../../Dmel/dm3_ref var2=4 var3=../../../Dmel/dm3_ref.l50_m4.mappability.txt qsub gem.sh

#var1=../../../Dmel/dm3_alt_line_40 var2=0 var3=../../../Dmel/dm3_alt_line_40.l36_m0.mappability.txt qsub gem.sh
#var1=../../../Dmel/dm3_alt_line_40 var2=1 var3=../../../Dmel/dm3_alt_line_40.l36_m1.mappability.txt qsub gem.sh
#var1=../../../Dmel/dm3_alt_line_40 var2=2 var3=../../../Dmel/dm3_alt_line_40.l36_m2.mappability.txt qsub gem.sh
#var1=../../../Dmel/dm3_alt_line_40 var2=3 var3=../../../Dmel/dm3_alt_line_40.l36_m3.mappability.txt qsub gem.sh



# Bowtie -n

#var1=../../tiled/line_40/DGRP_line_40_SNPs_const.txt var2=../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_n3_e161_l36_m1.pileup.txt var3=36 var4=../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_n3_e161_l36_m1.SNPs.txt qsub pileup_SNP_ASE.sh

# GSNAP

#var1=../../tiled/line_40/DGRP_line_40_SNPs_const.txt var2=../../../GSNAP/constExons_single_bp36_error0_tiled_line_40.gsnap.pileup.txt var3=36 var4=../../../GSNAP/constExons_single_bp36_error0_tiled_line_40.gsnap.SNPs.txt qsub pileup_SNP_ASE.sh

#var1=../../tiled/line_40/DGRP_line_40_SNPs_const.txt var2=../../../GSNAP/constExons_single_bp50_error0_tiled_line_40.gsnap.pileup.txt var3=50 var4=../../../GSNAP/constExons_single_bp50_error0_tiled_line_40.gsnap.SNPs.txt qsub pileup_SNP_ASE.sh


#var1=../../../Dmel/dm3_ref.fa var2=../../../GSNAP/constExons_single_bp36_error0_tiled_line_40.gsnap_unique.sam qsub samtoolsPipe.sh

#var1=../../tiled/line_40/DGRP_line_40_SNPs_const.txt var2=../../../GSNAP/constExons_single_bp36_error0_tiled_line_40.gsnap_unique.pileup.txt var3=36 var4=../../../GSNAP/constExons_single_bp36_error0_tiled_line_40.gsnap_unique.SNPs.txt qsub pileup_SNP_ASE.sh