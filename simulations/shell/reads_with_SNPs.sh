#!/bin/sh

# v0_m1

perl ../perl/makeBED_reads.pl ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v0_m1.max ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v0_m1.max.bed

intersectBed -a ../../../DGRP/DGRP_line_40_SNPs.bed -b ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v0_m1.max.bed -f 1 -wb > ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v0_m1.max.SNPs.txt

perl ../perl/makeBED_reads.pl ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v0_m1.un ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v0_m1.un.bed

intersectBed -a ../../../DGRP/DGRP_line_40_SNPs.bed -b ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v0_m1.un.bed -f 1 -wb > ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v0_m1.un.SNPs.txt

# v1_m1

perl ../perl/makeBED_reads.pl ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.max ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.max.bed

intersectBed -a ../../../DGRP/DGRP_line_40_SNPs.bed -b ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.max.bed -f 1 -wb > ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.max.SNPs.txt

perl ../perl/makeBED_reads.pl ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.un ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.un.bed

intersectBed -a ../../../DGRP/DGRP_line_40_SNPs.bed -b ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.un.bed -f 1 -wb > ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.un.SNPs.txt

# v2_m1

perl ../perl/makeBED_reads.pl ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.max ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.max.bed

intersectBed -a ../../../DGRP/DGRP_line_40_SNPs.bed -b ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.max.bed -f 1 -wb > ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.max.SNPs.txt

perl ../perl/makeBED_reads.pl ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.un ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.un.bed

intersectBed -a ../../../DGRP/DGRP_line_40_SNPs.bed -b ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.un.bed -f 1 -wb > ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.un.SNPs.txt

# v3_m1

perl ../perl/makeBED_reads.pl ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.max ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.max.bed

intersectBed -a ../../../DGRP/DGRP_line_40_SNPs.bed -b ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.max.bed -f 1 -wb > ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.max.SNPs.txt

perl ../perl/makeBED_reads.pl ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.un ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.un.bed

intersectBed -a ../../../DGRP/DGRP_line_40_SNPs.bed -b ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.un.bed -f 1 -wb > ../../tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.un.SNPs.txt