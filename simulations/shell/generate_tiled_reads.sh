#!/bin/sh

perl ../perl/tiled_reads.pl single 36 0 ../../../Dmel/chromFa/chr2L.fa ../../../Dmel/chr2L_alt_line_40.fa ../../sim_regs_chr2L.bed ../../tiled/line_40

perl ../perl/tiled_reads.pl single 36 0 ../../../Dmel/chromFa/chr2R.fa ../../../Dmel/chr2R_alt_line_40.fa ../../sim_regs_chr2R.bed ../../tiled/line_40

perl ../perl/tiled_reads.pl single 36 0 ../../../Dmel/chromFa/chr3L.fa ../../../Dmel/chr3L_alt_line_40.fa ../../sim_regs_chr3L.bed ../../tiled/line_40

perl ../perl/tiled_reads.pl single 36 0 ../../../Dmel/chromFa/chr3R.fa ../../../Dmel/chr3R_alt_line_40.fa ../../sim_regs_chr3R.bed ../../tiled/line_40

perl ../perl/tiled_reads.pl single 36 0 ../../../Dmel/chromFa/chrX.fa ../../../Dmel/chrX_alt_line_40.fa ../../sim_regs_chrX.bed ../../tiled/line_40

#cat ../../tiled/*.fa > ../../tiled/constExons_single_bp50_error0_tiled.fa

perl ../perl/tiled_reads.pl single 100 0 ../../../Dmel/chromFa/chr2L.fa ../../../Dmel/chr2L_alt_line_40.fa ../../sim_regs_chr2L.bed ../../tiled/line_40

perl ../perl/tiled_reads.pl single 100 0 ../../../Dmel/chromFa/chr2R.fa ../../../Dmel/chr2R_alt_line_40.fa ../../sim_regs_chr2R.bed ../../tiled/line_40

perl ../perl/tiled_reads.pl single 100 0 ../../../Dmel/chromFa/chr3L.fa ../../../Dmel/chr3L_alt_line_40.fa ../../sim_regs_chr3L.bed ../../tiled/line_40

perl ../perl/tiled_reads.pl single 100 0 ../../../Dmel/chromFa/chr3R.fa ../../../Dmel/chr3R_alt_line_40.fa ../../sim_regs_chr3R.bed ../../tiled/line_40

perl ../perl/tiled_reads.pl single 100 0 ../../../Dmel/chromFa/chrX.fa ../../../Dmel/chrX_alt_line_40.fa ../../sim_regs_chrX.bed ../../tiled/line_40

#cat ../../tiled/*.fa > ../../tiled/constExons_single_bp50_error0_tiled.fa