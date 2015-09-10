#!/bin/sh

perl ../perl/equal_total_reads.pl single 50 0 ../../../Dmel/chromFa/chr2L.fa ../../../Dmel/chromFa/chr2L_alt.fa ../../sim_regs_chr2L.bed ../../zhr_z30_exons_expression.txt ../../equal_total

perl ../perl/equal_total_reads.pl single 50 0 ../../../Dmel/chromFa/chr2R.fa ../../../Dmel/chromFa/chr2R_alt.fa ../../sim_regs_chr2R.bed ../../zhr_z30_exons_expression.txt ../../equal_total

perl ../perl/equal_total_reads.pl single 50 0 ../../../Dmel/chromFa/chr3L.fa ../../../Dmel/chromFa/chr3L_alt.fa ../../sim_regs_chr3L.bed ../../zhr_z30_exons_expression.txt ../../equal_total

perl ../perl/equal_total_reads.pl single 50 0 ../../../Dmel/chromFa/chr3R.fa ../../../Dmel/chromFa/chr3R_alt.fa ../../sim_regs_chr3R.bed ../../zhr_z30_exons_expression.txt ../../equal_total

perl ../perl/equal_total_reads.pl single 50 0 ../../../Dmel/chromFa/chrX.fa ../../../Dmel/chromFa/chrX_alt.fa ../../sim_regs_chrX.bed ../../zhr_z30_exons_expression.txt ../../equal_total

cat ../../equal_total/*.fa > ../../equal_total/constExons_single_bp50_error0_equal_total.fa