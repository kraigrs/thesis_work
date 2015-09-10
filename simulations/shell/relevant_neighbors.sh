

intersectBed -a sim_regs.36b_windows.bed -b ../DGRP/DGRP_line_40_SNPs.bed -c > sim_regs.36b_windows_SNPs.bed
intersectBed -a sim_regs.36b_windows_SNPs.bed -b ../DGRP/DGRP_line_40_SNPs.bed -wa -wb > sim_regs.36b_windows_SNPs.txt

intersectBed -a sim_regs.50b_windows.bed -b ../DGRP/DGRP_line_40_SNPs.bed -c > sim_regs.50b_windows_SNPs.bed
intersectBed -a sim_regs.50b_windows_SNPs.bed -b ../DGRP/DGRP_line_40_SNPs.bed -wa -wb > sim_regs.50b_windows_SNPs.txt


intersectBed -a berlin-updated-exonic-regions.36b_windows.bed -b berlin-updated-exonic-regions_fsa.SNPs.bed -c > berlin-updated-exonic-regions_fsa.36b_windows_SNPs.bed

intersectBed -a berlin-updated-exonic-regions_fsa.36b_windows_SNPs.bed -b berlin-updated-exonic-regions_fsa.SNPs.bed -wa -wb > berlin-updated-exonic-regions_fsa.36b_windows_SNPs.txt