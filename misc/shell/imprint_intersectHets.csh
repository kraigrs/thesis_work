#!/bin/csh -f

foreach chr (2L 2LHet 2R 3L 3LHet 3R 4 U X XHet YHet dmel_mitochondrion_genome)

intersectBed -a ../mel_mel_data/imprinting/chr$chr/chr$chr\_windows_dist.bed -b ../mel_mel_data/HetSites4imprinting/hetSNPs.bed -c > ../mel_mel_data/imprinting/chr$chr/chr$chr\_hetSNPs.bed

intersectBed -a ../mel_mel_data/imprinting/chr$chr/chr$chr\_windows_dist.bed -b ../mel_mel_data/HetSites4imprinting/zhr_rsq2_goodsnps2_hetindels.bed -c > ../mel_mel_data/imprinting/chr$chr/chr$chr\_hetINDELs.bed

intersectBed -a ../mel_mel_data/imprinting/chr$chr/chr$chr\_windows_dist.bed -b ../mel_mel_data/HetSites4imprinting/hetSNPs.bed -c > ../mel_mel_data/imprinting/chr$chr/chr$chr\_hetSNPs.bed

intersectBed -a ../mel_mel_data/imprinting/chr$chr/chr$chr\_windows_dist.bed -b ../mel_mel_data/HetSites4imprinting/z30_rsq2_goodsnps2_hetindels.bed -c > ../mel_mel_data/imprinting/chr$chr/chr$chr\_hetINDELs.bed

end
