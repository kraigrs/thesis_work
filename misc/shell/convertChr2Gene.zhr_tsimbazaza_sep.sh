#!/bin/sh

#zhr
perl ../perlScripts/convertChr2Gene.0.6.pl  ../mel_sim_data/mRNA-Seq/zhr/zhr.mate1.zhr.mosaik.bed.lifted ../McManus/constitutive_regions_overlap_filtered.bed
perl ../perlScripts/convertChr2Gene.0.6.pl  ../mel_sim_data/mRNA-Seq/zhr/zhr.mate2.zhr.mosaik.bed.lifted ../McManus/constitutive_regions_overlap_filtered.bed
perl ../perlScripts/convertChr2Gene.0.6.pl  ../mel_sim_data/mRNA-Seq/zhr/zhr.mate1.tsimbazaza.mosaik.bed.lifted ../McManus/constitutive_regions_overlap_filtered.bed
perl ../perlScripts/convertChr2Gene.0.6.pl  ../mel_sim_data/mRNA-Seq/zhr/zhr.mate2.tsimbazaza.mosaik.bed.lifted ../McManus/constitutive_regions_overlap_filtered.bed

#tsimbazaza
perl ../perlScripts/convertChr2Gene.0.6.pl  ../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate1.zhr.mosaik.bed.lifted ../McManus/constitutive_regions_overlap_filtered.bed
perl ../perlScripts/convertChr2Gene.0.6.pl  ../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate2.zhr.mosaik.bed.lifted ../McManus/constitutive_regions_overlap_filtered.bed
perl ../perlScripts/convertChr2Gene.0.6.pl  ../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate1.tsimbazaza.mosaik.bed.lifted ../McManus/constitutive_regions_overlap_filtered.bed
perl ../perlScripts/convertChr2Gene.0.6.pl  ../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate2.tsimbazaza.mosaik.bed.lifted ../McManus/constitutive_regions_overlap_filtered.bed
