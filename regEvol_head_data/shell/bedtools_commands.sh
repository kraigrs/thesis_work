#!/bin/sh

# exon information for regulatory evolution 

# find overlapping portions of exons

intersectBed -a ../../references/dm3_exons.bed -b ../../references/dm3_exons.bed -wa -wb > ../../references/dm3_exons.bed.intersected

# get exons that overlap within the same gene

perl ../perl/get_overlaps.pl ../../references/dm3_exons.bed > ../../references/dm3_exons_overlap_within_gene.bed

sort -u ../../references/dm3_exons_overlap_within_gene.bed > ../../references/temp
mv ../../references/temp ../../references/dm3_exons_overlap_within_gene.bed

mergeBed -s -nms -i dm3_exons_overlap_within_gene.bed > dm3_exons_overlap_within_gene_merged.bed