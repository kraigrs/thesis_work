#!/bin/sh

# Totals for dm3
# introns --> 59376
# genes with introns --> 10132
# genes --> 13913

# get introns that have whole genes in them 
# introns --> 1145
# genes with introns --> 617

intersectBed -a dm3_genes.bed -b dm3_introns.bed -wa -wb -f 1 | cut -f7-12 | sort -u > dm3_introns_with_genes.bed
cut -f4 dm3_introns_with_genes.bed > file


# get introns without genes in them
# introns --> 58227
# genes with introns --> 10116



# remove overlapping genes from intron list
# 5351/13913 annotated genes have some amount of overlap

intersectBed -a dm3_genes.bed -b dm3_genes.bed -wa -wb | cut -f4,10 | awk '$1 != $2' | cut -f1 | sort -u | wc -l
intersectBed -a dm3_genes.bed -b dm3_genes.bed -wa -wb | cut -f4,10 | awk '$1 != $2' | cut -f1 | sort -u > remove_genes.txt



# get intergenic regions only from candidate genes



# subtract coding information from introns

subtractBed -a dm3_introns.bed -b dm3_exons.bed > dm3_introns_no_exons.bed


# exon information for regulatory evolution 

