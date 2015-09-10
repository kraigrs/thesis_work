#!/bin/sh

# use this on a SAMtools pileup file to calculate average depth of coverage

#sh depth_of_coverage <samtools pileup>

awk '{ count++ ; SUM += $4 } END { print "Total: " SUM "\t" "Nucleotides: " count "\t" "Average_coverage: " SUM/count }' $1