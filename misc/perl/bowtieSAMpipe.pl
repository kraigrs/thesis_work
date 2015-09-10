#!/usr/bin/perl

########################################################################
# 
# 10/08/2009
#
# bowtieSAMpipe.pl
#
# Purpose: align mate paired-end short sequence reads using Bowtie and
#          collect SNP statistics using SAM
# 
# Input: mate paired ends to bowtie, then a SAM output to samtools
#
# Output: SNPs and their frequency
#
########################################################################

use strict;
use warnings;

#my$start = time;

print "\n\n";
print "################################################################\n";
print "#                                                              #\n";
print "# Welcome to Joe and Kraig's fabulous bioinformatics pipeline! #\n";
print "#                                                              #\n";
print "# This pipeline is meant to analyze high-throughput sequencing #\n";
print "# data from various species of Drosophila. It takes a set of   #\n";
print "# reads and aligns them to the Drosophila reference species of #\n";
print "# your choice. The reads of interest in developing this        #\n";
print "# are those derived from hybrid progeny of highly inbred       #\n";
print "# parental strains of the same species. In this case, each     #\n";
print "# offspring should have a maternal and paternal copy of each   #\n";
print "# species. By aligning these reads from hybrids, we can select #\n";
print "# trans-acting mutations out and only consider cis-acting      #\n";
print "# mutations. The relative expression levels can be gathered    #\n";
print "# according to the number of reads that align to a particular  #\n";
print "# segment of the reference genome. This allows us to identify  #\n";
print "# cis-regulatory differences between the two species.          #\n";
print "#                                                              #\n";
print "# The pipeline is organized in the following way:              #\n";
print "#                                                              #\n";
print "# Step 1: perform alignment using Bowtie (Langmead et al.)     #\n";
print "# Step 2: find variations (in this case, SNPs) using SAMtools  #\n";
print "# (Li et al.)                                                  #\n";
print "#                                                              #\n";
print "# You will now be prompted to answer a few questions in order  #\n";
print "# to tailor the alignment to your specific needs.              #\n";
print "#                                                              #\n";
print "#                       Good luck!                             #\n";
print "#                                                              #\n";
print "################################################################\n";
print "\n";

print "\nWhich species would you like to align your reads to? e.g. dmel, dsim, etc. --> ";
my$species = <STDIN>;

print "\nWhere is the reference file? --> ";
my$fasta = <STDIN>;
my$ebwt;
if($fasta =~ /\(\S+\/*\S+)\.fast[aq]$/){$ebwt = $1;} # takes the fasta file and finds the text for the index build
#print "\n$ebwt\n\n";

print "\nWould you like to use paired-end reads? (yes/no) --> ";
my$answer = <STDIN>;
if($answer =~ /[yesYES]/)
{
  print "\nWhere is the first mate-paired end reads file?  --> ";
  my$mate1 = <STDIN>;
  print "Where is the second mate-paired end reads file? --> ";
  my$mate2 = <STDIN>;
}
else
{
  print "\nWhere is the reads file? --> ";
  my$reads = <STDIN>;
}

print "\nWhat suffix would you like to use for all your files? i.e. s_2_sequence --> ";
my$output = <STDIN>;

print "The following questions are regarding parameters in Bowtie that help";
print "\nto optimize the settings in order to incorporate the maximum number";
print "\nof reads in the alignment.";

#system "./bowtie-build -f $fasta $ebwt";

#system "./bowtie -q -k 1 -n 3 -X 400 --best --solexa1.3-quals --sam $ebwt -1 $mate1 -2 $mate2 --al $output\.$species.al --un $output\.$species.un > $output\.$species\.sam";

#system "./samtools view -bS -o $output\.$species\.bam $output\.$species\.sam";

#system "./samtools sort $output\.$species\.bam $output\.$species\.sorted";

#system "./samtools pileup -cv -f $path2ebwt.fasta $output\.$species\.sorted.bam > $output\.$species\.SNPs.txt";

#printf ("\nTime elapsed for $output: %d\n\n",time-$start);
