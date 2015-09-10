#!/usr/bin/perl

########################################################################
# 
# 10/08/2009
#
# bowtiePipe0.2.pl
#
# Purpose: align mate paired-end short sequence reads using Bowtie 
# 
# Input: mate paired ends to bowtie
#
# Output: along with a bunch of files, a list of genes and exons and their read counts
#
# Syntax: perl bowtiePipe.pl Dmel Dsim ../Dmel/Dmel_genome.dm3 ../Dsec/Dsec_genome.droSec1 

########################################################################

use strict;
use warnings;

my$start = time;

print "\n\n";
print "################################################################\n";
print "#                                                              #\n";
print "# Welcome to Joe and Kraig's fabulous bioinformatics pipeline! #\n";
print "#                                                              #\n";
print "# This pipeline is meant to analyze high-throughput sequencing #\n";
print "# data from various species of Drosophila. It takes a set of   #\n";
print "# reads and aligns them to two Drosophila reference species of #\n";
print "# your choice. The two reference species should be those       #\n";
print "# that are used to make the hybrid cross of interest. This     #\n";
print "# pipeline can also analyze the parental strains.              #\n";
print "#                                                              #\n";
print "# The main purpose is to quantify expression values using      #\n";
print "# mRNA-Seq data. This is done by aligning the short reads to a #\n";
print "# reference genome to characterize the transcriptome and to    #\n";
print "# count the number of reads that align to particular alleles.  #\n";
print "#                                                              #\n";
print "# The pipeline is organized in the following way:              #\n";
print "#                                                              #\n";
print "# Step 1: perform alignment using Bowtie (Langmead et al.)     #\n";
print "# Step 2: combine mate pair information (if applicable)        #\n";
print "# Step 3: convert the .sam file(s) into .map and .bed files    #\n";
print "# Step 4: lift the .bed file from the worse-characterized      #\n";
print "#         species to the better-characterized species          #\n";
print "# Step 5: convert the genomic positions into genic locations   #\n";
print "#         (including exons) and removing any gaps between the  #\n";
print "#         two species                                          #\n";
print "# Step 6: classify the reads based on certain criteria and     #\n";
print "#         gather together the genes and exons for total counts #\n";
print "#                                                              #\n";
print "# You will now be prompted to answer a few questions in order  #\n";
print "# to tailor the alignment to your specific needs.              #\n";
print "#                                                              #\n";
print "#                       Good luck!                             #\n";
print "#                                                              #\n";
print "################################################################\n";
print "\n";

print "\nWhere is the first reference file? i.e. the one that all reads will be mapped back to? --> ";
my$fasta1 = <STDIN>;
my$ref1;
my$ebwt1;
if($fasta1 =~ /(\S+\/{1}\w+)\.fast[aq]$/){$ebwt1 = $1;} # takes the fasta file and finds the text for the index build
if($fasta1 =~ /\S+\/{1}(\w+)\.fa[staq]$/){$ref1 = $1;} 
#print "\nFasta: $fasta1\tRef: $ref1\tEBWT: $ebwt1\n\n";

print "\nWhere is the second reference file? --> ";
my$fasta2 = <STDIN>;
my$ref2;
my$ebwt2;
if($fasta2 =~ /(\S+\/{1}\w+)\.fast[aq]$/){$ebwt2 = $1;} # takes the fasta file and finds the text for the index build
if($fasta2 =~ /\S+\/{1}(\w+)\.fa[staq]$/){$ref2 = $1;} 
#print "\nFasta: $fasta2\tRef: $ref2\tEBWT: $ebwt2\n\n";

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

#Bowtie files
system "./bowtie -k 1 -v 0 --best --sam $ebwt1 $mate1 > $name1\.$ref1\.sam";
system "./bowtie -k 1 -v 0 --best --sam $ebwt1 $mate2 > $name2\.$ref1\.sam";
system "./bowtie -k 1 -v 0 --best --sam $ebwt2 $mate1 > $name1\.$ref2\.sam";
system "./bowtie -k 1 -v 0 --best --sam $ebwt2 $mate2 > $name2\.$ref2\.sam";

#Single-reference SNP calling
#system "./samtools view -bS -o $output\.$species\.bam $output\.$species\.sam";
#system "./samtools sort $output\.$species\.bam $output\.$species\.sorted";
#system "./samtools pileup -cv -f $path2ebwt.fasta $output\.$species\.sorted.bam > $output\.$species\.SNPs.txt";

$name = <>;
#reformatting step, convert the .sam files to .map and .bed files, combine mate files
system "perl reformatSAM2BED $name1\.$ref1\.sam $name2\.$ref1\.sam $name1\.$ref2\.sam $name2\.$ref2\.sam ";

my$chain = $ARGV[6]; 
#liftover the coordinates into Dmel space
system "./liftOver.MacOSX.ppc $secbed $chain $secbed\.lifted $secbed\.lifted.un";

printf ("\nTime elapsed: %d\n\n",time-$start);
