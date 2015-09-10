#!/usr/bin/perl

########################################################################
# 
# 03/30/2012
#
# samtoolsPipe.0.2.pl
#
# Purpose: run samtools on sam files
# 
# Input: set of reads
#
# Output: alignments
#
# Syntax: perl samtoolsPipe.0.2.pl <reference> <SAM> 
#
########################################################################

use strict;
use warnings;

my$ref = $ARGV[0];
my$SAM = $ARGV[1];

my$filename;

my$pathebwt;
if($ref =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$pathebwt = $1;}
#print "\n$pathebwt\n";

my$label;
if($ref =~ /^\S+\/(\S+)\.(fastq|fasta|fa|fq)$/){$label = $1;}
#print "\n$label\n";

my$prefix;
if($SAM =~ /^(\S+)\.sam$/){$prefix = $1;}
#print "\n$prefix\n";

# convert SAM to BAM
$filename = "$prefix\.bam";
#print "\n$filename\n"; exit;
unless(-e $filename){system "samtools view -S -b -T $ref $SAM > $prefix\.bam";}
print "\nsamtools view -S -b -T $ref $SAM > $prefix\.bam\n";

# sort BAM
$filename = "$prefix\.sorted.bam";
#print "\n$filename\n"; exit;
unless(-e $filename){system "samtools sort $prefix\.bam $prefix\.sorted";}
print "\nsamtools sort $prefix\.bam $prefix\.sorted\n";

# index BAM
#$filename = "$prefix\.sorted.bam.bai";
#print "\n$filename\n"; exit;
#unless(-e $filename){system "samtools index $prefix\.sorted.bam";}
#print "\nsamtools index $prefix\.sorted.bam\n";

# create pileup
$filename = "$prefix\.pileup.txt";
#print "\n$filename\n"; exit;
unless(-e $filename){system "samtools mpileup -f $ref $prefix\.sorted.bam > $prefix\.pileup\.txt";}
print "\nsamtools mpileup -f $ref $prefix\.sorted.bam > $prefix\.pileup\.txt\n";
