#!/usr/bin/perl

########################################################################
# 
# 05/23/2013
#
# bowtie_iterate.pl
#
# Purpose: 
# 
# Input: 
#
# Output: 
#
# Syntax: perl bowtie_iterate.pl <reference> <reads> 
#
########################################################################

use strict;
use warnings;

my$ref1 = $ARGV[0];
my$ref2 = $ARGV[1];
my$mate = $ARGV[2];

my$filename; my$pathebwt1; my$pathebwt2; my$label1; my$label2; my$matlab;

if($ref1 =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$pathebwt1 = $1;}
if($ref2 =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$pathebwt2 = $1;}

if($ref1 =~ /^\S+\/(\S+)\.(fastq|fasta|fa|fq)$/){$label1 = $1;}
if($ref2 =~ /^\S+\/(\S+)\.(fastq|fasta|fa|fq)$/){$label2 = $1;}

# build references
$filename = "$pathebwt1\.1\.ebwt";
#print "\n$filename\n"; exit;
unless(-e $filename){system "bowtie-build $ref1 $pathebwt1";}
#print "\nbowtie-build $ref $pathebwt\n";

$filename = "$pathebwt2\.1\.ebwt";
#print "\n$filename\n"; exit;
unless(-e $filename){system "bowtie-build $ref2 $pathebwt2";}
#print "\nbowtie-build $ref $pathebwt\n";

if($mate =~ /^(\S+)\.(fastq|fq)$/)
{
  $matlab = $1;

  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --best --sam $pathebwt1 $mate --al $matlab\.$label1\.bowtie_v0_m1\.al --max $matlab\.$label1\.bowtie_v0_m1\.max --un $matlab\.$label1\.bowtie_v0_m1\.un > $matlab\.$label1\.bowtie_v0_m1\.sam";
  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --best --sam $pathebwt2 $matlab\.$label1\.bowtie_v0_m1\.un --al $matlab\.$label2\.bowtie_v0_m1\.al --max $matlab\.$label2\.bowtie_v0_m1\.max --un $matlab\.$label2\.bowtie_v0_m1\.un > $matlab\.$label2\.bowtie_v0_m1\.sam";

  #sort $matlab\.$label1\.bowtie_v0_m1\.un > $matlab\.$label1\.bowtie_v0_m1\.un.sorted
  #sort $matlab\.$label2\.bowtie_v0_m1\.un > $matlab\.$label2\.bowtie_v0_m1\.un.sorted
  #comm -12 $matlab\.$label1\.bowtie_v0_m1\.un.sorted $matlab\.$label2\.bowtie_v0_m1\.un.sorted | sed -r 's/^\t//'

  # this comm pulls out reads that didnt align to either genome and should be trimmed and added to the next step

  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --trim3 24 --best --sam $pathebwt1 $matlab\.$label2\.bowtie_v0_m1\.un --al $matlab\.$label1\.bowtie_v0_m1_76\.al --max $matlab\.$label1\.bowtie_v0_m1_76\.max --un $matlab\.$label1\.bowtie_v0_m1_76\.un > $matlab\.$label1\.bowtie_v0_m1_76\.sam";
  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --trim3 24 --best --sam $pathebwt2 $matlab\.$label1\.bowtie_v0_m1_76\.un --al $matlab\.$label2\.bowtie_v0_m1_76\.al --max $matlab\.$label2\.bowtie_v0_m1_76\.max --un $matlab\.$label2\.bowtie_v0_m1_76\.un > $matlab\.$label2\.bowtie_v0_m1_76\.sam";

  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --trim3 37 --best --sam $pathebwt1 $matlab\.$label2\.bowtie_v0_m1_76\.un --al $matlab\.$label1\.bowtie_v0_m1_63\.al --max $matlab\.$label1\.bowtie_v0_m1_63\.max --un $matlab\.$label1\.bowtie_v0_m1_63\.un > $matlab\.$label1\.bowtie_v0_m1_63\.sam";
  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --trim3 37 --best --sam $pathebwt2 $matlab\.$label1\.bowtie_v0_m1_63\.un --al $matlab\.$label2\.bowtie_v0_m1_63\.al --max $matlab\.$label2\.bowtie_v0_m1_63\.max --un $matlab\.$label2\.bowtie_v0_m1_63\.un > $matlab\.$label2\.bowtie_v0_m1_63\.sam";

  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --trim3 50 --best --sam $pathebwt1 $matlab\.$label2\.bowtie_v0_m1_63\.un --al $matlab\.$label1\.bowtie_v0_m1_50\.al --max $matlab\.$label1\.bowtie_v0_m1_50\.max --un $matlab\.$label1\.bowtie_v0_m1_50\.un > $matlab\.$label1\.bowtie_v0_m1_50\.sam";
  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --trim3 50 --best --sam $pathebwt2 $matlab\.$label1\.bowtie_v0_m1_50\.un --al $matlab\.$label2\.bowtie_v0_m1_50\.al --max $matlab\.$label2\.bowtie_v0_m1_50\.max --un $matlab\.$label2\.bowtie_v0_m1_50\.un > $matlab\.$label2\.bowtie_v0_m1_50\.sam";

  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --trim3 63 --best --sam $pathebwt1 $matlab\.$label2\.bowtie_v0_m1_50\.un --al $matlab\.$label1\.bowtie_v0_m1_37\.al --max $matlab\.$label1\.bowtie_v0_m1_37\.max --un $matlab\.$label1\.bowtie_v0_m1_37\.un > $matlab\.$label1\.bowtie_v0_m1_37\.sam";
  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --trim3 63 --best --sam $pathebwt2 $matlab\.$label1\.bowtie_v0_m1_37\.un --al $matlab\.$label2\.bowtie_v0_m1_37\.al --max $matlab\.$label2\.bowtie_v0_m1_37\.max --un $matlab\.$label2\.bowtie_v0_m1_37\.un > $matlab\.$label2\.bowtie_v0_m1_37\.sam";
}
