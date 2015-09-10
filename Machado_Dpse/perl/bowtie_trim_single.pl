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

my$ref = $ARGV[0];
my$mate = $ARGV[1];

my$filename; my$pathebwt; my$label; my$matlab;

if($ref =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$pathebwt = $1;}
#print "\n$pathebwt\n";

if($ref =~ /^\S+\/(\S+)\.(fastq|fasta|fa|fq)$/){$label = $1;}
#print "\n$label\n";

# build references
$filename = "$pathebwt\.1\.ebwt";
#print "\n$filename\n"; exit;
unless(-e $filename){system "bowtie-build $ref $pathebwt";}
#print "\nbowtie-build $ref $pathebwt\n";

if($mate =~ /^(\S+)\.(fastq|fq)$/)
{
  $matlab = $1;

  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --best --sam $pathebwt $mate --al $matlab\.$label\.bowtie_v0_m1\.al --max $matlab\.$label\.bowtie_v0_m1\.max --un $matlab\.$label\.bowtie_v0_m1\.un > $matlab\.$label\.bowtie_v0_m1\.sam";

  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --trim3 24 --best --sam $pathebwt $matlab\.$label\.bowtie_v0_m1\.un --al $matlab\.$label\.bowtie_v0_m1_76\.al --max $matlab\.$label\.bowtie_v0_m1_76\.max --un $matlab\.$label\.bowtie_v0_m1_76\.un > $matlab\.$label\.bowtie_v0_m1_76\.sam";

  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --trim3 37 --best --sam $pathebwt $matlab\.$label\.bowtie_v0_m1_76\.un --al $matlab\.$label\.bowtie_v0_m1_63\.al --max $matlab\.$label\.bowtie_v0_m1_63\.max --un $matlab\.$label\.bowtie_v0_m1_63\.un > $matlab\.$label\.bowtie_v0_m1_63\.sam";

  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --trim3 50 --best --sam $pathebwt $matlab\.$label\.bowtie_v0_m1_63\.un --al $matlab\.$label\.bowtie_v0_m1_50\.al --max $matlab\.$label\.bowtie_v0_m1_50\.max --un $matlab\.$label\.bowtie_v0_m1_50\.un > $matlab\.$label\.bowtie_v0_m1_50\.sam";

  system "bowtie -q -v 0 -m 1 -p 2 --trim5 1 --trim3 63 --best --sam $pathebwt $matlab\.$label\.bowtie_v0_m1_50\.un --al $matlab\.$label\.bowtie_v0_m1_37\.al --max $matlab\.$label\.bowtie_v0_m1_37\.max --un $matlab\.$label\.bowtie_v0_m1_37\.un > $matlab\.$label\.bowtie_v0_m1_37\.sam";
}
