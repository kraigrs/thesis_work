#!/usr/bin/perl

########################################################################
# 
# 02/20/2012
#
# bowtiePipe.0.2.pl
#
# Purpose: use Bowtie to align reads to reference
# 
# Input: set of reads
#
# Output: alignments
#
# Syntax: perl bowtiePipe.0.2.pl <reference> <reads> 
#
########################################################################

use strict;
use warnings;

my$ref = $ARGV[0];
my$mate = $ARGV[1];

my$filename;

my$pathebwt;
if($ref =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$pathebwt = $1;}
#print "\n$pathebwt\n";

my$label;
if($ref =~ /^\S+\/(\S+)\.(fastq|fasta|fa|fq)$/){$label = $1;}
#print "\n$label\n";

# build references
$filename = "$pathebwt\.1\.ebwt";
#print "\n$filename\n"; exit;
unless(-e $filename){system "bowtie-build $ref $pathebwt";}
#print "\nbowtie-build $ref $pathebwt\n";

# mate
my$matlab;
if($mate =~ /^(\S+)\.(fastq|fq)$/)
{
  $matlab = $1;

  #system "bowtie -q -v 0 -m 1 -p 8 --best --sam $pathebwt $mate --al $matlab\.$label\.bowtie_v0_m1\.al --max $matlab\.$label\.bowtie_v0_m1\.max --un $matlab\.$label\.bowtie_v0_m1\.un > $matlab\.$label\.bowtie_v0_m1\.sam";
  #system "perl ../perl/samtoolsPipe.0.2.pl $ref $matlab\.$label\.bowtie_v0_m1\.sam";
   
  # Degner approach
  #system "bowtie -q -v 1 -m 1 -p 8 --best --sam $pathebwt $mate --al $matlab\.$label\.bowtie_v1_m1\.al --max $matlab\.$label\.bowtie_v1_m1\.max --un $matlab\.$label\.bowtie_v1_m1\.un > $matlab\.$label\.bowtie_v1_m1\.sam";
  #system "perl ../perl/samtoolsPipe.0.2.pl $ref $matlab\.$label\.bowtie_v1_m1\.sam";

  #system "bowtie -q -v 2 -m 1 -p 8 --best --sam $pathebwt $mate --al $matlab\.$label\.bowtie_v2_m1\.al --max $matlab\.$label\.bowtie_v2_m1\.max --un $matlab\.$label\.bowtie_v2_m1\.un > $matlab\.$label\.bowtie_v2_m1\.sam";
  #system "perl ../perl/samtoolsPipe.0.2.pl $ref $matlab\.$label\.bowtie_v2_m1\.sam";

  #system "bowtie -q -v 3 -m 1 -p 8 --best --sam $pathebwt $mate --al $matlab\.$label\.bowtie_v3_m1\.al --max $matlab\.$label\.bowtie_v3_m1\.max --un $matlab\.$label\.bowtie_v3_m1\.un > $matlab\.$label\.bowtie_v3_m1\.sam";
  #system "perl ../perl/samtoolsPipe.0.2.pl $ref $matlab\.$label\.bowtie_v3_m1\.sam";
}
elsif($mate =~ /^(\S+)\.(fasta|fa)$/)
{
  $matlab = $1;

  #system "bowtie -f -v 0 -m 1 -p 8 --best --sam $pathebwt $mate --al $matlab\.$label\.bowtie_v0_m1\.al --max $matlab\.$label\.bowtie_v0_m1\.max --un $matlab\.$label\.bowtie_v0_m1\.un > $matlab\.$label\.bowtie_v0_m1\.sam";
  #system "perl ../perl/samtoolsPipe.0.2.pl $ref $matlab\.$label\.bowtie_v0_m1\.sam";

  # Degner approach

  #system "bowtie -f -v 1 -m 1 -p 8 --best --sam $pathebwt $mate --al $matlab\.$label\.bowtie_v1_m1\.al --max $matlab\.$label\.bowtie_v1_m1\.max --un $matlab\.$label\.bowtie_v1_m1\.un > $matlab\.$label\.bowtie_v1_m1\.sam";
  #system "perl ../perl/samtoolsPipe.0.2.pl $ref $matlab\.$label\.bowtie_v1_m1\.sam";

  #system "bowtie -f -v 2 -m 1 -p 8 --best --sam $pathebwt $mate --al $matlab\.$label\.bowtie_v2_m1\.al --max $matlab\.$label\.bowtie_v2_m1\.max --un $matlab\.$label\.bowtie_v2_m1\.un > $matlab\.$label\.bowtie_v2_m1\.sam";
  #system "perl ../perl/samtoolsPipe.0.2.pl $ref $matlab\.$label\.bowtie_v2_m1\.sam";

  #system "bowtie -f -v 3 -m 1 -p 8 --best --sam $pathebwt $mate --al $matlab\.$label\.bowtie_v3_m1\.al --max $matlab\.$label\.bowtie_v3_m1\.max --un $matlab\.$label\.bowtie_v3_m1\.un > $matlab\.$label\.bowtie_v3_m1\.sam";
  #system "perl ../perl/samtoolsPipe.0.2.pl $ref $matlab\.$label\.bowtie_v3_m1\.sam";

  system "bowtie -f -n 3 -e 161 -l 36 -m 1 -p 8 --best --sam $pathebwt $mate --al $matlab\.$label\.bowtie_n3_e161_l36_m1\.al --max $matlab\.$label\.bowtie_n3_e161_l36_m1\.max --un $matlab\.$label\.bowtie_n3_e161_l36_m1\.un > $matlab\.$label\.bowtie_n3_e161_l36_m1\.sam";
  #system "perl ../perl/samtoolsPipe.0.2.pl $ref $matlab\.$label\.bowtie_n3_e161_l36_m1\.sam";
}
