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

my$pathebwt;
if($ref =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$pathebwt = $1;}
print "$pathebwt\n";

my$label;
if($ref =~ /^\S+\/(\S+)\.(fastq|fasta|fa|fq)$/){$label = $1;}
print "$label\n";

# build references
system "bowtie-build $ref $pathebwt";

# mate
my$matlab;
if($mate =~ /^(\S+)\.(fastq|fq)$/)
{
  $matlab = $1;
  system "bowtie -q -v 0 -m 1 --sam $pathebwt $mate --al $matlab\.$label\.bowtie\.sam --max $matlab\.$label\.bowtie\.multiple.sam --un $matlab\.$label\.bowtie\.un.sam";
}
elsif($mate1 =~ /^(\S+)\.(fasta|fa)$/)
{
  $matlab = $1;
  system "bowtie -f -v 0 -m 1 --sam $pathebwt $mate --al $matlab\.$label\.bowtie\.sam --max $matlab\.$label\.bowtie\.multiple.sam --un $matlab\.$label\.bowtie\.un.sam";
}
