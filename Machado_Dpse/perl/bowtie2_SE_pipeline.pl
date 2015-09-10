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
$filename = "$pathebwt\.1\.bt2";
#print "\n$filename\n"; exit;
unless(-e $filename){system "bowtie2-build $ref $pathebwt";}
#print "\nbowtie-build $ref $pathebwt\n";

# mate
my$matlab;
if($mate =~ /^(\S+)\.(fastq|fq)$/)
{
  $matlab = $1;

  #gDNA is Phred+64

  system "bowtie2 -q --phred64 -p 12 --very-sensitive -x $pathebwt -U $mate -S $matlab\.$label\.bowtie2_very_sensitive\.sam --un $matlab\.$label\.bowtie2_very_sensitive\.un";

  system "bowtie2 -q --phred64 -p 12 --trim3 13 --very-sensitive -x $pathebwt -U $matlab\.$label\.bowtie2_very_sensitive\.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_13\.sam --un $matlab\.$label\.bowtie2_very_sensitive_trim3_13\.un";

  system "bowtie2 -q --phred64 -p 12 --trim3 26 --very-sensitive -x $pathebwt -U $matlab\.$label\.bowtie2_very_sensitive_trim3_13\.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_26\.sam --un $matlab\.$label\.bowtie2_very_sensitive_trim3_26\.un";

  system "bowtie2 -q --phred64 -p 12 --trim3 39 --very-sensitive -x $pathebwt -U $matlab\.$label\.bowtie2_very_sensitive_trim3_26\.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_39\.sam --un $matlab\.$label\.bowtie2_very_sensitive_trim3_39\.un";

  #RNA is Phred+33

  #system "bowtie2 -q --phred33 -p 2 --very-sensitive -x $pathebwt -U $mate -S $matlab\.$label\.bowtie2_very_sensitive\.sam --un $matlab\.$label\.bowtie2_very_sensitive\.un";

  #system "bowtie2 -q --phred33 -p 2 --trim3 13 --very-sensitive -x $pathebwt -U $matlab\.$label\.bowtie2_very_sensitive\.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_13\.sam --un $matlab\.$label\.bowtie2_very_sensitive_trim3_13\.un";

  #system "bowtie2 -q --phred33 -p 2 --trim3 26 --very-sensitive -x $pathebwt -U $matlab\.$label\.bowtie2_very_sensitive_trim3_13\.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_26\.sam --un $matlab\.$label\.bowtie2_very_sensitive_trim3_26\.un";

  #system "bowtie2 -q --phred33 -p 2 --trim3 39 --very-sensitive -x $pathebwt -U $matlab\.$label\.bowtie2_very_sensitive_trim3_26\.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_39\.sam --un $matlab\.$label\.bowtie2_very_sensitive_trim3_39\.un";
}
elsif($mate =~ /^(\S+)\.(fasta|fa)$/)
{
  $matlab = $1;

  system "bowtie2 -f -p 12 --very-sensitive -x $pathebwt -U $mate -S $matlab\.$label\.bowtie2_very_sensitive\.sam --un $matlab\.$label\.bowtie2_very_sensitive\.un";

  system "bowtie2 -f -p 12 --trim3 13 --very-sensitive -x $pathebwt -U $matlab\.$label\.bowtie2_very_sensitive\.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_13\.sam --un $matlab\.$label\.bowtie2_very_sensitive_trim3_13\.un";

  system "bowtie2 -f -p 12 --trim3 26 --very-sensitive -x $pathebwt -U $matlab\.$label\.bowtie2_very_sensitive_trim3_13\.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_26\.sam --un $matlab\.$label\.bowtie2_very_sensitive_trim3_26\.un";

  system "bowtie2 -f -p 12 --trim3 39 --very-sensitive -x $pathebwt -U $matlab\.$label\.bowtie2_very_sensitive_trim3_26\.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_39\.sam --un $matlab\.$label\.bowtie2_very_sensitive_trim3_39\.un";
}
