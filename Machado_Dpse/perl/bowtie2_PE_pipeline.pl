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
my$mate1 = $ARGV[1];
my$mate2 = $ARGV[2];

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
if($mate1 =~ /^(\S+)\_R[12]{1}\.(fastq|fq)$/)
{
  $matlab = $1;

  system "bowtie2 -q --phred64 -p 12 --very-sensitive --no-mixed --no-discordant --no-dovetail --no-contain --no-overlap -x $pathebwt -1 $mate1 -2 $mate2 -S $matlab\.$label\.bowtie2_very_sensitive\.sam --un-conc $matlab\.$label\.bowtie2_very_sensitive.un";

  system "bowtie2 -q --phred64 -p 12 --trim3 13 --very-sensitive --no-mixed --no-discordant --no-dovetail --no-contain --no-overlap -x $pathebwt -1 $matlab\.$label\.bowtie2_very_sensitive.1.un -2 $matlab\.$label\.bowtie2_very_sensitive.2.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_13\.sam --un-conc $matlab\.$label\.bowtie2_very_sensitive_trim3_13.un";

  system "bowtie2 -q --phred64 -p 12 --trim3 26 --very-sensitive --no-mixed --no-discordant --no-dovetail --no-contain --no-overlap -x $pathebwt -1 $matlab\.$label\.bowtie2_very_sensitive_trim3_13.1.un -2 $matlab\.$label\.bowtie2_very_sensitive_trim3_13.2.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_26\.sam --un-conc $matlab\.$label\.bowtie2_very_sensitive_trim3_26.un";

  system "bowtie2 -q --phred64 -p 12 --trim3 39 --very-sensitive --no-mixed --no-discordant --no-dovetail --no-contain --no-overlap -x $pathebwt -1 $matlab\.$label\.bowtie2_very_sensitive_trim3_26.1.un -2 $matlab\.$label\.bowtie2_very_sensitive_trim3_26.2.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_39\.sam --un-conc $matlab\.$label\.bowtie2_very_sensitive_trim3_39.un";
}
elsif($mate1 =~/^(\S+)\_R[12]{1}\.(fasta|fa)$/)
{
  $matlab = $1;

  system "bowtie2 -f -p 12 --very-sensitive --no-mixed --no-discordant --no-dovetail --no-contain --no-overlap -x $pathebwt -1 $mate1 -2 $mate2 -S $matlab\.$label\.bowtie2_very_sensitive\.sam --un-conc $matlab\.$label\.bowtie2_very_sensitive.un";

  system "bowtie2 -f -p 12 --trim3 13 --very-sensitive --no-mixed --no-discordant --no-dovetail --no-contain --no-overlap -x $pathebwt -1 $matlab\.$label\.bowtie2_very_sensitive.1.un -2 $matlab\.$label\.bowtie2_very_sensitive.2.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_13\.sam --un-conc $matlab\.$label\.bowtie2_very_sensitive_trim3_13.un";

  system "bowtie2 -f -p 12 --trim3 26 --very-sensitive --no-mixed --no-discordant --no-dovetail --no-contain --no-overlap -x $pathebwt -1 $matlab\.$label\.bowtie2_very_sensitive_trim3_13.1.un -2 $matlab\.$label\.bowtie2_very_sensitive_trim3_13.2.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_26\.sam --un-conc $matlab\.$label\.bowtie2_very_sensitive_trim3_26.un";

  system "bowtie2 -f -p 12 --trim3 39 --very-sensitive --no-mixed --no-discordant --no-dovetail --no-contain --no-overlap -x $pathebwt -1 $matlab\.$label\.bowtie2_very_sensitive_trim3_26.1.un -2 $matlab\.$label\.bowtie2_very_sensitive_trim3_26.2.un -S $matlab\.$label\.bowtie2_very_sensitive_trim3_39\.sam --un-conc $matlab\.$label\.bowtie2_very_sensitive_trim3_39.un";
}
