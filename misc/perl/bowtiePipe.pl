#!/usr/bin/perl

########################################################################
# 
# 05/25/2010
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
# Syntax: perl bowtiePipe.pl ../Dmel/dm3.fa ../Dsec/droSim1.fa ../mel_sim_data/mRNA-Seq/Dmel-zhr_cDNA/Dmel.mate1.fa ../mel_sim_data/mRNA-Seq/Dmel-zhr_cDNA/Dmel.mate2.fa 
#
########################################################################

use strict;
use warnings;

my$ref1 = $ARGV[0];
my$ref2 = $ARGV[1];
my$mate1 = $ARGV[2];
my$mate2 = $ARGV[3];

my$pathebwt1;
if($ref1 =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$pathebwt1 = $1;}
print "$pathebwt1\n";
my$pathebwt2;
if($ref2 =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$pathebwt2 = $1;}
print "$pathebwt2\n";

my$label1;
if($ref1 =~ /^\S+\/(\S+)\.(fastq|fasta|fa|fq)$/){$label1 = $1;}
print "$label1\n";
my$label2;
if($ref2 =~ /^\S+\/(\S+)\.(fastq|fasta|fa|fq)$/){$label2 = $1;}
print "$label2\n";

# build references
#system "bowtie-build $ref1 $pathebwt1";
#system "bowtie-build $ref2 $pathebwt2";

# mate 1
my$matlab1;
if($mate1 =~ /^(\S+)\.(fastq|fq)$/)
{
  $matlab1 = $1;
  system "bowtie -q -k 1 -v 0 --best --sam $pathebwt1 $mate1 --un $matlab1\.$label1\.bowtie\.sam.un > $matlab1\.$label1\.bowtie\.sam";
  system "bowtie -q -k 1 -v 0 --best --sam $pathebwt2 $mate1 --un $matlab1\.$label2\.bowtie\.sam.un > $matlab1\.$label2\.bowtie\.sam";
}
elsif($mate1 =~ /^(\S+)\.(fasta|fa)$/)
{
  $matlab1 = $1;
  system "bowtie -f -k 1 -v 0 --best --sam $pathebwt1 $mate1 --un $matlab1\.$label1\.bowtie\.sam.un > $matlab1\.$label1\.bowtie\.sam";
  system "bowtie -f -k 1 -v 0 --best --sam $pathebwt2 $mate1 --un $matlab1\.$label2\.bowtie\.sam.un > $matlab1\.$label2\.bowtie\.sam";
}

# mate 2
my$matlab2;
if($mate2 =~ /^(\S+)\.(fastq|fq)$/)
{
  $matlab2 = $1;
  system "bowtie -q -k 1 -v 0 --best --sam $pathebwt1 $mate2 --un $matlab2\.$label1\.bowtie\.sam.un > $matlab2\.$label1\.bowtie\.sam";
  system "bowtie -q -k 1 -v 0 --best --sam $pathebwt2 $mate2 --un $matlab2\.$label2\.bowtie\.sam.un > $matlab2\.$label2\.bowtie\.sam";
}
elsif($mate2 =~ /^(\S+)\.(fasta|fa)$/)
{
  $matlab2 = $1;
  system "bowtie -f -k 1 -v 0 --best --sam $pathebwt1 $mate2 --un $matlab2\.$label1\.bowtie\.sam.un > $matlab2\.$label1\.bowtie\.sam";
  system "bowtie -f -k 1 -v 0 --best --sam $pathebwt2 $mate2 --un $matlab2\.$label2\.bowtie\.sam.un > $matlab2\.$label2\.bowtie\.sam";
}
