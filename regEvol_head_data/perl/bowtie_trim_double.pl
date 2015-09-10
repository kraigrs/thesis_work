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
my$ref3 = $ARGV[2];
my$ref4 = $ARGV[3];
my$mate = $ARGV[4];

my$filename; my$pathebwt1; my$pathebwt2; my$label1; my$label2; my$matlab;
my$pathebwt3; my$pathebwt4; my$label3; my$label4;

if($ref1 =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$pathebwt1 = $1;}
if($ref2 =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$pathebwt2 = $1;}
if($ref3 =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$pathebwt3 = $1;}
if($ref4 =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$pathebwt4 = $1;}

if($ref1 =~ /^\S+\/(\S+)\.(fastq|fasta|fa|fq)$/){$label1 = $1;}
if($ref2 =~ /^\S+\/(\S+)\.(fastq|fasta|fa|fq)$/){$label2 = $1;}
if($ref3 =~ /^\S+\/(\S+)\.(fastq|fasta|fa|fq)$/){$label3 = $1;}
if($ref4 =~ /^\S+\/(\S+)\.(fastq|fasta|fa|fq)$/){$label4 = $1;}

# build references
$filename = "$pathebwt1\.1\.ebwt";
#print "\n$filename\n"; exit;
unless(-e $filename){system "bowtie-build $ref1 $pathebwt1";}
#print "\nbowtie-build $ref $pathebwt\n";

$filename = "$pathebwt2\.1\.ebwt";
#print "\n$filename\n"; exit;
unless(-e $filename){system "bowtie-build $ref2 $pathebwt2";}
#print "\nbowtie-build $ref $pathebwt\n";

$filename = "$pathebwt3\.1\.ebwt";
#print "\n$filename\n"; exit;
unless(-e $filename){system "bowtie-build $ref3 $pathebwt3";}
#print "\nbowtie-build $ref $pathebwt\n";

$filename = "$pathebwt4\.1\.ebwt";
#print "\n$filename\n"; exit;
unless(-e $filename){system "bowtie-build $ref4 $pathebwt4";}
#print "\nbowtie-build $ref $pathebwt\n";

if($mate =~ /^(\S+)\.(fastq|fq)$/)
{
  $matlab = $1;

  system "bowtie -q -p 1 -v 0 -m 1 --sam $pathebwt1 $mate --un $matlab\.$label1\.bowtie_v0_m1\.un > $matlab\.$label1\.bowtie_v0_m1\.sam";
  system "bowtie -q -p 1 -v 0 -m 1 --sam $pathebwt2 $matlab\.$label1\.bowtie_v0_m1\.un --un $matlab\.$label2\.bowtie_v0_m1\.un > $matlab\.$label2\.bowtie_v0_m1\.sam";
  system "bowtie -q -p 1 -v 0 -m 1 --sam $pathebwt3 $mate --un $matlab\.$label3\.bowtie_v0_m1\.un > $matlab\.$label3\.bowtie_v0_m1\.sam";
  system "bowtie -q -p 1 -v 0 -m 1 --sam $pathebwt4 $matlab\.$label3\.bowtie_v0_m1\.un --un $matlab\.$label4\.bowtie_v0_m1\.un > $matlab\.$label4\.bowtie_v0_m1\.sam";

  system "perl ../perl/get_trim_reads.pl $matlab\.$label2\.bowtie_v0_m1\.un $matlab\.$label4\.bowtie_v0_m1\.un $mate $matlab\.un";


  system "bowtie -q -p 1 -v 0 -m 1 --trim3 13 --sam $pathebwt1 $matlab\.un --un $matlab\.$label1\.bowtie_v0_m1_trim3_13\.un > $matlab\.$label1\.bowtie_v0_m1_trim3_13\.sam";
  system "bowtie -q -p 1 -v 0 -m 1 --trim3 13 --sam $pathebwt2 $matlab\.$label1\.bowtie_v0_m1_trim3_13\.un --un $matlab\.$label2\.bowtie_v0_m1_trim3_13\.un > $matlab\.$label2\.bowtie_v0_m1_trim3_13\.sam";

  system "bowtie -q -p 1 -v 0 -m 1 --trim3 13 --sam $pathebwt3 $matlab\.un --un $matlab\.$label3\.bowtie_v0_m1_trim3_13\.un > $matlab\.$label3\.bowtie_v0_m1_trim3_13\.sam";
  system "bowtie -q -p 1 -v 0 -m 1 --trim3 13 --sam $pathebwt4 $matlab\.$label3\.bowtie_v0_m1_trim3_13\.un --un $matlab\.$label4\.bowtie_v0_m1_trim3_13\.un > $matlab\.$label4\.bowtie_v0_m1_trim3_13\.sam";

  system "perl ../perl/get_trim_reads.pl $matlab\.$label2\.bowtie_v0_m1_trim3_13\.un $matlab\.$label4\.bowtie_v0_m1_trim3_13\.un $matlab\.un $matlab\_trim3_13\.un";


  system "bowtie -q -p 1 -v 0 -m 1 --trim3 26 --sam $pathebwt1 $matlab\_trim3_13\.un --un $matlab\.$label1\.bowtie_v0_m1_trim3_26\.un > $matlab\.$label1\.bowtie_v0_m1_trim3_26\.sam";
  system "bowtie -q -p 1 -v 0 -m 1 --trim3 26 --sam $pathebwt2 $matlab\.$label1\.bowtie_v0_m1_trim3_26\.un --un $matlab\.$label2\.bowtie_v0_m1_trim3_26\.un > $matlab\.$label2\.bowtie_v0_m1_trim3_26\.sam";

  system "bowtie -q -p 1 -v 0 -m 1 --trim3 26 --sam $pathebwt3 $matlab\_trim3_13\.un --un $matlab\.$label3\.bowtie_v0_m1_trim3_26\.un > $matlab\.$label3\.bowtie_v0_m1_trim3_26\.sam";
  system "bowtie -q -p 1 -v 0 -m 1 --trim3 26 --sam $pathebwt4 $matlab\.$label3\.bowtie_v0_m1_trim3_26\.un --un $matlab\.$label4\.bowtie_v0_m1_trim3_26\.un > $matlab\.$label4\.bowtie_v0_m1_trim3_26\.sam";

  system "perl ../perl/get_trim_reads.pl $matlab\.$label2\.bowtie_v0_m1_trim3_26\.un $matlab\.$label4\.bowtie_v0_m1_trim3_26\.un $matlab\_trim3_13\.un $matlab\_trim3_26\.un";


  system "bowtie -q -p 1 -v 0 -m 1 --trim3 39 --sam $pathebwt1 $matlab\_trim3_26\.un --un $matlab\.$label1\.bowtie_v0_m1_trim3_39\.un > $matlab\.$label1\.bowtie_v0_m1_trim3_39\.sam";
  system "bowtie -q -p 1 -v 0 -m 1 --trim3 39 --sam $pathebwt2 $matlab\.$label1\.bowtie_v0_m1_trim3_39\.un --un $matlab\.$label2\.bowtie_v0_m1_trim3_39\.un > $matlab\.$label2\.bowtie_v0_m1_trim3_39\.sam";

  system "bowtie -q -p 1 -v 0 -m 1 --trim3 39 --sam $pathebwt3 $matlab\_trim3_26\.un --un $matlab\.$label3\.bowtie_v0_m1_trim3_39\.un > $matlab\.$label3\.bowtie_v0_m1_trim3_39\.sam";
  system "bowtie -q -p 1 -v 0 -m 1 --trim3 39 --sam $pathebwt4 $matlab\.$label3\.bowtie_v0_m1_trim3_39\.un --un $matlab\.$label4\.bowtie_v0_m1_trim3_39\.un > $matlab\.$label4\.bowtie_v0_m1_trim3_39\.sam";

  system "perl ../perl/get_trim_reads.pl $matlab\.$label2\.bowtie_v0_m1_trim3_39\.un $matlab\.$label4\.bowtie_v0_m1_trim3_39\.un $matlab\_trim3_26\.un $matlab\_trim3_39\.un";


  system "bowtie -q -p 1 -v 0 -m 1 --trim3 52 --sam $pathebwt1 $matlab\_trim3_39\.un --un $matlab\.$label1\.bowtie_v0_m1_trim3_52\.un > $matlab\.$label1\.bowtie_v0_m1_trim3_52\.sam";
  system "bowtie -q -p 1 -v 0 -m 1 --trim3 52 --sam $pathebwt2 $matlab\.$label1\.bowtie_v0_m1_trim3_52\.un --un $matlab\.$label2\.bowtie_v0_m1_trim3_52\.un > $matlab\.$label2\.bowtie_v0_m1_trim3_52\.sam";

  system "bowtie -q -p 1 -v 0 -m 1 --trim3 52 --sam $pathebwt3 $matlab\_trim3_39\.un --un $matlab\.$label3\.bowtie_v0_m1_trim3_52\.un > $matlab\.$label3\.bowtie_v0_m1_trim3_52\.sam";
  system "bowtie -q -p 1 -v 0 -m 1 --trim3 52 --sam $pathebwt4 $matlab\.$label3\.bowtie_v0_m1_trim3_52\.un --un $matlab\.$label4\.bowtie_v0_m1_trim3_52\.un > $matlab\.$label4\.bowtie_v0_m1_trim3_52\.sam";

  system "perl ../perl/get_trim_reads.pl $matlab\.$label2\.bowtie_v0_m1_trim3_52\.un $matlab\.$label4\.bowtie_v0_m1_trim3_52\.un $matlab\_trim3_39\.un $matlab\_trim3_52\.un";


  system "bowtie -q -p 1 -v 0 -m 1 --trim3 65 --sam $pathebwt1 $matlab\_trim3_52\.un --un $matlab\.$label1\.bowtie_v0_m1_trim3_65\.un > $matlab\.$label1\.bowtie_v0_m1_trim3_65\.sam";
  system "bowtie -q -p 1 -v 0 -m 1 --trim3 65 --sam $pathebwt2 $matlab\.$label1\.bowtie_v0_m1_trim3_65\.un --un $matlab\.$label2\.bowtie_v0_m1_trim3_65\.un > $matlab\.$label2\.bowtie_v0_m1_trim3_65\.sam";

  system "bowtie -q -p 1 -v 0 -m 1 --trim3 65 --sam $pathebwt3 $matlab\_trim3_52\.un --un $matlab\.$label3\.bowtie_v0_m1_trim3_65\.un > $matlab\.$label3\.bowtie_v0_m1_trim3_65\.sam";
  system "bowtie -q -p 1 -v 0 -m 1 --trim3 65 --sam $pathebwt4 $matlab\.$label3\.bowtie_v0_m1_trim3_65\.un --un $matlab\.$label4\.bowtie_v0_m1_trim3_65\.un > $matlab\.$label4\.bowtie_v0_m1_trim3_65\.sam";
}
