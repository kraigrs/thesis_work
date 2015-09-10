#!/usr/bin/perl

########################################################################
# 
# 05/25/2010
#
# mosaikPipe.pl
#
# Purpose: align mate paired-end short sequence reads using MOSAIK 
# 
# Input: mate paired ends to MOSAIK
#
# Output: .sam alignment files (hopefully)
#
# Syntax: perl mosaikPipe.pl ../mel_sim_data/Resequencing/resequencing-assembly/zhr/zhr_reseq_may_all.fa ../mel_sim_data/Resequencing/resequencing-assembly/z30/z30_reseq_may_all.fa ../mel_sim_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.fastq ../mel_sim_data/mRNA-Seq/zhrXz30/zhrXz30.mate2.fastq 
#
########################################################################

use strict;
use warnings;

my$ref1 = $ARGV[0];
my$ref2 = $ARGV[1];
my$mate1 = $ARGV[2];
my$mate2 = $ARGV[3];

my$label1; my$refname1;
if($ref1 =~ /^(\S+\/(\S+))\.(fastq|fasta|fa|fq)$/){$label1 = $1; $refname1 = $2;}
#print "$label1\n$refname1\n";

my$label2; my$refname2;
if($ref2 =~ /^(\S+\/(\S+))\.(fastq|fasta|fa|fq)$/){$label2 = $1; $refname2 = $2;}
#print "$label2\n$refname2\n";

# build references if necessary
#system "MosaikBuild -fr $ref1 -oa $label1\.dat";
#system "MosaikBuild -fr $ref2 -oa $label2\.dat";

# make jump database to reduce memory footprint
#system "MosaikJump -ia $label1\.dat -out $label1\_15 -hs 15";
#system "MosaikJump -ia $label2\.dat -out $label2\_15 -hs 15";

# build reads (either fastq or fasta)
my$matlab1;
if($mate1 =~ /^(\S+)\.(fastq|fq)$/)
{
  $matlab1 = $1;
  #system "MosaikBuild -q $mate1 -out $matlab1\.dat -st illumina";
}
elsif($mate1 =~ /^(\S+)\.(fasta|fa)$/)
{
  $matlab1 = $1;
  #system "MosaikBuild -fr $mate1 -assignQual 40 -out $matlab1\.dat -st illumina";
}

my$matlab2;
if($mate2 =~ /^(\S+)\.(fastq|fq)$/)
{
  $matlab2 = $1;
  #system "MosaikBuild -q $mate2 -out $matlab2\.dat -st illumina";
}
elsif($mate2 =~ /^(\S+)\.(fasta|fa)$/)
{
  $matlab2 = $1;
  #system "MosaikBuild -fr $mate2 -assignQual 40 -out $matlab2\.dat -st illumina";
}

################################
# align reads to the reference #
################################

# mate 1
system "MosaikAligner -in $matlab1\.dat -out $matlab1\.$refname1\.dat -ia $label1\.dat -hs 15 -mm 0 -m unique -j $label1\_15 -p 8 > $matlab1\.$refname1\.output";
system "MosaikAligner -in $matlab1\.dat -out $matlab1\.$refname2\.dat -ia $label2\.dat -hs 15 -mm 0 -m unique -j $label2\_15 -p 8 > $matlab1\.$refname2\.output";

# mate 2
system "MosaikAligner -in $matlab2\.dat -out $matlab2\.$refname1\.dat -ia $label1\.dat -hs 15 -mm 0 -m unique -j $label1\_15 -p 8 > $matlab2\.$refname1\.output";
system "MosaikAligner -in $matlab2\.dat -out $matlab2\.$refname2\.dat -ia $label2\.dat -hs 15 -mm 0 -m unique -j $label2\_15 -p 8 > $matlab2\.$refname2\.output";

# convert files to sam format (remove non-unique reads with -u)

system "MosaikText -in $matlab1\.$label1\.dat -u -sam $matlab1\.$label1\.mosaik\.sam";
system "MosaikText -in $matlab1\.$label2\.dat -u -sam $matlab1\.$label2\.mosaik\.sam";

system "MosaikText -in $matlab2\.$label1\.dat -u -sam $matlab2\.$label1\.mosaik\.sam";
system "MosaikText -in $matlab2\.$label2\.dat -u -sam $matlab2\.$label2\.mosaik\.sam";

