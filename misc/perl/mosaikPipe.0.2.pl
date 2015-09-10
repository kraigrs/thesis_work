#!/usr/bin/perl

########################################################################
# 
# 08/31/2010
#
# mosaikPipe.0.2.pl
#
# Purpose: align short sequence reads using MOSAIK 
# 
# Input: reads
#
# Output: .sam alignment files (hopefully)
#
# Syntax: perl mosaikPipe.0.2.pl ../mel_sim_data/Resequencing/resequencing-assembly/zhr/zhr_reseq_may_all.fa ../mel_sim_data/mRNA-Seq/zhrXz30/zhrXz30.mate1.fastq  
#
########################################################################

use strict;
use warnings;

my$ref = $ARGV[0];
my$mate = $ARGV[1];

our($label,$refname,$filename,$matlab);

if($ref =~ /^(\S+\/(\S+))\.(fastq|fasta|fa|fq)$/){$label = $1; $refname = $2;}
#print "label:$label\trefname:$refname\n"; exit;

# build references if necessary
$filename = "$label\.dat";
unless(-e $filename){system "MosaikBuild -fr $ref -oa $label\.dat";}
#print "MosaikBuild -fr $ref -oa $label\.dat";

# make jump database to reduce memory footprint
$filename = "$label\_14_meta.jmp";
unless(-e $filename){system "MosaikJump -ia $label\.dat -out $label\_14 -hs 14 -iupac";}
#print "MosaikJump -ia $label\.dat -out $label\_15 -hs 15";

# build reads (either fastq or fasta)
if($mate =~ /^(\S+)\.(fastq|fq)$/)
{
  $matlab = $1;
  $filename = "$matlab\.dat";
  unless(-e $filename){system "MosaikBuild -q $mate -out $matlab\.dat -st illumina";}
}
elsif($mate =~ /^(\S+)\.(fasta|fa)$/)
{
  $matlab = $1;
  $filename = "$matlab\.dat";
  unless(-e $filename){system "MosaikBuild -fr $mate -out $matlab\.dat -st illumina -assignQual 60";}
  #print "MosaikBuild -fr $mate -out $matlab\.dat -st illumina -assignQual 60";
}

################################
# align reads to the reference #
################################


#print "\n\nMosaikAligner -in $matlab\.dat -out $matlab\.$refname\.dat -ia $label\.dat -hs 14 -mm 0 -m unique -j $label\_14 -p 2 -annse 2.1.26.se.100.005.ann -annpe 2.1.26.pe.100.0065.ann > $matlab\.$refname\.output\n\n";
#exit;
system "MosaikAligner -in $matlab\.dat -out $matlab\.$refname\.dat -ia $label\.dat -hs 14 -mm 0 -m unique -j $label\_14 -p 2 -annse 2.1.26.se.100.005.ann -annpe 2.1.26.pe.100.0065.ann > $matlab\.$refname\.output";

#print "MosaikAligner -in $matlab\.dat -out $matlab\.$refname\.dat -ia $label\.dat -hs 15 -mm 0 -m unique -j $label\_15 -p 2 > $matlab\.$refname\.output";
#system "MosaikAligner -in $matlab\.dat -out $matlab\.$refname\.dat -ia $label\.dat -hs 15 -mm 0 -m unique -j $label\_15 -p 2 > $matlab\.$refname\.output";

#system "MosaikAligner -in $matlab\.dat -out $matlab\.$refname\.dat -ia $label\.dat -hs 15 -mm 0 -m unique -j $label\_15 -p 2";

# convert files to .bed format

#system "MosaikText -in $matlab\.$refname\.dat -bed $matlab\.$refname\.mosaik\.bed";



