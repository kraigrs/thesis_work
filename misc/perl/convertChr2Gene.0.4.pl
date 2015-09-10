#!/usr/bin/perl

###################################################################################
#
# 04/28/2010
#
# convertChr2Gene.0.4.pl
#
# Purpose: query a gap file to see if mapped read falls in a gap and, if not, 
#          in which gene and exon it happens to fall in 
# 
# Input: .map file, list of gaps, list of constitutive exons
#
# Output: for every read, an indication of being in a gap, constitutive exon, or neither 
#
# Usage: perl convertChr2Gene.0.4.pl <bedfile> <gaps> <constExons> <output>  
#        <bedfile>    ==> the .bed file to be converted 
#        <gaps>       ==> a file containing the gap regions
#        <constExons> ==> a file containing a list of constitutive exons
#        <output>     ==> the output file to write converted reads to
#
# Improvement: This uses a perl module called Tree::R which allows logarithmic searching for gaps and exons
#
###################################################################################

use strict;
use warnings;
use Tree::R;

main ();
sub main 
{ 
  my@elements;
  my@bbox;
  my@args;
  my@results;
  my@rect;
  my%exons;
  my%gaps;
  
  my$bedfile = $ARGV[0];  
  my$gaps_file = $ARGV[1];
  my$constEx_file = $ARGV[2]; 
  my$read_out = $ARGV[3];

  my$object;
  my$line; 
  my$begin; 
  my$end;
  my$chrom; 
  my$start;
  my$stop;  
  my$read;
  my$gene;

  my$gap_chr2L = new Tree::R;
  my$gap_chr2LHet = new Tree::R;
  my$gap_chr2R = new Tree::R;
  my$gap_chr2RHet = new Tree::R;
  my$gap_chr3L = new Tree::R;
  my$gap_chr3LHet = new Tree::R;
  my$gap_chr3R = new Tree::R;
  my$gap_chr3RHet = new Tree::R;
  my$gap_chr4 = new Tree::R;
  my$gap_chrM = new Tree::R;
  my$gap_chrU = new Tree::R;
  my$gap_chrUextra = new Tree::R;
  my$gap_chrX = new Tree::R;
  my$gap_chrXHet = new Tree::R;
  my$gap_chrYHet = new Tree::R;

  my$exon_chr2L = new Tree::R;
  my$exon_chr2LHet = new Tree::R;
  my$exon_chr2R = new Tree::R;
  my$exon_chr2RHet = new Tree::R;
  my$exon_chr3L = new Tree::R;
  my$exon_chr3LHet = new Tree::R;
  my$exon_chr3R = new Tree::R;
  my$exon_chr3RHet = new Tree::R;
  my$exon_chr4 = new Tree::R;
  my$exon_chrM = new Tree::R;
  my$exon_chrU = new Tree::R;
  my$exon_chrUextra = new Tree::R;
  my$exon_chrX = new Tree::R;
  my$exon_chrXHet = new Tree::R;
  my$exon_chrYHet = new Tree::R;

  sub build_gap_tree
  {
    #arguments: chrom,begin,end
    @args = shift(@_);
    $chrom = $args[0];
    $begin = $args[1];
    $end = $args[2];
    
    @bbox = ($begin,0,$end,0);                #force box to be a line (20,0,30,0)
    $object = $chrom."_".$begin."_".$end;     #looks like "chrX_123_456"

    if($chrom eq "chr2L"){$gap_chr2L->insert($object,@bbox);}
    elsif($chrom eq "chr2LHet"){$gap_chr2LHet->insert($object,@bbox);}
    elsif($chrom eq "chr2R"){$gap_chr2R->insert($object,@bbox);}
    elsif($chrom eq "chr2RHet"){$gap_chr2RHet->insert($object,@bbox);}
    elsif($chrom eq "chr3L"){$gap_chr3L->insert($object,@bbox);}
    elsif($chrom eq "chr3LHet"){$gap_chr3LHet->insert($object,@bbox);}
    elsif($chrom eq "chr3R"){$gap_chr3R->insert($object,@bbox);}
    elsif($chrom eq "chr3RHet"){$gap_chr3RHet->insert($object,@bbox);}
    elsif($chrom eq "chr4"){$gap_chr4->insert($object,@bbox);}
    elsif($chrom eq "chrM"){$gap_chrM->insert($object,@bbox);}
    elsif($chrom eq "chrU"){$gap_chrU->insert($object,@bbox);}
    elsif($chrom eq "chrUextra"){$gap_chrUextra->insert($object,@bbox);}
    elsif($chrom eq "chrX"){$gap_chrX->insert($object,@bbox);}
    elsif($chrom eq "chrXHet"){$gap_chrXHet->insert($object,@bbox);}
    elsif($chrom eq "chrYHet"){$gap_chrYHet->insert($object,@bbox);}
    else{print "Gap on $chrom from $begin to $end does not exist!";}
  }

  sub build_exon_tree
  {
    #arguments: chrom,begin,end,gene
    @args = shift(@_);
    $chrom = $args[0];
    $begin = $args[1];
    $end = $args[2];
    $gene = $args[3];
    
    @bbox = ($begin,0,$end,0);                #force box to be a line (20,0,30,0)
    $object = $gene."_".$begin."_".$end;      #looks like "CG9876_123_456"

    if($chrom eq "chr2L"){$exon_chr2L->insert($object,@bbox);}
    elsif($chrom eq "chr2LHet"){$exon_chr2LHet->insert($object,@bbox);}
    elsif($chrom eq "chr2R"){$exon_chr2R->insert($object,@bbox);}
    elsif($chrom eq "chr2RHet"){$exon_chr2RHet->insert($object,@bbox);}
    elsif($chrom eq "chr3L"){$exon_chr3L->insert($object,@bbox);}
    elsif($chrom eq "chr3LHet"){$exon_chr3LHet->insert($object,@bbox);}
    elsif($chrom eq "chr3R"){$exon_chr3R->insert($object,@bbox);}
    elsif($chrom eq "chr3RHet"){$exon_chr3RHet->insert($object,@bbox);}
    elsif($chrom eq "chr4"){$exon_chr4->insert($object,@bbox);}
    elsif($chrom eq "chrM"){$exon_chrM->insert($object,@bbox);}
    elsif($chrom eq "chrU"){$exon_chrU->insert($object,@bbox);}
    elsif($chrom eq "chrUextra"){$exon_chrUextra->insert($object,@bbox);}
    elsif($chrom eq "chrX"){$exon_chrX->insert($object,@bbox);}
    elsif($chrom eq "chrXHet"){$exon_chrXHet->insert($object,@bbox);}
    elsif($chrom eq "chrYHet"){$exon_chrYHet->insert($object,@bbox);}
    else{print "Exon on $chrom from $begin to $end does not exist!";}
  }

  sub query_gaps
  {
    #arguments: chrom,begin,end
    @args = shift(@_);
    $chrom = $args[0];
    $begin = $args[1];
    $end = $args[2];
    
    @rect = ($begin,0,$end,0);                #force box to be a line (20,0,30,0)

    if($chrom eq "chr2L"){$gap_chr2L->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr2LHet"){$gap_chr2LHet->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr2R"){$gap_chr2R->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr2RHet"){$gap_chr2RHet->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr3L"){$gap_chr3L->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr3LHet"){$gap_chr3LHet->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr3R"){$gap_chr3R->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr3RHet"){$gap_chr3RHet->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr4"){$gap_chr4->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chrM"){$gap_chrM->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chrU"){$gap_chrU->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chrUextra"){$gap_chrUextra->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chrX"){$gap_chrX->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chrXHet"){$gap_chrXHet->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chrYHet"){$gap_chrYHet->query_partly_within_rect(@rect,\@results);}
    return @results;
  }

  sub query_exons
  {
    #arguments: chrom,begin,end
    @args = shift(@_);
    $chrom = $args[0];
    $begin = $args[1];
    $end = $args[2];
    
    @rect = ($begin,0,$end,0);                #force box to be a line (20,0,30,0)

    if($chrom eq "chr2L"){$exon_chr2L->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr2LHet"){$exon_chr2LHet->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr2R"){$exon_chr2R->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr2RHet"){$exon_chr2RHet->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr3L"){$exon_chr3L->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr3LHet"){$exon_chr3LHet->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr3R"){$exon_chr3R->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr3RHet"){$exon_chr3RHet->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chr4"){$exon_chr4->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chrM"){$exon_chrM->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chrU"){$exon_chrU->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chrUextra"){$exon_chrUextra->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chrX"){$exon_chrX->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chrXHet"){$exon_chrXHet->query_partly_within_rect(@rect,\@results);}
    elsif($chrom eq "chrYHet"){$exon_chrYHet->query_partly_within_rect(@rect,\@results);}
    return @results;
  }

  open (GAPS,"$gaps_file") or die "Can't open $gaps_file for reading!\n";
  while (<GAPS>) 
  {
    chomp;
    unless(/^track.+/)
    {
      @elements = split(/\s+/,$_);

      #$gaps{$elements[0]} = 1;

      $chrom = $elements[0];
      $begin = $elements[1];
      $end = $elements[2];

      build_gap_tree($chrom,$begin,$end);
    }
  }
  close GAPS;

  open (EXONS,"$constEx_file") or die "Can't open $constEx_file for reading!\n";
  while (<EXONS>) 
  {
    chomp;
    unless (/^track.+/)
    {
      @elements = split(/\s+/,$_); 
      #$exons{$elements[0]} = 1;

      $chrom = $elements[0];
      $begin = $elements[1];
      $end = $elements[2];
      $gene = $elements[3];

      build_exon_tree($chrom,$begin,$end,$gene)
    }
  }
  close EXONS;

  #print "gap chromosomes\n";
  #foreach $chrom (sort keys %gaps){print "$chrom\n";}

  #print "constitutive exon chromosomes\n";
  #foreach $chrom (sort keys %exons){print "$chrom\n";}

  open(BED,"$bedfile") or die "Can't open $bedfile for reading!\n";
  #open(OUT,">$read_out") or die "Error writing to $read_out!\n";
  
  LINE:while(<BED>)
  {
    chomp;
    $line = $_;
    
    if ($line =~ /^track.+/){next LINE;} 
    else 
    {
      @elements = split(/\s+/,$line); 
      $chrom = $elements[0];
      $start = $elements[1];
      $stop = $elements[2];
      $read = $elements[3];

      @results = query_gaps($chrom,$start,$stop);
	    
      if(@results) 
      {
        for $object (@results) {print "$object\n";}
        #print OUT "gap_$chrom\t$begin\t$end\t$read\n";
        #next LINE; 
      }

      @results = query_exons($chrom,$start,$stop);
	    
      if(@results) 
      {
        for $object (@results) {print "$object\n";}
        #print OUT "gap_$chrom\t$begin\t$end\t$read\n";
        #next LINE; 
      }

      #print OUT "noGap_noExon\t\*\t\*\t$read\n";
      ##This will only happen if neither the gap or exon requirements are met       
    } 
  }
  close BED;
  #close OUT;
}
