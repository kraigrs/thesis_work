#!/usr/bin/perl

###################################################################################
#
# 02/28/2013
#
# remove_genes.pl
#
###################################################################################

use strict;
use warnings;

my$remove = $ARGV[0];
my$introns = $ARGV[1];
my$intergenics = $ARGV[2];

my$gene1; my$gene2; my$intron1; my$intron2; my$inter; my$name; my$line;
my@elements; my@meta;
my%list;

open(REM,"$remove") or die "Can't open $remove for reading!";
while(<REM>) 
{
  chomp;
  #print "[[[$_\]]]\n"; exit;
  $list{$_} = 1;
}
close REM;

open(INTRON,"$introns") or die "Can't open $introns for reading!";
while(<INTRON>) 
{
  chomp;
  $line = $_;
  @elements = split(/\s/,$line);
  $name = $elements[3];

  @meta = split(/\_/,$name);
  $intron1 = $meta[1];
  $intron2 = $meta[2];
  if($intron1 =~ /(\w+)\:\d+/){$gene1 = $1;}
  if($intron2 =~ /(\w+)\:\d+/){$gene2 = $1;}

  if($gene1 eq $gene2 && $intron1 ne $intron2)
  {
    unless($list{$gene1}){print "$line\n";}
  }
}
close INTRON;

open(INTER,"$intergenics") or die "Can't open $intergenics for reading!";
HERE: while(<INTER>) 
{
  chomp;
  $line = $_;
  @elements = split(/\s/,$line);
  $name = $elements[3];

  @meta = split(/\_/,$name);

  if(scalar(@meta) < 3)
  {
    if($meta[0] =~ /inter/ && $meta[1] =~ /CG\d+/)
    {
      $gene1 = $meta[1];
      unless($list{$gene1}){print "$line\n"; next HERE;}
    }
    elsif($meta[0] =~ /CG\d+/ && $meta[1] =~ /inter/)
    {
      $gene1 = $meta[0];
      unless($list{$gene1}){print "$line\n"; next HERE;}
    }
  }
  else
  {
    $gene1 = $meta[0];
    $gene2 = $meta[2];
    
    if($gene1 =~ /CG\d+/ || $gene2 =~ /CG\d+/)
    {
      if($gene1 =~ /CG\d+/ && !$list{$gene1}){print "$line\n"; next HERE;}
      elsif($gene2 =~ /CG\d+/ && !$list{$gene2}){print "$line\n"; next HERE;}
    }
  }
}
close INTER;
