#!/usr/bin/perl

###################################################################################
#
# 01/13/2010
#
# parseHTML.pl
#
# Purpose: go through an HTML file and pull out links for gene names and make a hash 
#          for all those genes
# 
# Input: the HTML file containing CG- names and their known aliases 
#
###################################################################################

use strict;
use warnings;

my$starttime = time;

my$HTML = $ARGV[0];

my@elements;
my@genes;
my$i;
my%geneHash;

open(HTML,"$HTML") or die "\nError opening $HTML\n";
while(<HTML>)
{
  chomp;
  @elements = split(/></,$_);
  foreach(@elements)
  {
    #print "$_\t";
    if(/>(\S+)</){push(@genes,$1);}   
  }
  #foreach(@genes){print "$_\t";}
  #print "\n";
  for($i=1;$i<@genes;$i++)
  {
    #print "\n$genes[0]\t$genes[$i]\n\n"; 
    $geneHash{$genes[$i]} = $genes[0];
  }
  splice @genes;
}
close HTML;

foreach(keys %geneHash){print "Gene: $geneHash{$_} Alias: $_\n";}

