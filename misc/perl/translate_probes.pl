#!/usr/bin/perl

use strict;
use warnings;

my$key = $ARGV[0];
my$link = $ARGV[1];
my$probe_list = $ARGV[2];
my$out = $ARGV[3];

our(@elements);
our(%FB2gene,%probe2FB);

open(KEY,"$key") or die "\nError opening $key\n";
while(<KEY>)    
{
  chomp;
  @elements = split(/\t/,$_); 
  $FB2gene{$elements[0]} = $elements[2]; # FB2gene {submitted FB id} = gene name 
}
close KEY;

open(LINK,"$link") or die "\nError opening $link\n";
while(<LINK>)    
{
  chomp;
  @elements = split(/\t/,$_);
  unless(@elements < 2){$probe2FB{$elements[0]} = $elements[1];} # probe2FB {probe name} = submitted FB id 
}
close LINK;

open(OUT,">$out") or die "\nError opening $out\n";
print OUT "probe\tgene";

open(LIST,"$probe_list") or die "\nError opening $probe_list\n";
while(<LIST>)    
{
  chomp;
  if($probe2FB{$_})
  {
    if($FB2gene{$probe2FB{$_}})
    {
      #print "\n$_\t$FB2gene{$probe2FB{$_}}";
      print OUT "\n$_\t$FB2gene{$probe2FB{$_}}";
    }
  } 
}
close LIST;

close OUT;
