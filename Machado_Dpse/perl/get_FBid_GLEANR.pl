#!/usr/bin/perl

use strict;
use warnings;

my$list = $ARGV[0];

print "gene\tGLEANR\n";

open(LIST,"$list") or die "\nError opening $list\n";
while(<LIST>)    
{
  chomp;
  if(/^\>.+FlyBase\:(\w+)\,.+GLEANR\:(\w+)\,/){print "$1\t$2\n";}
}
close LIST;
