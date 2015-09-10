#!/usr/bin/perl

use strict;
use warnings;
use POSIX;

our($i,$floor,$integer,$ceil);

for ($i = 0; $i <= 7; $i += 0.25)
{
  print "\ni = $i\t";
  printf("printf = %.0f\t",$i);
  $floor = floor($i);
  print "floor = $floor\t";
  $ceil = ceil($i);
  print "ceil = $ceil\t";
  $integer = int($i+0.5);
  print "int = $integer";
}
print "\n\n";
