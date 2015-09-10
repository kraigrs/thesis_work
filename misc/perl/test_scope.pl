#!usr/bin/perl

use strict;
use warnings;

my$test = 0;
my$counter = 0;

while($counter < 10)
{
  $test += 5;
  $counter += 1;
  print "Counter = $counter\nTest = $test\n\n";
}

print "Final_counter = $counter\nFinal_test = $test\n\n";
