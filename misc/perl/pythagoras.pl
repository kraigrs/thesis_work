#!/usr/bin/perl

use strict;
use warnings;

main();
sub main
{
  print "\nWhat is the length of the bottom of the triangle? --> ";
  my$s1 = <STDIN>;

  print "\nWhat is the length of the side of the triangle? --> ";
  my$s2 = <STDIN>;

  my$hyp = sqrt($s1**2 + $s2**2); # calculates the hypotenuse
  print "\nThe length of the hypotenuse is $hyp\n\n";
}
