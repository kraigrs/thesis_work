#!/usr/bin/perl

########################################################################
# 
# 04/09/2010
#
# liftMap.pl
#
# Purpose: using the .map file of sechellia reads and the .bed liftover file,
#          compile complete list of reads and their liftover coordinates
# 
# Input: 
#
# Output: 
#
########################################################################

use strict;
#use warnings;

my$secmap = $ARGV[0];
my$secbed = $ARGV[1];
my$seclifted = $ARGV[2];

my$mapline; my$mapread; my$mapcontig; 
my$bedline; my$bedread; my$bedcontig;

my@mapfields;  
my@bedfields; 
 
open(MAP,"$secmap") or die "Can't open $secmap for reading!";
open(BED,"$secbed") or die "Can't open $secbed for reading!";

open(OUT,">$seclifted") or die "Error writing to $seclifted\n";

$bedline = <BED>;
chomp $bedline;
@bedfields = split(/\s+/,$bedline);  
$bedcontig = $bedfields[0];
$bedread = $bedfields[3];

HERE:while(<MAP>)
{
  chomp;
  $mapline = $_;
  @mapfields = split(/\s+/,$mapline);  
  $mapcontig = $mapfields[0];
  $mapread = $mapfields[3];

  if($mapcontig eq "*")
  {
    print OUT "$mapline\n"; 
    next HERE;    
  }
  elsif($mapread eq $bedread)
  {
    print OUT "$bedline\n"; 
    
    # maybe add chomp!
    $bedline = <BED>;
    chomp $bedline;
    @bedfields = split(/\s+/,$bedline);  
    $bedcontig = $bedfields[0];
    $bedread = $bedfields[3];  
    next HERE;  
  }
  else
  {
    print OUT "no_lift\t\*\t\*\t$mapread\n"; 
    next HERE;
  }
}
close OUT;
close MAP;



