#! /usr/bin/perl

# This script is to eliminate duplicates in files after intersectBed step. If two lines are exactly the same,
# then we keep only one of them. If two lines are only same in metadata column, we exclude both of them. 

use strict;
use warnings;
#use 5.010;

my($inputfile);
my(@info);
my(@sameMeta);
my($line);
my($lastMeta, $lastStart, $lastEnd, $lastLine, $repeat, $llStart, $llEnd, $llMeta);
my $First_Read = 1;
my $Kick_Out = 0;


$inputfile = $ARGV[0]; #get the input filehandle.

if(!open INPUT, "<$inputfile"){
  die "Cannot open file $inputfile";
}
open OUT, ">$ARGV[1]";

while(<INPUT>){
  #print 1;
  chomp($line = $_);
  @info = split("@",$line);
  if($First_Read){                   # If it is first read, push it into the stack. If next read is same with the first, then it does not matter. While not,
    $lastMeta = $info[0];            # $Kick_Out is still zero, so it will just shift the first read and print it. 
    $lastStart = $info[1];
    $lastEnd = $info[2];
    $lastLine = $line;
    push @sameMeta, $line;
    $First_Read = 0;
    next;
  }
  
  if($info[0] ne $lastMeta){         # If read does not have same as previous one, print last line out if stack is empty, or using information in $Kick_Out to 
    if($#sameMeta == 0){             # print out first read in stack or just throw all of them away.
      print OUT $lastLine."\n";
      #print $lastLine."\n";
    }
    else{
      if($Kick_Out){
	while($#sameMeta > 0){
	  pop @sameMeta;
	  
	}
	$Kick_Out = 0;
      }
      else{
	$repeat = pop @sameMeta;
	print OUT $repeat."\n";
	#print $repeat."\n";
	#print $#sameMeta."\n";
	while($#sameMeta > 0){
	  pop @sameMeta;
	}
        #push @sameMeta, $line;
      }
    }
  }
  
  else{
    push @sameMeta, $line;
    if(!$Kick_Out && ($info[1] != $lastStart || $info[2] != $lastEnd)){
      $Kick_Out = 1;
    }
  }
  
  $llMeta = $lastMeta;
  $llStart = $lastStart;
  $llEnd = $lastEnd;
  $lastLine = $line;
  $lastMeta = $info[0];
  $lastStart = $info[1];
  $lastEnd = $info[2];
}

#print "\nKickout = $Kick_Out\n";
if($#sameMeta){                       
  if(!$Kick_Out){                     
    $repeat = pop @sameMeta;      # changed from shift and it worked 
    print OUT $repeat."\n";                
  }
}
else{
  print OUT $line."\n";
}

  
  
  
	

    
      
  


