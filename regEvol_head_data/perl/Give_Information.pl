#! /usr/bin/perl

#use 5.010;
use strict;
use warnings;

my($line);
my(@info);

if(!open INPUT, "<$ARGV[0]"){
  die "Cannon open file $ARGV[0]";
}


while(<INPUT>){
  chomp;
  if($_ =~ /^\s*$/){
    print "No\n";
  }
  else{
    @info = split("@",$_);
    print $info[1]."_".$info[2].":".$info[3]."\n";
  }
}
