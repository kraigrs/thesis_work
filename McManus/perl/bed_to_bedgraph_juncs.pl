#!/usr/bin/perl
#
# bed_to_wig_for_li.pl
# Li 10/26/09

# Convert bed output file (including jctn alignments)
# to wig (wiggle) representation

# This version is supposed to be easy on the RAM, writing temporary per-chromosome bracket files
# which are sorted and then processed to generate the wig file. (to be done)

# ("wig" plots are actually bed-graph plots that have starts that are 0-based and ends that are 1-based)

# *****IT IS ASSUMED THAT THE bed ALIGNMENTS WERE DONE WITH -B 1 (1-based coordinate output)*****


# *****ALSO BE SURE TO PARSE THE READING OF JCTN READS ($chromo field)******* (lines 86-91 --- check $f[ ] indices matchup with jctn interval 
# used in jctn db hearder format)
# This is typically different for different species (depends on jctn db header line format)

# Invoke the script via:  "perl bed_to_wig_human2.pl"
# Output file will be filename.wig



use strict;
use FileHandle;

my @bed_output_file;

#my $bed_out_path = "/data/pipeline/SPA/FOR_CARMICHAEL/bed/UNTRIMMED/K_1_M_10_V_2";           #bed alignment files specified here
my $bed_out_path = "./";
push(@bed_output_file, "$bed_out_path$ARGV[0]");

my @f;
my @ff;
my $line;
my $chromo;
my $start;
my $end;
my $start1;
my $end1;
my $start2;
my $end2;

my $read_length;

my $wig_name;

#---------------------LOOP OVER bed OUTPUT FILES----------------------------------------------------------------

for (my $k=0; $k <= $#bed_output_file; $k++){

  open(bed_OUT, "<$bed_output_file[$k]") or die("can't open $bed_output_file[$k]"); #OPEN bed ALIGNMENT OUTPUT FILE

  $wig_name = $bed_output_file[$k]; 
 # $wig_name =~ s/bed/WIG/;
  $wig_name =~ s/\.bed/\.wig/i;

  my $WIG = new FileHandle ">$wig_name" or die("can't open $wig_name"); #WIGGLE BROWSER TRAK FILE 

  @f = split /\//,$wig_name;
  my $trak_name = $f[-1];

  #wig file header line 

  $WIG->printf("track type=wiggle_0 name=\"$trak_name"."\" description=\"$trak_name\" color=51,0,255 autoScale=on visibility=full\n");     #purple counts



  my %read_boundary;        # position of [ or ] boundary
  my %read_boundary_type;   #  +1 => [   -1 => ]   

  my $number_of_reads = 0;

  my %temp_chromo_bracket_files;   #temporary chromo-specific bracket files

  while ($line = <bed_OUT>){    #read bed alignment data lines
    $number_of_reads++;

    chomp($line);
    @f = split /\s+/,$line;
    $chromo   = $f[0];
    if (index($chromo,"_") != -1) {                  #if jctn read (detected via appearence of "_" somewhere in the chromosome field)
      @ff = split /_/,$chromo;                       #chromo for jctns is: chromo1_start1_end1_chromo2_start2_end2
      $chromo = $ff[0];                              #                      [0]      [1]   [2]   [3]    [4]    [5]

      $start1 = $ff[1]+$f[1];   #NOTE: coords are all assumed to be true-coords
      $end1 = $ff[2];
      $start2 = $ff[3];
      $end2 = $start2+($f[2]-$f[1])-($end1-$start1-2);   #TRICKY! (everything is 1-based)

      push(@{$read_boundary{$chromo}},$start1);   #chunk1
      push(@{$read_boundary_type{$chromo}},1);
      push(@{$read_boundary{$chromo}},$end1);
      push(@{$read_boundary_type{$chromo}},-1);
      push(@{$read_boundary{$chromo}},$start2);   #chunk2
      push(@{$read_boundary_type{$chromo}},1);
      push(@{$read_boundary{$chromo}},$end2);
      push(@{$read_boundary_type{$chromo}},-1);
    }
    else{                                            #normal (non-junction) read
      $start = $f[1]+1;
      $end   = $f[2];
      push(@{$read_boundary{$chromo}},$start);   #normal read
      push(@{$read_boundary_type{$chromo}},1);
      push(@{$read_boundary{$chromo}},$end);
      push(@{$read_boundary_type{$chromo}},-1);

    }
  }



  foreach my $chromo (keys %read_boundary) {

    printf("Sorting by position and type for chromosome $chromo \n");

# sort read_boundary (and read_boundary_type) ----- ascending order in pos, descending in type (i.e., [ comes before ] if position is common)

    @{$read_boundary{$chromo}} = @{$read_boundary{$chromo}}[
  		my @idx = sort {$read_boundary{$chromo}[$a] <=> $read_boundary{$chromo}[$b]  ||
  		$read_boundary_type{$chromo}[$b] <=> $read_boundary_type{$chromo}[$a]} 
        0..$#{@{$read_boundary{$chromo}}}];

    @{$read_boundary_type{$chromo}}   = @{$read_boundary_type{$chromo}}[@idx];


# Ok, so all read brackets are sorted within each chromo
# THIS IS THE MAIN ALGORITHM:
#------------------------------------------------------------------------------------------------


    my $count=0;
    my $start;
    my $end;
    my $bracket;        #bracket position
    my $bracket_type;   # [ => 1      ] => -1
    my $old_bracket=0;
    my $old_bracket_type=0;
    my $start_zero_based;

    for (my $k=0; $k<= $#{@{$read_boundary{$chromo}}}; $k++){                      #traverse bracket array, writing to .wig file
    
      $bracket = $read_boundary{$chromo}[$k];

      $bracket_type = $read_boundary_type{$chromo}[$k];


      if( $old_bracket_type == 1){                                                 #start
        $start= $old_bracket;                                    
      }
      else{
       $start= $old_bracket+1;
      }

      if ($bracket != $old_bracket  || $bracket_type != $old_bracket_type){         #end
        if ($bracket_type == 1 && $count != 0  && $start != $bracket){              #last bit is hack for abutting reads
          $start--;                                                                 #hack to account for BED records being zero-based half-open
          $end = $bracket-1;
          $WIG->printf("$chromo\t$start\t$end\t$count\n");
        }
        if ($bracket_type == -1){
          $start--;                                                                 #hack to account for BED records being zero-based half-open
	  $end = $bracket;
          $WIG->printf("$chromo\t$start\t$end\t$count\n");
        }
      }



      $count = $count + $bracket_type;
      $old_bracket = $bracket;
      $old_bracket_type = $bracket_type;

    }  #next bracket




#---------------------------------------------------------------------------------------------

 
  } #next chromo

  printf("Number of reads = $number_of_reads\n");

}  #next bed alignment file


