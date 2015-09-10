#!/bin/sh

#####################################################################
#
# Purpose: this code is for getting the allel specific information of
#          exon. The argument you should give is the file produced by
#          intersectBed between alignment and consecutive_region. It 
#          give you a file contain the counts of each reads that is in
#          any of the four files.
#
# Syntax: Get_allel_specific_information.csh your four argument files
#         The files should be in the order of StrainA_mate_1 StrainA_
#         mate_2 StrainB_mate_1 StrainB_mate_2
# 
# Author:ypauling@umich.edu
#
#####################################################################



if [ $# != 5 ]
then
    echo "You are required to give five files!"
    exit 1
fi

if [[ $1 =~ ([[:print:]]*)/[^/]*$ ]]
then
    work_dir=${BASH_REMATCH[1]}
fi

FINAL=$work_dir/final_information.txt

comm $2 $1 | cut -f3 >$work_dir/tmp
perl ../perlScripts/Give_Information.pl $work_dir/tmp >$work_dir/tmp1
paste $1 $work_dir/tmp1 >$work_dir/final1
rm -f $work_dir/tmp $work_dir/tmp1

comm $3 $1 | cut -f3 >$work_dir/tmp
perl ../perlScripts/Give_Information.pl $work_dir/tmp >$work_dir/tmp1
paste $work_dir/final1 $work_dir/tmp1 >$work_dir/final2
rm -f $work_dir/tmp $work_dir/tmp1

comm $4 $1 | cut -f3 >$work_dir/tmp
perl ../perlScripts/Give_Information.pl $work_dir/tmp >$work_dir/tmp1
paste $work_dir/final2 $work_dir/tmp1 >$work_dir/final3
rm -f $work_dir/tmp $work_dir/tmp1

comm $5 $1 | cut -f3 >$work_dir/tmp
perl ../perlScripts/Give_Information.pl $work_dir/tmp >$work_dir/tmp1
paste $work_dir/final3 $work_dir/tmp1 >$FINAL
rm -f $work_dir/tmp $work_dir/tmp1
