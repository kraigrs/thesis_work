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



if [ $# != 3 ]
then
    echo "You are required to give three files!"
    exit 1
fi

FINAL=$1.final_information.txt

comm $2 $1.wholelist.txt | cut -f3 >$1.tmp
perl ../perl/Give_Information.pl $1.tmp >$1.tmp1
paste $1.wholelist.txt $1.tmp1 >$1.final1
rm -f $1.tmp $1.tmp1

comm $3 $1.wholelist.txt | cut -f3 >$1.tmp
perl ../perl/Give_Information.pl $1.tmp >$1.tmp1
paste $1.final1 $1.tmp1 >$FINAL
rm -f $1.tmp $1.tmp1 $1.final1
