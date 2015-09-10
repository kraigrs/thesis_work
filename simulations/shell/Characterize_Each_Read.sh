#!/bin/sh

#PBS -S /bin/sh
#PBS -N Characterize_Each_Read
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=1,walltime=1:00:00
#PBS -m abe
#PBS -j oe
#PBS -V

####  End PBS preamble
#  Include the next three lines always

if [ "x${PBS_NODEFILE}" != "x" ] ; then
   cat $PBS_NODEFILE   # contains a list of the CPUs you were using if run with PBS
fi

cd $PBS_O_WORKDIR

# sort file files according to metadata
sh sortFile.sh $var1 $var2

# remove duplicates
#perl ../perl/Eliminate_Duplicates.pl $var1.sorted $var1.sorted.undup
#perl ../perl/Eliminate_Duplicates.pl $var2.sorted $var2.sorted.undup

# produce the big file
#sh combineFiles.sh $1.sort.undup $2.sort.undup
sh combineFiles.sh $var3 $var1.sorted $var2.sorted

# get the read specific information
#sh readsExist.sh $work_dir/wholelist.txt $1.sort.undup $2.sort.undup
sh readsExist.sh $var3 $var1.sorted $var2.sorted

# get the final information
perl ../perl/Remove_Duplicate_Metadata.pl $var3.final_information.txt $var3.final_info_rmDup.txt

#rm $var3.final_information.txt $var3.wholelist.txt
