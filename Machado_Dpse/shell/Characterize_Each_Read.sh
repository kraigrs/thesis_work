#!/bin/sh

#PBS -S /bin/sh
#PBS -N Characterize_Each_Read
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=1,walltime=1:00:00
#PBS -j oe
#PBS -V

####  End PBS preamble
#  Include the next three lines always

if [ "x${PBS_NODEFILE}" != "x" ] ; then
   cat $PBS_NODEFILE   # contains a list of the CPUs you were using if run with PBS
fi

cd $PBS_O_WORKDIR

# sort file files according to metadata

sh sortFile.sh $var1 $var2 $var3 $var4

# remove duplicates

#perl ../perl/Eliminate_Duplicates.pl $var1.sorted $var1.sorted.undup
#perl ../perl/Eliminate_Duplicates.pl $var2.sorted $var2.sorted.undup
#perl ../perl/Eliminate_Duplicates.pl $var3.sorted $var3.sorted.undup
#perl ../perl/Eliminate_Duplicates.pl $var4.sorted $var4.sorted.undup

# produce the big file

sh combineFiles.sh $var5 $var1.sorted $var2.sorted $var3.sorted $var4.sorted

# get the read specific information

sh readsExist.sh $var5 $var1.sorted $var2.sorted $var3.sorted $var4.sorted

# get the final information

perl ../perl/Remove_Duplicate_Metadata.pl $var5.final_information.txt $var5.final_info_rmDup.txt

#rm $var5.final_information.txt $var5.wholelist.txt
