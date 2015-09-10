#!/bin/sh

#PBS -S /bin/sh
#PBS -N characterize
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l procs=1,pmem=4000mb,walltime=5:00:00
#PBS -j oe
#PBS -V

####  End PBS preamble
#  Include the next three lines always

if [ "x${PBS_NODEFILE}" != "x" ] ; then
   cat $PBS_NODEFILE   # contains a list of the CPUs you were using if run with PBS
fi

cd $PBS_O_WORKDIR

echo DIR = $var1

# sort file files according to metadata

#sh sortFile.sh $1 $2 $3 $4

sh sortFile.sh $var1\_m1.$var2.bowtie_v0_m1.merged.bed.lifted.converted $var1\_m1.$var3.bowtie_v0_m1.merged.bed.lifted.converted $var1\_m2.$var2.bowtie_v0_m1.merged.bed.lifted.converted $var1\_m2.$var3.bowtie_v0_m1.merged.bed.lifted.converted


# remove duplicates

#perl ../perl/Eliminate_Duplicates.pl $1.sorted $1.sorted.undup
#perl ../perl/Eliminate_Duplicates.pl $2.sorted $2.sorted.undup
#perl ../perl/Eliminate_Duplicates.pl $3.sorted $3.sorted.undup
#perl ../perl/Eliminate_Duplicates.pl $4.sorted $4.sorted.undup

#perl ../perl/Eliminate_Duplicates.pl $var1\_m1.$var2.bowtie_v0_m1.merged.bed.lifted.converted.sorted $var1\_m1.$var2.bowtie_v0_m1.merged.bed.lifted.converted.sorted.undup
#perl ../perl/Eliminate_Duplicates.pl $var1\_m1.$var3.bowtie_v0_m1.merged.bed.lifted.converted.sorted $var1\_m1.$var3.bowtie_v0_m1.merged.bed.lifted.converted.sorted.undup
#perl ../perl/Eliminate_Duplicates.pl $var1\_m2.$var2.bowtie_v0_m1.merged.bed.lifted.converted.sorted $var1\_m2.$var2.bowtie_v0_m1.merged.bed.lifted.converted.sorted.undup
#perl ../perl/Eliminate_Duplicates.pl $var1\_m2.$var3.bowtie_v0_m1.merged.bed.lifted.converted.sorted $var1\_m2.$var3.bowtie_v0_m1.merged.bed.lifted.converted.sorted.undup

# produce the big file

#sh combineFiles.sh $1.sorted.undup $2.sorted.undup $3.sorted.undup $4.sorted.undup

sh combineFiles.sh $var1\_m1.$var2.bowtie_v0_m1.merged.bed.lifted.converted.sorted $var1\_m1.$var3.bowtie_v0_m1.merged.bed.lifted.converted.sorted $var1\_m2.$var2.bowtie_v0_m1.merged.bed.lifted.converted.sorted $var1\_m2.$var3.bowtie_v0_m1.merged.bed.lifted.converted.sorted $var1

# get the read specific information

#sh readsExist.sh $work_dir/wholelist.txt $1.sorted.undup $2.sorted.undup $3.sorted.undup $4.sorted.undup

sh readsExist.sh $var1 $var1\_m1.$var2.bowtie_v0_m1.merged.bed.lifted.converted.sorted $var1\_m1.$var3.bowtie_v0_m1.merged.bed.lifted.converted.sorted $var1\_m2.$var2.bowtie_v0_m1.merged.bed.lifted.converted.sorted $var1\_m2.$var3.bowtie_v0_m1.merged.bed.lifted.converted.sorted


# get the final information

#perl ../perl/Remove_Duplicate_Metadata.pl $work_dir/final_information.txt $work_dir/final_info_rmDup.txt

perl ../perl/Remove_Duplicate_Metadata.pl $var1.final_information.txt $var1.final_info_rmDup.txt

#rm $var1.final1 $var1.final2 $var1.final3 $var1.final_information.txt $var1.wholelist.txt
#rm $var1*sort*