#! /bin/sh

if [ $# != 4 ]
then
	echo "You are required to provide 4 files!"
	exit 1
fi

# get the folder name
if [[ $1 =~ ([[:print:]]*)/[^/]*$ ]]
then
	work_dir=${BASH_REMATCH[1]}
fi

# sort file files according to metadata
sh sortFile.sh $1 $2 $3 $4

# remove duplicates
perl ../perlScripts/Eliminate_Duplicates.pl $1.sort $1.sort.undup
perl ../perlScripts/Eliminate_Duplicates.pl $2.sort $2.sort.undup
perl ../perlScripts/Eliminate_Duplicates.pl $3.sort $3.sort.undup
perl ../perlScripts/Eliminate_Duplicates.pl $4.sort $4.sort.undup

# produce the big file
sh combineFiles.sh $1.sort.undup $2.sort.undup $3.sort.undup $4.sort.undup

# get the read specific information
sh readsExist.sh $work_dir/wholelist.txt $1.sort.undup $2.sort.undup $3.sort.undup $4.sort.undup

# get the final information
perl ../perlScripts/Remove_Duplicate_Metadata.pl $work_dir/final_information.txt $work_dir/final_info_rmDup.txt