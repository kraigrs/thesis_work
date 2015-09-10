#! /bin/sh

if [ $# != 4 ]
then
  echo "You are required to provide 4 files!"
fi

File1=$1
File2=$2
File3=$3
File4=$4


if [[ $File1 =~ ([[:print:]]*)/[^/]*$ ]]
then
  work_dir=${BASH_REMATCH[1]}
fi

cat $File1 $File2 $File3 $File4 | sort -t@ | uniq >$work_dir/wholelist.txt