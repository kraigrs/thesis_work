#!/bin/sh

if [ $# != 5 ]
then
  echo "You are required to provide 5 files!"
fi

File1=$1
File2=$2
File3=$3
File4=$4
File5=$5

cat $File1 $File2 $File3 $File4 | sort -t@ | uniq >$File5.wholelist.txt