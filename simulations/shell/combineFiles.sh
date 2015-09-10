#! /bin/sh

if [ $# != 3 ]
then
  echo "You are required to provide 3 files!"
fi

File1=$1
File2=$2
File3=$3

cat $File2 $File3 | sort -t@ | uniq >$File1.wholelist.txt