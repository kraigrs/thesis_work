#!/bin/sh

FILES=./*.pdf
for pdf in $FILES
do
  pdf2ps $pdf
  tmp=${pdf%*.pdf}
  #echo $tmp.ps
  ps2eps $tmp.ps
done