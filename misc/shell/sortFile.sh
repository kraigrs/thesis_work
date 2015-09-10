#!/bin/sh

awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $1 | sed 's/\/[12]//g' | sort >$1.sort
awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $2 | sed 's/\/[12]//g' | sort >$2.sort
awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $3 | sed 's/\/[12]//g' | sort >$3.sort
awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $4 | sed 's/\/[12]//g' | sort >$4.sort
