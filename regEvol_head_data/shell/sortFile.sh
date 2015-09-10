#!/bin/sh

awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $1 | sort >$1.sorted
awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $2 | sort >$2.sorted
awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $3 | sort >$3.sorted
awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $4 | sort >$4.sorted
