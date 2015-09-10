#!/bin/sh

#awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $1 | sort >$1.sorted
#awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $2 | sort >$2.sorted
#awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $3 | sort >$3.sorted
#awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $4 | sort >$4.sorted

awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $1 | sed 's/\/[12]//g' | sort >$1.sorted
awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $2 | sed 's/\/[12]//g' | sort >$2.sorted
awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $3 | sed 's/\/[12]//g' | sort >$3.sorted
awk 'BEGIN{FS="\t";OFS="@"}{print $4,$2,$3,$1}' $4 | sed 's/\/[12]//g' | sort >$4.sorted
