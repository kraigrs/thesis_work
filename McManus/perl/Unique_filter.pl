#!/usr/bin/perl
## This script condenses the gap-filtered alignment results in preparation for gene mapping.  This uses several scripts from the
## scriptome project (condenser, IDcounter, UnqID_picker, and Merge_Unq_ID).  Condenser eliminates repeats, so reads that mapped 
## to multiple places on the sechellia genome, but corresponded to only one place on the melanogaster genome after liftover are 
## only present once.  IDcounter counts the number of times each read is present.  UnqID_picker chooses readIDs that are present
## only one time.  Merge_Unq_ID chooses alignment data for the unique reads from the condensed file created by condenser.pl.

system 'perl condenser.pl Hyb_r1_mel_SNP_nogap.BED > Hyb_r1_mel_SNP_nogap_condensed.BED';
system 'perl IDcounter.pl Hyb_r1_mel_SNP_nogap_condensed.BED > Hyb_r1_mel_SNP_nogap_IDcount.txt';
system 'perl UnqID_picker.pl Hyb_r1_mel_SNP_nogap_IDcount.txt > Hyb_r1_mel_SNP_nogap_UnqIDs.txt';
system 'perl Merge_Unq_ID.pl Hyb_r1_mel_SNP_nogap_UnqIDs.txt Hyb_r1_mel_SNP_nogap_condensed.BED > Hyb_r1_mel_SNP_nogap_unq.txt';

system 'perl condenser.pl Hyb_r1_mel_link_nogap.BED > Hyb_r1_mel_link_nogap_condensed.BED';
system 'perl IDcounter.pl Hyb_r1_mel_link_nogap_condensed.BED > Hyb_r1_mel_link_nogap_IDcount.txt';
system 'perl UnqID_picker.pl Hyb_r1_mel_link_nogap_IDcount.txt > Hyb_r1_mel_link_nogap_UnqIDs.txt';
system 'perl Merge_Unq_ID.pl Hyb_r1_mel_link_nogap_UnqIDs.txt Hyb_r1_mel_link_nogap_condensed.BED > Hyb_r1_mel_link_nogap_unq.txt';

system 'perl condenser.pl Hyb_r2_mel_SNP_nogap.BED > Hyb_r2_mel_SNP_nogap_condensed.BED';
system 'perl IDcounter.pl Hyb_r2_mel_SNP_nogap_condensed.BED > Hyb_r2_mel_SNP_nogap_IDcount.txt';
system 'perl UnqID_picker.pl Hyb_r2_mel_SNP_nogap_IDcount.txt > Hyb_r2_mel_SNP_nogap_UnqIDs.txt';
system 'perl Merge_Unq_ID.pl Hyb_r2_mel_SNP_nogap_UnqIDs.txt Hyb_r2_mel_SNP_nogap_condensed.BED > Hyb_r2_mel_SNP_nogap_unq.txt';

system 'perl condenser.pl Hyb_r2_mel_link_nogap.BED > Hyb_r2_mel_link_nogap_condensed.BED';
system 'perl IDcounter.pl Hyb_r2_mel_link_nogap_condensed.BED > Hyb_r2_mel_link_nogap_IDcount.txt';
system 'perl UnqID_picker.pl Hyb_r2_mel_link_nogap_IDcount.txt > Hyb_r2_mel_link_nogap_UnqIDs.txt';
system 'perl Merge_Unq_ID.pl Hyb_r2_mel_link_nogap_UnqIDs.txt Hyb_r2_mel_link_nogap_condensed.BED > Hyb_r2_mel_link_nogap_unq.txt';

system 'perl condenser.pl Mix_r1_mel_SNP_nogap.BED > Mix_r1_mel_SNP_nogap_condensed.BED';
system 'perl IDcounter.pl Mix_r1_mel_SNP_nogap_condensed.BED > Mix_r1_mel_SNP_nogap_IDcount.txt';
system 'perl UnqID_picker.pl Mix_r1_mel_SNP_nogap_IDcount.txt > Mix_r1_mel_SNP_nogap_UnqIDs.txt';
system 'perl Merge_Unq_ID.pl Mix_r1_mel_SNP_nogap_UnqIDs.txt Mix_r1_mel_SNP_nogap_condensed.BED > Mix_r1_mel_SNP_nogap_unq.txt';

system 'perl condenser.pl Mix_r1_mel_link_nogap.BED > Mix_r1_mel_link_nogap_condensed.BED';
system 'perl IDcounter.pl Mix_r1_mel_link_nogap_condensed.BED > Mix_r1_mel_link_nogap_IDcount.txt';
system 'perl UnqID_picker.pl Mix_r1_mel_link_nogap_IDcount.txt > Mix_r1_mel_link_nogap_UnqIDs.txt';
system 'perl Merge_Unq_ID.pl Mix_r1_mel_link_nogap_UnqIDs.txt Mix_r1_mel_link_nogap_condensed.BED > Mix_r1_mel_link_nogap_unq.txt';

system 'perl condenser.pl Mix_r2_mel_SNP_nogap.BED > Mix_r2_mel_SNP_nogap_condensed.BED';
system 'perl IDcounter.pl Mix_r2_mel_SNP_nogap_condensed.BED > Mix_r2_mel_SNP_nogap_IDcount.txt';
system 'perl UnqID_picker.pl Mix_r2_mel_SNP_nogap_IDcount.txt > Mix_r2_mel_SNP_nogap_UnqIDs.txt';
system 'perl Merge_Unq_ID.pl Mix_r2_mel_SNP_nogap_UnqIDs.txt Mix_r2_mel_SNP_nogap_condensed.BED > Mix_r2_mel_SNP_nogap_unq.txt';

system 'perl condenser.pl Mix_r2_mel_link_nogap.BED > Mix_r2_mel_link_nogap_condensed.BED';
system 'perl IDcounter.pl Mix_r2_mel_link_nogap_condensed.BED > Mix_r2_mel_link_nogap_IDcount.txt';
system 'perl UnqID_picker.pl Mix_r2_mel_link_nogap_IDcount.txt > Mix_r2_mel_link_nogap_UnqIDs.txt';
system 'perl Merge_Unq_ID.pl Mix_r2_mel_link_nogap_UnqIDs.txt Mix_r2_mel_link_nogap_condensed.BED > Mix_r2_mel_link_nogap_unq.txt';


