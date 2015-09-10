################################################
#
# Calculate memory requirements for Velvet
#
# based on http://listserver.ebi.ac.uk/pipermail/velvet-users/2009-July/000474.html
#
################################################


GenomeSize <- 120; # in Mb
ReadSize <- 100; # length of reads
NumReads <- 30; # in millions of reads
K <- 31; # k-mer length

#velvetg = (-109635 + 18977*ReadSize + 86326*GenomeSize + 233353*NumReads - 51092*K)/1048576; # Gb of memory required
(-109635 + 18977*ReadSize + 86326*GenomeSize + 233353*NumReads - 51092*K)/1048576;