# Usage: R --slave --args <n> <file> <dest> < imprint_draw.r > ./trash 

n <- commandArgs()[4];
#n <- 130;

file <- commandArgs()[5];
dest <- commandArgs()[6];

all <- read.table(file);

sub <- sample(all[,1],n);
sub <- data.frame(sub);

write.table(sub,file=paste(dest,"/temp_imprint.txt",sep=""),quote=FALSE,sep="",row.names=FALSE,col.names=FALSE);
