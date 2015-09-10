data <- read.table("/Users/kraigrs/Desktop/imp.txt");

perm <- 1000;
i = 0;

while(i < perm)
{
	data <- cbind(data,sample(data[,2],nrow(data)));
	i = i + 1;
}

write.table(data,file="/Users/kraigrs/Desktop/permuted_imprinting.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE);