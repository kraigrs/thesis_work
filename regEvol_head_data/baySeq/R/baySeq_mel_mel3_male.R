ptm <- proc.time()

####################################################################################################
####################################################################################################
#data analyses for mel-mel hybrid 1 vs 2 male
####################################################################################################
####################################################################################################

setwd("/scratch/lsa_flux/kraigrs/regEvol_head_data/baySeq/mel_mel")
#setwd("/Users/kraigrs/Wittkopp/regEvol_head_data/baySeq/mel_mel")

####################################################################################################
#model generation
####################################################################################################

#Updated correct counting function
N <- function(n) {
	if(n == 0 | n == 1) {
		return(1)
	}
	out <- choose(n-1,0)
	I <- 1:(n-1)
	for(i in I){
		out <- out + choose(n-1,i)*2^(i-1)
	}     
	return(out)
}

#Convert integers to binary
int.to.bin <-
function(x, b=2){
	xi <- as.integer(x)
	N <- length(x)
	xMax <- max(x)	
	ndigits <- (floor(logb(xMax, base=2))+1)
	Base.b <- array(NA, dim=c(N, ndigits))
	for(i in 1:ndigits){#i <- 1
		Base.b[, ndigits-i+1] <- (x %% b)
		x <- (x %/% b)
	}
	if(N ==1) Base.b[1, ] else Base.b
}

#Generate groupings
M <- function(n) {
	if(n == 2) {
		return(matrix(c(0,0,0,1), nrow = 2))
	}
	x <- int.to.bin(0:((2^(n-1))-1))
	I <- 1:dim(x)[1]
	for(i in I) {
		sum <- sum(x[i,])
		if(sum > 1) {
			#R <- int.to.bin(0:((2^(sum-1))-1))
			R <- M(sum)
			#v <- vector(length = 1, mode = "integer") 
			#R <- cbind(v,R)
			J <- 2:dim(R)[1]
			K <- 1:dim(x)[2] 
			for(j in J) {
				l <- 1
				h <- vector(length = length(K), mode = "integer")
				for(k in K) {
					if(x[i,k] == 1) {
						h[k] <- x[i,k] + R[j,l]
						l <- l + 1
					}else {
						h[k] <- x[i,k]
					}
				}
				x <- rbind(x,h)
			}
		}
	}
	s <- vector(length = dim(x)[1], mode = "integer")
	x <- cbind(s,x)
	return(x)
}

model1<-M(4)

model2<-cbind(model1[,1],model1[,1],model1[,1],model1[,2],model1[,2],model1[,2],model1[,3],model1[,3],model1[,3],model1[,4],model1[,4],model1[,4])
model2<-model2+1

models <- vector("list", nrow(model2))
for (i in 1:nrow(model2)) {
    models[[i]] <- model2[i,]
	}

rm(int.to.bin, M, N, i,model1,model2)


####################################################################################################
#Mapping counts summed for each allele in each sample
####################################################################################################

mapcounts <- 1:12 ; dim(mapcounts) <- c(12,1)
mapcounts[,1]=0

#mel-mel

i=1
data1<-read.delim("zhrXz30_M_r1.genes.txt",header = T)
mapcounts[i,1]=sum(data1[,2])
i=i+1
mapcounts[i,1]=sum(data1[,3])
i=i+1
data1<-read.delim("zhrXz30_M_r2.genes.txt",header = T)
mapcounts[i,1]=sum(data1[,2])
i=i+1
mapcounts[i,1]=sum(data1[,3])
i=i+1
data1<-read.delim("zhrXz30_M_r3.genes.txt",header = T)
mapcounts[i,1]=sum(data1[,2])
i=i+1
mapcounts[i,1]=sum(data1[,3])
i=i+1

data1<-read.delim("z30Xzhr_M_r1.genes.txt",header = T)
mapcounts[i,1]=sum(data1[,2])
i=i+1
mapcounts[i,1]=sum(data1[,3])
i=i+1
data1<-read.delim("z30Xzhr_M_r2.genes.txt",header = T)
mapcounts[i,1]=sum(data1[,2])
i=i+1
mapcounts[i,1]=sum(data1[,3])
i=i+1
data1<-read.delim("z30Xzhr_M_r3.genes.txt",header = T)
mapcounts[i,1]=sum(data1[,2])
i=i+1
mapcounts[i,1]=sum(data1[,3])
i=i+1

mapcounts<-as.vector(mapcounts)

rm(i, data1)


####################################################################################################
#Construct data file
####################################################################################################

data1<-read.delim("zhrXz30_M_r1.genes.txt",header = T)
data2<-read.delim("zhrXz30_M_r2.genes.txt",header = T)
data3<-read.delim("zhrXz30_M_r3.genes.txt",header = T)


data4<-read.delim("z30Xzhr_M_r1.genes.txt",header = T)
data5<-read.delim("z30Xzhr_M_r2.genes.txt",header = T)
data6<-read.delim("z30Xzhr_M_r3.genes.txt",header = T)


melmel<-cbind(data1[,2],data2[,2],data3[,2],data1[,3],data2[,3],data3[,3],data4[,2],data5[,2],data6[,2],data4[,3],data5[,3],data6[,3])
names <- data1[,1]
melmel <- cbind(melmel,as.data.frame(names))

#melmel<-subset(melmel, ((melmel[,1]+melmel[,4]) >= 20) & ((melmel[,2]+melmel[,5]) >= 20) & ((melmel[,3]+melmel[,6]) >= 20) & ((melmel[,7]+melmel[,10]) >= 20) & ((melmel[,8]+melmel[,11]) >= 20) & ((melmel[,9]+melmel[,12]) >= 20))

list <- read.table("genes.txt",header=TRUE)
temp <- merge(melmel,list,by.x="names",by.y="gene")

genes <- temp[,1]
melmel <- temp[,2:13]
names(melmel)<-NULL

#melmel<-as.numeric(melmel)
melmel<-as.matrix(melmel)
storage.mode(melmel) <- "double"

rm(data1,data2,data3,data4,data5,data6)

####################################################################################################
#baySeq
####################################################################################################

library(baySeq)
library(snow)

#cl<-makeCluster(12)
cl <- makeMPIcluster(12)

#mel-mel

#data
data<-as.matrix(melmel)

#define replicate structure
replicates <-c(1,1,1,2,2,2,3,3,3,4,4,4)

#define hypotheses (list of groupings--get from model generator)
groups <- models

#define library sizes (vector of library sizes)
libsizes <- as.numeric(mapcounts)

#construct 'countData' object
CD <- new("countData", data = data, replicates = replicates, libsizes = libsizes, groups = groups)
CD@annotation <- as.data.frame(genes)

#estimate priors (negative binomial distributional assumption with ".NB")
CD.prior <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", cl = cl)
CD.post <- getLikelihoods.NB(CD.prior, pET = 'BIC', cl = cl)

#CD.prior <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", cl = NULL)
#CD.post <- getLikelihoods.NB(CD.prior, pET = 'BIC', cl = NULL)

#stopCluster(cl)

save(CD.post, file = "mel_mel_H1_H2_M_posteriors.bayseq")
save(CD.prior, file = "mel_mel_H1_H2_M_priors.bayseq")

write.table(cbind(CD.post@annotation,CD.post@posteriors), file = "mel_mel_H1_H2_M_posteriors.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

proc.time()-ptm
