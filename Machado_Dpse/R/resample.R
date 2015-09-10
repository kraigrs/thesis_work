########################################################################
# this part can be disregarded, was my attempt at creating a mock dataset

genes <- 20000;

PUF <- as.numeric(runif(genes)<0.10);
ALG <- as.numeric(runif(genes)<0.01);

data <- cbind(PUF,ALG);

nrow(subset(data,PUF==1&ALG==0));
nrow(subset(data,PUF==0&ALG==1));
observe <- nrow(subset(data,PUF==1&ALG==1));
nrow(subset(data,PUF==1&ALG==1));
nrow(subset(data,PUF==0&ALG==0));

# the data then look like this:

# gene     PUF     ALG
# 1        0       1
# 2        1       1
# 3        1       0
# 4        1       1
# ...
# n        0       1
#
# in a toy example of 20000 genes: 1908 bound by PUF only
#                                   170 bound by ALG only
#                                    16 bound by both PUF and ALG
#                                 17906 bound by neither
#         p-value = 0.7191
#
########################################################################

resample <- function(x,y)
{
	sample1 <- sample(x);
	sample2 <- sample(y);
	temp <- cbind(sample1,sample2);
	val <- nrow(subset(temp,sample1==1&sample2==1));
	return(val);
}

n <- 10000; # permutations to perform

null <- replicate(n,resample(PUF,ALG)); # creates null distribution of overlapping set
                                        # to be tested against your observed value

pval <- sum(null >= observe)/n; # number of values in the null set that are >= observed overlap
                                # divided by number of permutations