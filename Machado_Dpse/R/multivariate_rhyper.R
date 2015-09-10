library(BiasedUrn)

m1 <- rpois(500,lambda=10);
m2 <- rpois(1000,lambda=100);
m3 <- rpois(2000,lambda=500);
m4 <- rpois(3000,lambda=1000);

m <- c(m1,m2,m3,m4);

m <- rpois(100000,lambda=10);

rMFNCHypergeo(1, m, sum(m)/4, odds=rep(1,length(m)), precision = 1E-7);