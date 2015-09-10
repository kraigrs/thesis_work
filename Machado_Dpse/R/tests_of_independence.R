# tests of independence among rows and columns

library(Deducer); # G test
library(lsr); # http://en.wikipedia.org/wiki/Cram%C3%A9r's_V

# MOI

# carcass

mat <- matrix(c(866,1607,2467,864,311,1751,3127,1985,1166,373,632,1183),nrow=2,ncol=6,byrow=TRUE);
Gtest <- likelihood.test(mat,conservative=TRUE);
Gtest$statistic;
Gtest$p.value;
cramersV(mat);

# females

mat <- matrix(c(866,1607,2467,864,311,1751,4244,701,519,266,47,254),nrow=2,ncol=6,byrow=TRUE);
Gtest <- likelihood.test(mat,conservative=TRUE);
Gtest$statistic;
Gtest$p.value;
cramersV(mat);

# males

mat <- matrix(c(3127,1985,1166,373,632,1183,3429,2358,1064,439,384,866),nrow=2,ncol=6,byrow=TRUE);
Gtest <- likelihood.test(mat,conservative=TRUE);
Gtest$statistic;
Gtest$p.value;
cramersV(mat);

# gonads

mat <- matrix(c(4244,701,519,266,47,254,3429,2358,1064,439,384,866),nrow=2,ncol=6,byrow=TRUE);
Gtest <- likelihood.test(mat,conservative=TRUE);
Gtest$statistic;
Gtest$p.value;
cramersV(mat);



# regdiv

# carcass

mat <- matrix(c(337,1956,620,2966,1481,353,1228,422,2951,860),nrow=2,ncol=5,byrow=TRUE);
Gtest <- likelihood.test(mat,conservative=TRUE);
Gtest$statistic;
cramersV(mat);

# females

mat <- matrix(c(337,1956,620,2966,1481,592,293,154,3484,1264),nrow=2,ncol=5,byrow=TRUE);
Gtest <- likelihood.test(mat,conservative=TRUE);
Gtest$statistic;
cramersV(mat);

# males

mat <- matrix(c(353,1228,422,2951,860,988,316,281,2706,1071),nrow=2,ncol=5,byrow=TRUE);
Gtest <- likelihood.test(mat,conservative=TRUE);
Gtest$statistic;
cramersV(mat);

# gonads

mat <- matrix(c(592,293,154,3484,1264,988,316,281,2706,1071),nrow=2,ncol=5,byrow=TRUE);
Gtest <- likelihood.test(mat,conservative=TRUE);
Gtest$statistic;
cramersV(mat);
