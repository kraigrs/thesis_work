#x <- rpois(100,30);
#y <- rpois(100,35);

#xSort <- sort(x);
#ySort <- sort(y);

#plot(xSort,ySort);

#m <- cbind(x,y);

############
#
# quantNorm - this function takes a matrix whereby each column
# is a separate distribution, i.e. for n arrays of length p, the 
# matrix is p x n where each array is a column. It then performs 
# quantile normalization by taking the mean quantile and 
# substituting it as a value of the data item in the original
# dataset
#
# Bolstad et al. 2003 "A comparison of normalization methods
# for high density oligonucleotide array data based on 
# variance and bias"
#
############

quantNorm <- function(m)
{
  #declare variables
  sorted <- NULL;
  orders <- NULL;
  means  <- NULL;
  norm   <- NULL;
  
  for(j in 1:ncol(m))
  {
    temp <- sort(m[,j],index.return=TRUE);
    sorted <- cbind(sorted,temp$x);
    orders <- cbind(orders,temp$ix);  	
  }
  
  for(i in 1:nrow(sorted))
  {
    means <- c(means,mean(sorted[i,]));	
  }
  
  for(j in 1:ncol(orders))
  {
    temp <- cbind(orders[,j],means);
    norm <- cbind(norm,temp[order(temp[,1]),2]);
  }
  return(norm);  	
}







