knn <-function(m,k=max(dim(m)[1]*0.01,2),na.rm=TRUE,nan.rm=TRUE,inf.rm=TRUE, correlation=FALSE, dist.bound=FALSE)
{
  if(is.matrix(m)==FALSE)
    stop(message="not a valid matrix object")
 
  n<-dim(m)[1]

### at least 2 to compute the mean 
  if(k<2)
    stop(message="k should be bigger than 1")

### at most n-1 neighboors
  if(k>n)
    k<-n-1
    
  nb.col<-dim(m)[2]
  nb.row<-dim(m)[1]

###code when linking to C
  missing.code<--9999999
  
  vector<-as.double(t(m))
  tmp<-vector
  
###replace the missing values by -9999999 (C code) 
  vector[is.finite(vector)==FALSE]<-missing.code
  
  result<-.C("knnc",
            vector=as.double(vector),
            nb.col=as.integer(nb.col),
            nb.row=as.integer(nb.row),
            k=as.integer(k),
            as.integer(correlation),
            distance=double(nb.row),
            as.double(dist.bound),PACKAGE="EMV")

  vector<-result$vector
### Still missing values if complete row of missing values
  vector[vector==missing.code]<-tmp[vector==missing.code]
  
### Remove the non-missing rows for the distances
  distance<-result$distance[result$distance!=missing.code & result$distance!=-missing.code]
  row<-(1:nb.row)[result$distance!=missing.code & result$distance!=-missing.code]
  distance<-cbind(row,distance)
  
  if(na.rm==FALSE)
    vector[is.na(tmp)==TRUE & is.nan(tmp)==FALSE]<-NA
  
  if(inf.rm==FALSE)
    {
      index<-is.finite(tmp)==FALSE & is.na(tmp)==FALSE & is.nan(tmp)==FALSE
      vector[index]<-tmp[index]
    }
  
  if(nan.rm==FALSE)
    vector[is.nan(tmp)==TRUE]<-NaN
  
  ##coerce vector back into the matrix
  newdata<-matrix(vector , nrow = nb.row, ncol = nb.col, byrow = TRUE)
  list(data=newdata,distance=distance)
}
