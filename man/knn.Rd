\name{knn}
\alias{knn}
\title{Estimate the missing values of a matrix}
\description{
  \code{knn} estimates the missing values of a matrix based on a k-th
  neighboors algorithm. Missing values can be either -Inf,Inf, NA, NaN.
}
\usage{
knn(m,k=max(dim(m)[1]*0.01,2),na.rm=TRUE,nan.rm=TRUE,inf.rm=TRUE,correlation=FALSE, dist.bound=FALSE)
}
\arguments{
  \item{m}{a numeric matrix that contains the missing values to be estimated}
  \item{k}{the number of neighboors (rows) to estimate the missing values}
  \item{na.rm}{a logical value indicating whether `NA' values should be estimated.}
  \item{nan.rm}{a logical value indicating whether `NaN' values should
    be estimated.}  
  \item{inf.rm}{a logical value indicating whether `Inf' and '-Inf'
    values should be estimated.}
  \item{correlation}{a logical value, if TRUE the selection of the
    neighboors is based on the sample correlation. The neighboors with
    the highest correlations are selected.}
  \item{dist.bound}{A bound for the distance, if correlation is FALSE,
    the algorithm will only use a neighboor if the Euclidean distance is
    less than dist.bound. If correlation is TRUE, the algorithm will
    only use a neighboor if the sample correlation is
    greater than dist.bound (in this case, between -1 and 1). If dist.bound=FALSE, all the neighboors are
    used.}
}
\details{
  Based on the Euclidian distance, the algorithm selects the k-th nearest
  rows (that do not contain any missing values) to the one containing
  at least one missing value, based on the Euclidian distance or the
  sample correlation. Then the missing values are replaced by the average of the neighboors. Note
  that if a row only contains missing values then the estimation is not
  possible.
}
\value{
  \item{data}{
  The  data matrix, the missing values being
  replaced by their estimates (when possible).}
  \item{distance}{
  The average of the neighboor's distances used for the estimation.}
}
\references{
  
  Missing Value estimation methods for DNA microarrays.
  O.Troyanskaya, M. Cantor, G. Sherlock, P. Brown, T. Hastie,
  R. Tibshirani, D. Botstein, & R. B. Altman.
  Bioinformatics 17(6):520-525, 2001. 
  }
\author{Raphael Gottardo \email{raph@lanl.gov}}
\seealso{
\code{\link{NA}},
\code{\link{NaN}},
\code{\link{Inf}}
}
\examples{
m<-matrix(rnorm(1000),100,10)
## Place some missing values NA
m[1:10,1]<-NA
m[50:52,10]<-NA
## Place some infinite values Inf, -Inf
m[1:10,3]<-(1/0)
m[70:73,10]<-(-1/0)
## Estimate the missing values and infinite values based on the Euclidean distance
m1<-knn(m,k=10)
## Estimate the missing values and infinite values based on the correlation distance
m2<-knn(m,k=10, correlation=TRUE)
}
\keyword{NA}





