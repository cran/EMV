\name{knn}
\alias{knn}
\title{Estimate the missing values of a matrix}
\description{
  \code{knn} estimates the missing values of a matrix based on a k-th
  neighboors algorithm. Missing values can be either -Inf,Inf, NA, NaN.
}
\usage{
knn(m,k=max(dim(data)[1]*0.01,2),na.rm=TRUE,nan.rm=TRUE,inf.rm=TRUE)
}
\arguments{
  \item{m}{a numeric matrix that contains the missing values to be estimated}
  \item{k}{the number of neighboors (rows) to estimate the missing values}
  \item{na.rm}{a logical value indicating whether `NA' values should be estimated.}
  \item{nan.rm}{a logical value indicating whether `NaN' values should
    be estimated.}
  \item{inf.rm}{a logical value indicating whether `Inf' and '-Inf'
    values should be estimated.}
  
}
\details{
  Based on the Euclidian distance, the algorithm select the k-th nearest
  rows (that does not contain any missing values) to the one containing
  at least one missing value, based on the Euclidian distance. Then the
  missing values are replaced by the average of the neighboors. Note
  that if a row only contains missing values then the estimation is not
  possible.
}
\value{
  The function returns the same matrix, the missing values being
  replaced by their estimates.
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
## Estimate the missing values and infinite values
m2<-knn(m,k=10)
}
\keyword{NA}





