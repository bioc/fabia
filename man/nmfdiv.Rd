\name{nmfdiv}
\alias{nmfdiv}
\title{Non-negative Matrix Factorization: Kullback-Leibler Divergence}
\description{

  \code{nmfdiv}: \R implementation of  \code{nmfdiv}.
}
\usage{

nmfdiv(X,p=5,cyc=100)


}
\arguments{
  \item{X}{the data matrix.}
  \item{p}{number of hidden factors = number of biclusters; default = 5.}
  \item{cyc}{maximal number of iterations; default = 100.}
}
\details{

  Non-negative Matrix Factorization represents positive matrix
  \eqn{X} by positive matrices \eqn{L} and \eqn{Z}.

   Objective for reconstruction is Kullback-Leibler divergence.

  Essentially the model is the sum of
  outer products of vectors:
   \deqn{X  =  \sum_{i=1}^{p} \lambda_i  z_i^T}
  where the number of summands  \eqn{p}
  is the number of biclusters.
  The matrix factorization is
  \deqn{X  = L  Z }

   Here \eqn{\lambda_i} are from \eqn{R^n}, \eqn{z_i} from
   \eqn{R^l}, \eqn{L} from \eqn{R^{n \times p}},
   \eqn{Z} from \eqn{R^{p \times l}}, and \eqn{X}
   from \eqn{R^{n \times l}}.
   

	 The model selection is performed
	 according  to  D. D. Lee and H. S. Seung, 1999, 2001.



  The code is implemented in \R.

}
\value{
 \item{}{object of the class \code{Factorization}. Containing
  \code{LZ} (estimated noise free data \eqn{L  Z}),
  \code{L} (loading \eqn{L}),
  \code{Z} (factors  \eqn{Z}),
  \code{U} (noise \eqn{X-LZ}),
  \code{X} (scaled data \eqn{X}).
}
}
\seealso{
\code{\link{fabia}},
\code{\link{fabias}},
\code{\link{fabiap}},
\code{\link{fabi}},
\code{\link{fabiasp}},
\code{\link{mfsc}},
\code{\link{nmfdiv}},
\code{\link{nmfeu}},
\code{\link{nmfsc}},
\code{\link{extractPlot}},
\code{\link{extractBic}},
\code{\link{plotBicluster}},
\code{\link{Factorization}},
\code{\link{projFuncPos}},
\code{\link{projFunc}},
\code{\link{estimateMode}},
\code{\link{makeFabiaData}},
\code{\link{makeFabiaDataBlocks}},
\code{\link{makeFabiaDataPos}},
\code{\link{makeFabiaDataBlocksPos}},
\code{\link{matrixImagePlot}},
\code{\link{fabiaDemo}},
\code{\link{fabiaVersion}}
}
\author{Sepp Hochreiter}
\examples{


#---------------
# TEST
#---------------

dat <- makeFabiaDataBlocks(n = 100,l= 50,p = 3,f1 = 5,f2 = 5,
  of1 = 5,of2 = 10,sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,
  sd_z = 1.0,sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

X <- dat[[1]]
Y <- dat[[2]]
X <- abs(X)


resEx <- nmfdiv(X,3)


\dontrun{
#---------------
# DEMO
#---------------

dat <- makeFabiaDataBlocks(n = 1000,l= 100,p = 10,f1 = 5,f2 = 5,
  of1 = 5,of2 = 10,sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,
  sd_z = 1.0,sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

X <- dat[[1]]
Y <- dat[[2]]
X <- abs(X)


resToy <- nmfdiv(X,13)

extractPlot(resToy,ti="NMFDIV",Y=Y)

}
}
\references{

  D. D. Lee and H. S. Seung,
  \sQuote{Algorithms for non-negative matrix factorization},
  In Advances in Neural Information Processing Systems 13, 556-562,
  2001.

  D. D. Lee and H. S. Seung,
  \sQuote{Learning the parts of objects by non-negative matrix
  factorization}, Nature, 401(6755):788-791, 1999.
}
\keyword{methods}
\concept{sparse coding}
\concept{non-negative matrix factorization}
