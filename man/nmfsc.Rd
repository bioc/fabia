\name{nmfsc}
\alias{nmfsc}
\title{Non-negative Sparse Matrix Factorization}
\description{

  \code{nmfsc}: \R implementation of  \code{nmfsc}.
}
\usage{

nmfsc(X,p=5,cyc=100,sL=0.6,sZ=0.6)


}
\arguments{
  \item{X}{the data matrix.}
  \item{p}{number of hidden factors = number of biclusters; default = 5.}
  \item{cyc}{maximal number of iterations; default = 100.}
  \item{sL}{sparseness loadings; default = 0.6.}
  \item{sZ}{sparseness factors; default = 0.6.}
}
\details{


     Non-negative Matrix Factorization represents positive matrix
  \eqn{X} by positive matrices \eqn{L} and \eqn{Z} that are sparse.

  Objective for reconstruction is Euclidean distance and sparseness
  constraints.



  Essentially the model is the sum of
  outer products of vectors:
   \deqn{X  =  \sum_{i=1}^{p} \lambda_i  z_i^T}
  where the number of summands  \eqn{p}
  is the number of biclusters.
  The matrix factorization is
  \deqn{X  =  L  Z }

   Here \eqn{\lambda_i} are from \eqn{R^n}, \eqn{z_i} from
   \eqn{R^l}, \eqn{L} from \eqn{R^{n \times p}},
   \eqn{Z} from \eqn{R^{p \times l}}, and \eqn{X}
   from \eqn{R^{n \times l}}.
   


  If the nonzero components of
  the sparse vectors are grouped together then the outer product results in
  a matrix with a nonzero block and zeros elsewhere.



	 The model selection is performed by a constraint optimization
	 according
  to  Hoyer, 2004. The Euclidean distance (the Frobenius norm) is
  minimized  subject to sparseness and non-negativity constraints.

  Model selection is done by gradient descent on the Euclidean
  objective and thereafter projection of single vectors of
  \eqn{L} and single vectors of \eqn{Z}
  to fulfill the sparseness and non-negativity constraints.

  The projection minimize the Euclidean distance
  to the original vector  given an
  \eqn{l_1}-norm and an \eqn{l_2}-norm and enforcing non-negativity.

  The projection is a convex quadratic problem which is solved
  iteratively where at each iteration at least one component is set to
  zero. Instead of the  \eqn{l_1}-norm a sparseness measurement is used
  which relates the \eqn{l_1}-norm to the \eqn{l_2}-norm.

  The code is implemented in \R.

}
\value{
 \item{}{object of the class \code{Factorization}. Containing
  \code{LZ} (estimated noise free data \eqn{L  Z}),
  \code{L} (loadings \eqn{L}),
  \code{Z} (factors \eqn{Z}),
  \code{U} (noise \eqn{X-LZ}),
  \code{X} (data \eqn{X}).
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
\code{\link{plot}},
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
\code{\link{summary}},
\code{\link{show}},
\code{\link{showSelected}},
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


resEx <- nmfsc(X,3,30,0.6,0.6)


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


resToy <- nmfsc(X,13,100,0.6,0.6)

extractPlot(resToy,ti="NMFSC",Y=Y)

}
}
\references{
  Patrik O. Hoyer,
  \sQuote{Non-negative Matrix Factorization with Sparseness Constraints},
  Journal of Machine Learning Research 5:1457-1469, 2004.

  D. D. Lee and H. S. Seung,
  \sQuote{Algorithms for non-negative matrix factorization},
  In Advances in Neural Information Processing Systems 13, 556-562,
  2001.
}
\keyword{methods}
\keyword{cluster}
\concept{sparse coding}
\concept{non-negative matrix factorization}
