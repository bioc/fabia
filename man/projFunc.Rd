\name{projFunc}
\alias{projFunc}
\title{Projection of a Vector to a Sparse Vector}
\description{

  \code{projFunc}: \R implementation of  \code{projFunc}.
}
\usage{

projFunc(s, k1, k2)

}
\arguments{
  \item{s}{data vector.}
  \item{k1}{sparseness, l1 norm constraint.}
  \item{k2}{l2 norm constraint.}
}
\details{

  The projection is done according to Hoyer, 2004: given an
  \eqn{l_1}-norm and an \eqn{l_2}-norm minimize the Euclidean distance
  to the original vector.
  The projection is a convex quadratic problem which is solved
  iteratively where at each iteration at least one component is set to
  zero.

  In the applications,
  instead of the  \eqn{l_1}-norm a sparseness measurement is used
  which relates the \eqn{l_1}-norm to the \eqn{l_2}-norm.


  Implementation in \R.

}
\value{
  \item{v}{sparse projected vector.}

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
# DEMO
#---------------

size <- 30
sparseness <- 0.7

s <- as.vector(rnorm(size))
sp <- sqrt(1.0*size)-(sqrt(1.0*size)-1.0)*sparseness

ss <- projFunc(s,k1=sp,k2=1)

s
ss

}
\references{
  Patrik O. Hoyer,
  \sQuote{Non-negative Matrix Factorization with Sparseness Constraints},
  Journal of Machine Learning Research 5:1457-1469, 2004.
}
\keyword{methods}
\concept{sparse coding}
\concept{sparse matrix factorization}
