\name{estimateMode}
\alias{estimateMode}
\title{Estimation of the modes of the rows of a matrix}
\description{

  \code{estimateMode}: \R implementation of  \code{estimateMode}.
}
\usage{

estimateMode(X,maxiter=50,tol=0.001,alpha=0.1,a1=4.0,G1=FALSE)

}
\arguments{
  \item{X}{matrix of which the modes of the rows are estimated.}
  \item{maxiter}{maximal number of iterations; default = 50.}
  \item{tol}{tolerance for stopping; default = 0.001.}
  \item{alpha}{learning rate; default = 0.1.}
  \item{a1}{parameter of the width of the given distribution; default = 4.}
  \item{G1}{kind of distribution,  \code{TRUE}: G1,
   \code{FALSE}: G2; default = \code{FALSE}.}
}
\details{

  The mode is estimated by contrast functions G1
  \deqn{(1/a_1) * \ln (\cosh (a1*x))}
  or G2
  \deqn{ - (1/a_1)*\exp(-1/2 * x*x)}
  The estimation is performed by gradient descent initialized by the
  median.
  
  
  Implementation in \R.

}
\value{
  \item{u}{the vector of estimated modes.}
  \item{xu}{\eqn{X-u} the mode centered data.}

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
# DEMO
#---------------

dat <- makeFabiaDataBlocksPos(n = 100,l= 50,p = 10,f1 = 5,f2 = 5,
       of1 = 5,of2 = 10,sd_noise = 2.0,sd_z_noise = 0.2,mean_z = 2.0,
       sd_z = 1.0,sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

     X <- dat[[1]]

modes <- estimateMode(X)

modes$u - apply(X, 1, median)


}
\references{
A. Hyvaerinen,
\sQuote{Fast and Robust Fixed-Point Algorithms for Independent
Component Analysis},
IEEE Transactions on Neural Networks 10(3):626-634,
1999.

}
\keyword{methods}
\concept{mode estimation}

