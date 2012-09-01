\name{matrixImagePlot}
\alias{matrixImagePlot}
\title{Plotting of a Matrix}
\description{

  \code{matrixImagePlot}: \R implementation of  \code{myImagePlot}.
}
\usage{

matrixImagePlot(x,xLabels=NULL, yLabels=NULL, zlim=NULL, title=NULL)


}
\arguments{
  \item{x}{the matrix.}
  \item{xLabels}{vector of strings to label the columns
    (default "colnames(x)").}
  \item{yLabels}{vector of strings to label the rows (default "rownames(x)").}
  \item{zlim}{vector containing a low and high value to use for the color scale.}
  \item{title}{title of the plot.}

}
\details{

  Plotting a table of numbers as an image using \R.

  The color scale is based on the highest and lowest values in the
  matrix.

  The original \R code has been obtained by
 \url{http://www.phaget4.org/R/myImagePlot.R}
  and then has been modified.

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
\examples{


#---------------
# TEST
#---------------

dat <- makeFabiaDataBlocks(n = 100,l= 50,p = 3,f1 = 5,f2 = 5,
  of1 = 5,of2 = 10,sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,
  sd_z = 1.0,sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

X <- dat[[1]]
Y <- dat[[2]]

matrixImagePlot(Y)
dev.new()
matrixImagePlot(X)


\dontrun{
#---------------
# DEMO
#---------------

dat <- makeFabiaDataBlocks(n = 1000,l= 100,p = 10,f1 = 5,f2 = 5,
  of1 = 5,of2 = 10,sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,
  sd_z = 1.0,sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

X <- dat[[1]]
Y <- dat[[2]]
X <- X- rowMeans(X)
XX <- (1/ncol(X))*tcrossprod(X)
dXX <- 1/sqrt(diag(XX)+0.001*as.vector(rep(1,nrow(X))))
X <- dXX*X

matrixImagePlot(X)

}
}
\references{

 \url{http://www.phaget4.org/R/myImagePlot.R}
}
\keyword{hplot}
\concept{matrix plot}

