\name{show-methods}
\docType{methods}
\alias{show}
\alias{show-methods}
\alias{show,Factorization-method}
\title{Plot Statistics of a Matrix Factorization}
\description{
  This function plots statistics on a matrix factorization
  which is stored as an instance of
  \code{\link[fabia]{Factorization-class}}.}
\section{Methods}{
\describe{
\item{\code{signature(object = "Factorization")}}{Display Statistics of a Matrix Factorization }
}}
\usage{\S4method{show}{Factorization}(object)}
\arguments{
   \item{object}{An instance of \code{\link[fabia]{Factorization-class}}
   .}
 }
\details{
  The following is plotted:
  \enumerate{
    \item{the information content of biclusters.}
    \item{the information content of samples.}
    \item{the loadings per bicluster.}
    \item{the factors per bicluster.}
    }
}
\value{
  no value.
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


resEx <- fabia(X,3,0.1,100,1.0,1.0)

show(resEx)

}
\keyword{manip}
\keyword{methods}

