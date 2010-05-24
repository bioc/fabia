\name{showSelected-methods}
\docType{methods}
\alias{showSelected}
\alias{showSelected-methods}
\alias{showSelected,Factorization-method}
\title{Plot Selected Statistics of a Matrix Factorization}
\description{
  This function plots selected statistics on a matrix factorization
  which is stored as an instance of
  \code{\link[fabia]{Factorization-class}}.}
\section{Methods}{
  \describe{
\item{\code{ signature(object="Factorization")}}{Display Statistics of a Matrix Factorization }
}}
\usage{\S4method{showSelected}{Factorization}(object, which=c(1,2,3,4))}
\arguments{
   \item{object}{An instance of \code{\link[fabia]{Factorization-class}}.}
   \item{which}{used to provide a
     list of which plots should be generated: 1=the
     information content of biclusters, 2=the information content of
     samples, 3=the loadings per bicluster, 4=the factors per bicluster,
     default c(1,2,3,4).}
  }
\details{

  The following is plotted depending on the display selection
  variable \code{which}:
  \enumerate{
    \item{the information content of biclusters.}
    \item{the information content of samples.}
    \item{the loadings per bicluster.}
    \item{the factors per bicluster.}
    }
  
}
\value{
  \item{}{no value.}
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

showSelected(resEx,which=1)
showSelected(resEx,which=2)

}
\keyword{manip}
\keyword{methods}

