\name{summary-methods}
\docType{methods}
\alias{summary}
\alias{summary-methods}
\alias{summary,Factorization-method}
\title{Summary of Matrix Factorization}
\description{
  This function gives information on a matrix factorization
  which is stored as an instance of
  \code{\link[fabia]{Factorization-class}}.
  }
\section{Methods}{
\describe{
\item{\code{signature(object = "Factorization")}}{ Summary of Matrix Factorization }
}}
\usage{\S4method{summary}{Factorization}(object, ...) }
\arguments{
   \item{object}{An instance of \code{\link[fabia]{Factorization-class}} .}
   \item{...}{further arguments.}
         }
\details{
  
  First, it gives the number or rows and columns of the original matrix.

  Then the number of clusters for rows and columns is given.

  Then for the row cluster the information content is given.

  Then for each column its information is given.

  Then for each column cluster a summary is given.
  
  Then for each row cluster a summary is given.
  
}
\value{ \item{}{no value.}
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


resEx <- fabia(X,3,0.1,100)

summary(resEx)

}
\keyword{manip}
\keyword{methods}

