\name{extractPlot}
\alias{extractPlot}
\title{Plotting of Biclustering Results}
\description{

  \code{extractPlot}: \R implementation of  \code{extractPlot}.
}
\usage{

extractPlot(fact,thresZ=0.5,ti="",thresL=NULL,Y=NULL,which=c(1,2,3,4,5,6))

}
\arguments{
  \item{fact}{object of the class \code{Factorization}.}
  \item{thresZ}{threshold for sample belonging to bicluster; default 0.5.}
  \item{thresL}{threshold for loading belonging to bicluster (estimated
    if not given).}
  \item{ti}{plot title; default "".}
  \item{Y}{noise free data matrix.}
  \item{which}{which plot is shown: 1=noise free data (if available),
    2=data, 3=reconstructed data, 4=error, 5=absolute factors,
    6=absolute loadings; default c(1,2,3,4,5,6).}

}
\details{

  Essentially the model is the sum of
  outer products of vectors:
   \deqn{X  =  \sum_{i=1}^{p} \lambda_i  z_i^T  +  U}
  where the number of summands  \eqn{p}
  is the number of biclusters.
  The matrix factorization is
  \deqn{X  =  L  Z  +  U}
  

   Here \eqn{\lambda_i} are from \eqn{R^n}, \eqn{z_i} from
   \eqn{R^l}, \eqn{L} from \eqn{R^{n \times p}},
   \eqn{Z} from \eqn{R^{p \times l}}, and \eqn{X}, \eqn{U}
   from \eqn{R^{n \times l}}.
   
   The hidden dimension \eqn{p} is used for kmeans clustering of
   \eqn{\lambda_i} and \eqn{z_i }.

   The \eqn{\lambda_i } and \eqn{z_i } are used to
   extract the bicluster \eqn{i}, where a threshold determines which
   observations and which samples belong the the bicluster.

   The method produces following plots depending what plots are chosen
   by the "which" variable:

   \dQuote{Y}: noise free data (if available),
   \dQuote{X}: data,
   \dQuote{LZ}: reconstructed data,
   \dQuote{LZ-X}: error,
   \dQuote{abs(Z)}: absolute factors,
   \dQuote{abs(L)}: absolute loadings.

  Implementation in \R.
}
\value{Returns corresponding plots}
\seealso{
\code{\link{fabia}},
\code{\link{fabias}},
\code{\link{fabiap}},
\code{\link{fabi}},
\code{\link{fabiasp}},
\code{\link{spfabia}},
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


resEx <- fabia(X,3,0.1,20)

extractPlot(resEx,ti="FABIA",Y=Y)



\dontrun{
#---------------
# DEMO1
#---------------

dat <- makeFabiaDataBlocks(n = 1000,l= 100,p = 10,f1 = 5,f2 = 5,
  of1 = 5,of2 = 10,sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,
  sd_z = 1.0,sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

X <- dat[[1]]
Y <- dat[[2]]


resToy <- fabia(X,13,0.01,200)

extractPlot(resToy,ti="FABIA",Y=Y)

#---------------
# DEMO2
#---------------

avail <- require(fabiaData)

if (!avail) {
    message("")
    message("")
    message("#####################################################")
    message("Package 'fabiaData' is not available: please install.")
    message("#####################################################")
} else {

data(Breast_A)

X <- as.matrix(XBreast)

resBreast <- fabia(X,5,0.1,200)

extractPlot(resBreast,ti="FABIA Breast cancer(Veer)")

#sorting of predefined labels
CBreast%*%rBreast$pmZ
}

}
}
\keyword{hplot}
\concept{biclustering}
\concept{sparse coding}
\concept{sparse matrix factorization}
\concept{non-negative matrix factorization}
