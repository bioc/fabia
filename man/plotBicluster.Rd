\name{plotBicluster}
\alias{plotBicluster}
\title{Plotting of a bicluster}
\description{

  \code{plotBicluster}: \R implementation of  \code{plotBicluster}.
}
\usage{
plotBicluster(r,p,opp=FALSE,zlim=NULL,title=NULL,which=c(1, 2))
}
\arguments{
  \item{r}{the result of extract_bic.}
  \item{p}{the bicluster to plot. }
  \item{opp}{plot opposite bicluster, default = FALSE.}
  \item{zlim}{vector containing a low and high value to use for the color scale.}
  \item{title}{title of the plot.}
  \item{which}{which plots are shown: 1=data matrix with
    bicluster on upper left, 2=plot of the bicluster; default c(1, 2).}
}
\details{
  
  One bicluster is visualized by two plots. The variable
  "which" indicates which plots should be shown.

   Plot1 (which=1):
  The data matrix is sorted such that the bicluster appear at the
  upper left corner.
  The bicluster is marked by a rectangle.

  Plot2 (which=2):
  Only the bicluster is plotted.
  
   Implementation in \R.
}
\value{Plotting of a bicluster}
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


resEx <- fabia(X,3,0.01,20)

rEx <- extractBic(resEx)

plotBicluster(rEx,p=1)


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



rToy <- extractBic(resToy)

plotBicluster(rToy,p=1)

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

rBreast <- extractBic(resBreast)

plotBicluster(rBreast,p=1)

}

}
}
\keyword{hplot}
\concept{biclustering}

