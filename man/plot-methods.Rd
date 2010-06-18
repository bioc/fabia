\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,Factorization,missing-method}
\title{Plot of a Matrix Factorization}
\description{
  Produces a (biplot) of a Matrix Factorization result
}
\section{Methods}{
\describe{

\item{\code{signature(x = "Factorization", y = "missing")}}{ Plot of a Matrix Factorization }


}
}


\usage{
\S4method{plot}{Factorization,missing}(x, Rm=NULL, Cm=NULL, dim = c(1, 2),
    zoom = rep(1, 2), col.group = NULL,
    colors = c("orange1", "red", rainbow(length(unique(col.group)),
               start=2/6, end=4/6)),
    col.areas = TRUE, col.symbols = c(1, rep(2, length(unique(col.group)))),
    sampleNames = TRUE, rot = rep(-1, length(dim)),
    labels = NULL, label.tol = 0.1, lab.size = 0.725, col.size = 10,
    row.size = 10, do.smoothScatter = FALSE, 
    do.plot = TRUE, ... )
}
\arguments{
  \item{x}{object of the class \code{Factorization}.}
  \item{Rm}{row weighting vector. If \code{NULL}, it defaults to
    \code{rep(1,nrow(x@L))}.}
  \item{Cm}{column weighting vector. If \code{NULL}, it defaults to
    \code{rep(1,ncol(x@Z))}.}
  \item{dim}{optional principal factors that are plotted along the
    horizontal and vertical axis. Defaults to \code{c(1,2)}.}
  \item{zoom}{optional zoom factor for row and column items. Defaults to
    \code{c(1,1)}.}
  \item{col.group}{optional vector (character or numeric) indicating the
    different groupings of the columns. Defaults to 1.}
  \item{colors}{vector specifying the colors for the annotation of the plot;
    the first two elements concern the rows; the third till the last element
    concern the columns; the first element will be used to color the unlabeled
    rows; the second element for the labeled rows and the remaining elements to
    give different colors to different groups of columns. Defaults to
    \code{c("orange1", "red", rainbow(length(unique(col.group)),
    start=2/6, end=4/6))}.}
  \item{col.areas}{logical value indicating whether columns should be
    plotted as squares with areas proportional to their marginal mean
    and colors representing the different groups (\code{TRUE}), or with
    symbols representing the groupings and identical size
    (\code{FALSE}). Defaults to \code{TRUE}.}
  \item{col.symbols}{vector of symbols when \code{col.areas=FALSE}
    corresponds to the \code{pch} argument of the function \code{plot}.
    Defaults to \code{c(1, rep(2, length(unique(col.group))))}.}
  \item{sampleNames}{either a logical vector of length one or a character vector
    of length equal to the number of samples in the dataset. If a
    logical is provided, sample names will be displayed on the plot
    (\code{TRUE}; default) or not (\code{FALSE}); if a character vector
    is provided, the names provided will be used to label the samples
    instead of the default column names.}
  \item{rot}{rotation of plot. Defaults to \code{c(-1,-1)}.}
  \item{labels}{character vector to be used for labeling points on the graph;
    if \code{NULL} (default), the row names of \code{x} are used instead.}
  \item{label.tol}{numerical value specifying either the percentile
    (\code{label.tol<=1})of rows or the number of rows
    (\code{label.tol>1}) most distant from the plot-center (0,0) that
    are labeled and are plotted as circles with area proportional to the
    marginal means of the original data. Defaults to \code{1}.}
  \item{lab.size}{size of identifying labels for row- and column-items
    as \code{cex} parameter of the \code{text} function. Defaults to
    \code{0.725}.}
  \item{col.size}{size of the column symbols in mm. Defaults to \code{10}.}
  \item{row.size}{size of the row symbols in mm. Defaults to \code{10}.}
  \item{do.smoothScatter}{use smoothScatter or not instead of plotting
    individual points. Defaults to \code{FALSE}.}
  \item{do.plot}{produce a plot or not. Defaults to \code{TRUE}.}
  \item{...}{further arguments are passed on to \code{eqscaleplotLoc}
    which draws the canvas for the plot; useful for adding a \code{main}
    or a custom \code{sub}.}
}

\details{
  
The function \code{plot} is based on
the function \code{plot.mpm} in the \R package \code{mpm}
(Version: 1.0-16, Date: 2009-08-26, Title: Multivariate Projection
Methods, Maintainer: Tobias Verbeke <tobias.verbeke@openanalytics.be>,
Author: Luc Wouters <wouters_luc@telenet.be>).


 Biclusters are found by sparse factor analysis where \emph{both} the factors
  and the loadings are sparse.

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
     

   For noise free projection like independent component analysis
   we set the noise term to zero: \eqn{U=0}.
   
   
   The argument \code{label.tol} can be used to select the
   most informative rows, i.e. rows that are most distant from the
   center of the plot
   (smaller 1: percentage of rows, larger 1: number of rows). 
 
   Only these row-items are then labeled and represented as circles with
   their areas proportional to the row weighting.
   
   If the column-items are grouped these groups can be visualized by
   colors given  by \code{col.group}.
   
}
\value{
  \item{Rows}{a list with the X and Y coordinates of the rows and
    an indication \code{Select} of whether the row was selected
    according to \code{label.tol}.}
  \item{Columns}{a list with the X and Y coordinates of the columns.}
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

n=200
l=100
p=4

dat <- makeFabiaDataBlocks(n = n,l= l,p = p,f1 = 5,f2 = 5,
  of1 = 5,of2 = 10,sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,
  sd_z = 1.0,sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

X <- dat[[1]]
ZC <- dat[[3]]
LC <- dat[[4]]


resEx <- fabia(X,p,0.1,400)


gclab <- rep.int(0,l)
gllab <- rep.int(0,n)
clab <- as.character(1:l)
llab <- as.character(1:n)
for (i in 1:p){
 for (j in ZC[i]){
     clab[j] <- paste(as.character(i),"_",clab[j],sep="")
 }
 for (j in LC[i]){
     llab[j] <- paste(as.character(i),"_",llab[j],sep="")
 }
 gclab[unlist(ZC[i])] <- gclab[unlist(ZC[i])] + p^i
 gllab[unlist(LC[i])] <- gllab[unlist(LC[i])] + p^i
}


groups <- gclab

colnames(resEx@X) <- clab

rownames(resEx@X) <- llab


plot(resEx,dim=c(1,2),label.tol=0.1,col.group = groups,lab.size=0.6)
plot(resEx,dim=c(1,3),label.tol=0.1,col.group = groups,lab.size=0.6)
plot(resEx,dim=c(2,3),label.tol=0.1,col.group = groups,lab.size=0.6)

}

\keyword{methods}
\keyword{multivariate}
\keyword{hplot}
\concept{biclustering}

