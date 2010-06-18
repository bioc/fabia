\name{fabiasp}
\alias{fabiasp}
\title{Factor Analysis for Bicluster Acquisition: Sparseness Projection (FABIASP)}
\description{

  \code{fabiasp}: \R implementation of  \code{fabias}, therefore it is \bold{slow}.
}
\usage{

fabiasp(X,p=5,alpha=0.6,cyc=500,spz=0.5,center=2,norm=1)

}
\arguments{
  \item{X}{the data matrix.}
  \item{p}{number of hidden factors = number of biclusters; default = 5.}
  \item{alpha}{sparseness loadings (0.1 - 1.0); default = 0.1.}
  \item{cyc}{number of iterations; default = 500.}
  \item{spz}{sparseness factors (0.5 - 2.0); default = 0.5 (Laplace).}
  \item{center}{data centering: 1 (mean), 2
    (median), > 2 (mode), 0 (no); default = 2.}
  \item{norm}{data normalization:   1 (0.75-0.25 quantile), >1 (var=1), 0 (no); default = 1.}
}
\details{

  Biclusters are found by sparse factor analysis where \emph{both} the factors
  and the loadings are sparse.

  Essentially the model is the sum of
  outer products of vectors:
   \deqn{X  =  \sum_{i=1}^{p} \lambda_i  z_i^T  +  U}
  where the number of summands   \eqn{p}
  is the number of biclusters.
  The matrix factorization is
  \deqn{X  =  L  Z  +  U}

    Here \eqn{\lambda_i} are from \eqn{R^n}, \eqn{z_i} from
   \eqn{R^l}, \eqn{L} from \eqn{R^{n \times p}},
   \eqn{Z} from \eqn{R^{p \times l}}, and \eqn{X}, \eqn{U}
   from \eqn{R^{n \times l}}.
   

   If the nonzero components of
  the sparse vectors are grouped together then the outer product results in
  a matrix with a nonzero block and zeros elsewhere.


  The model selection is performed by a variational approach according
  to Girolami 2001 and Palmer et al. 2006.

  The prior has finite support, therefore after each
  update of the loadings they are projected to the finite support.
  The projection is done according to Hoyer, 2004: given an
  \eqn{l_1}-norm and an \eqn{l_2}-norm minimize the Euclidean distance
  to the original vector (currently the \eqn{l_2}-norm is fixed to 1).
  The projection is a convex quadratic problem which is solved
  iteratively where at each iteration at least one component is set to
  zero. Instead of the  \eqn{l_1}-norm a sparseness measurement is used
  which relates the \eqn{l_1}-norm to the \eqn{l_2}-norm.

   The code is implemented in \R, therefore it is \bold{slow}.

}
\value{
  \item{}{object of the class \code{Factorization}. Containing
  \code{LZ} (estimated noise free data \eqn{L  Z}),
  \code{L} (loadings \eqn{L}),
  \code{Z} (factors \eqn{Z}),
  \code{U} (noise \eqn{X-LZ}),
  \code{center} (centering vector),
  \code{scaleData} (scaling vector),
  \code{X} (centered and scaled data \eqn{X}),
  \code{Psi} (noise variance \eqn{\sigma}),
  \code{lapla} (variational parameter),
  \code{avini} (the information which the factor
    \eqn{z_{ij}} contains about \eqn{x_j} averaged over \eqn{j})
  \code{xavini} (the information which the factor
    \eqn{z_{j}} contains about \eqn{x_j})
  \code{ini} (for each \eqn{j} the information which the factor
    \eqn{z_{ij}} contains about \eqn{x_j}).
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



resEx <- fabiasp(X,3,0.6,50)


\dontrun{
#---------------
# DEMO1
#---------------

dat <- makeFabiaDataBlocks(n = 1000,l= 100,p = 10,f1 = 5,f2 = 5,
  of1 = 5,of2 = 10,sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,
  sd_z = 1.0,sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

X <- dat[[1]]
Y <- dat[[2]]


resToy <- fabiasp(X,13,0.6,200)

rToy <- extractPlot(resToy,"ti=FABIASP",Y=Y)

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

resBreast <- fabiasp(X,5,0.6,200)

rBreast <- extractPlot(resBreast,ti="FABIASP Breast cancer(Veer)")

#sorting of predefined labels
CBreast%*%rBreast$pmZ

}

#---------------
# DEMO3
#---------------


avail <- require(fabiaData)

if (!avail) {
    message("")
    message("")
    message("#####################################################")
    message("Package 'fabiaData' is not available: please install.")
    message("#####################################################")
} else {

data(Multi_A)

X <- as.matrix(XMulti)

resMulti <- fabiasp(X,5,0.6,200)

rMulti <- extractPlot(resMulti,"ti=FABIASP Multiple tissues(Su)")

#sorting of predefined labels
CMulti%*%rMulti$pmZ

}


#---------------
# DEMO4
#---------------


avail <- require(fabiaData)

if (!avail) {
    message("")
    message("")
    message("#####################################################")
    message("Package 'fabiaData' is not available: please install.")
    message("#####################################################")
} else {

data(DLBCL_B)

X <- as.matrix(XDLBCL)


resDLBCL <- fabiasp(X,5,0.6,200)

rDLBCL <- extractPlot(resDLBCL,ti="FABIASP Lymphoma(Rosenwald)")

#sorting of predefined labels
CDLBCL%*%rDLBCL$pmZ
}


}
}
\references{
  S. Hochreiter et al.,
  \sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
  Bioinformatics, doi:10.1093/bioinformatics/btq227, 2010.
  http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btq227
  
  Mark Girolami,
  \sQuote{A Variational Method for Learning Sparse and
  Overcomplete Representations},
  Neural Computation 13(11): 2517-2532, 2001.

 J. Palmer, D. Wipf, K. Kreutz-Delgado, B. Rao,
  \sQuote{Variational EM algorithms for
     non-Gaussian latent variable models},
     Advances in Neural Information Processing
     Systems 18, pp. 1059-1066, 2006.
     
  Patrik O. Hoyer,
  \sQuote{Non-negative Matrix Factorization with Sparseness Constraints},
  Journal of Machine Learning Research 5:1457-1469, 2004.
}
\keyword{methods}
\keyword{multivariate}
\keyword{cluster}
\concept{biclustering}
\concept{factor analysis}
\concept{sparse coding}
\concept{Laplace distribution}
\concept{EM algorithm}
\concept{nonnegative matrix factorization}
\concept{multivariate analysis}
\concept{latent variables}


