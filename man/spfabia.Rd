\name{spfabia}
\alias{spfabia}
\title{Factor Analysis for Bicluster Acquisition: SPARSE FABIA}
\description{

  \code{spfabia}: C implementation of  \code{spfabia}.
}
\usage{

spfabia(X,p=5,alpha=0.1,cyc=500,spl=0,spz=0.5,non_negative=0,random=1.0,write_file=1,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0,samples=0,initL=0,iter=1,quant=0.001)

}
\arguments{
  \item{X}{the file name of the sparse matrix in sparse format.}
  \item{p}{number of hidden factors = number of biclusters; default = 5.}
  \item{alpha}{sparseness loadings (0 - 1.0); default = 0.1.}
  \item{cyc}{number of iterations; default = 500.}
  \item{spl}{sparseness prior loadings (0 - 2.0); default = 0 (Laplace).}
  \item{spz}{sparseness factors (0.5 - 2.0); default = 0.5 (Laplace).}
  \item{non_negative}{Non-negative factors and loadings if non_negative
    > 0; default = 0.}
  \item{random}{>0: random initialization of loadings in
    [0,random], <0: random initialization of loadings in
    [-random,random]; default = 1.0.}
\item{write_file}{>0: results are written to files (L in sparse format),
  default = 1.}
  \item{norm}{data normalization:  >0 (var=1),
    0 (no); default = 1.}
  \item{scale}{loading vectors are scaled in each iteration to the given
    variance. 0.0 indicates non scaling; default = 0.0.}
  \item{lap}{minimal value of the variational parameter; default = 1.0.} 
  \item{nL}{maximal number of biclusters at which a row element can participate; default
    = 0 (no limit).} 
  \item{lL}{maximal number of row elements per bicluster; default
    = 0 (no limit).} 
  \item{bL}{cycle at which the nL or lL maximum starts; default
    = 0 (start at the beginning).} 
  \item{samples}{vector of samples which should be included into the analysis; default
    = 0 (all samples)} 
  \item{initL}{vector of indices of the selected samples which are used to initialize L; default
    = 0 (random initialization).} 
  \item{iter}{number of iterations; default
    = 1.} 
  \item{quant}{qunatile of largest L values to remove in each iteration; default
    = 0.001.} 
}
\details{


  Version of fabia for a sparse data matrix. The data matrix is
  directly scanned by the C-code and must be in sparse matrix format.

  Sparse matrix format:
  *first line: numer of columns (the samples).
  *second line: number of rows (the features).
  *following lines: for each sample (column) three lines with

  I) number of nonzero row elements

  II) indices of the nonzero row elements

  III) values of the nonzero row elements

  
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
   
   If the nonzero components of
  the sparse vectors are grouped together then the outer product results in
  a matrix with a nonzero block and zeros elsewhere.


  The model selection is performed by a variational approach according
  to Girolami 2001 and Palmer et al. 2006.

  We included a prior on the parameters and minimize
  a lower bound on the posterior
  of the parameters given the data.
  The update of the loadings includes an additive term which pushes the
  loadings toward zero (Gaussian prior leads to an multiplicative factor).

  The code is implemented in C.

}
\value{
  \item{}{object of the class \code{Factorization}. Containing
  \code{L} (loadings  \eqn{L}),
  \code{Z} (factors  \eqn{Z}),
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
\code{\link{spfabia}},
\code{\link{readSamplesSpfabia}},
\code{\link{readSpfabiaResult}},
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

samples <- 20
features <- 200
sparseness <- 0.9
write(samples, file = "sparseFarmsTest.txt",ncolumns = features,append = FALSE, sep = " ")
write(features, file = "sparseFarmsTest.txt",ncolumns = features,append = TRUE, sep = " ")
for (i in 1:samples)
{
  ind <- which(runif(features)>sparseness)-1
  num <- length(ind)
  val <- abs(rnorm(num))
  write(num, file = "sparseFarmsTest.txt",ncolumns = features,append = TRUE, sep = " ")
  write(ind, file = "sparseFarmsTest.txt",ncolumns = features,append = TRUE, sep = " ")
  write(val, file = "sparseFarmsTest.txt",ncolumns = features,append = TRUE, sep = " ")
}

res <- spfabia("sparseFarmsTest",p=3,alpha=0.03,cyc=50,non_negative=1,write_file=0,norm=0)

unlink("sparseFarmsTest.txt")

plot(res,dim=c(1,2))
plot(res,dim=c(1,3))
plot(res,dim=c(2,3))



}


\references{
  S. Hochreiter et al.,
  \sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
  Bioinformatics 26(12):1520-1527, 2010.
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

}
\keyword{methods}
\keyword{multivariate}
\keyword{cluster}
\concept{biclustering}
\concept{factor analysis}
\concept{sparse coding}
\concept{Laplace distribution}
\concept{EM algorithm}
\concept{non-negative matrix factorization}
\concept{multivariate analysis}
\concept{latent variables}


