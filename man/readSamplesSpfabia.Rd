\name{readSamplesSpfabia}
\alias{readSamplesSpfabia}
\title{Factor Analysis for Bicluster Acquisition: Read Sparse Matrix Samples}
\description{

  \code{readSamplesSpfabia}: C implementation of  \code{readSamplesSpfabia}.
}
\usage{

readSamplesSpfabia(X,samples=0)

}
\arguments{
  \item{X}{the file name of the sparse matrix in sparse format.}
  \item{samples}{vector of samples which should be read; default
    = 0 (all samples)} 
}
\details{


  The data matrix is
  directly scanned by the C-code and must be in sparse matrix format.

  Sparse matrix format:
  *first line: numer of columns (the samples).
  *second line: number of rows (the features).
  *following lines: for each sample (column) three lines with

  I) number of nonzero row elements

  II) indices of the nonzero row elements

  III) values of the nonzero row elements


  The code is implemented in C.

}
\value{
  \code{X} (data of the given samples),
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



}


\references{
  S. Hochreiter et al.,
  \sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
  Bioinformatics 26(12):1520-1527, 2010.
  http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btq227
  

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

