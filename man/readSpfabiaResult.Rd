\name{readSpfabiaResult}
\alias{readSpfabiaResult}
\title{Factor Analysis for Bicluster Acquisition: Read Results of SpFabia}
\description{

  \code{readSpfabiaResult}: C implementation of  \code{readSpfabiaResult}.
}
\usage{

readSpfabiaResult(X)

}
\arguments{
  \item{X}{the file prefix name of the result files of spfabia.}
}
\details{

  Read the results of spfabia.

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


