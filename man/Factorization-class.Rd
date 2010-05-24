\name{Factorization-class}
\title{The Factorization Class}
\docType{class}
\alias{Factorization}
\alias{Factorization-class}
\description{Factorization is the class structure for results of matrix
  factorization. Especially it is designed for factor analysis used for
  biclustering.
  
  The \code{summary} method shows information about biclusters.
  The \code{show} method plots information about biclusters.
  The \code{plot} method produces biplots of biclusters.
}

\section{Objects from the Class}{
  Objects can be created by \code{fabia}, \code{fabias}, \code{fabiap},
  \code{fabiasp},  \code{mfsc}, \code{nmfsc},
  \code{nmfdiv}, and \code{nmfeu}.} 


\section{Slots}{
  Objects of class \code{Factorization} have the following slots:
  \describe{
    
    \item{\code{parameters}:}{Saves parameters of the factorization
    method in a list: ("method","number of cycles","sparseness
    weight","sparseness prior for loadings","sparseness prior for
    factors","number biclusters","projection sparseness loadings",
    "projection sparseness factors","initialization range","are loadings
    rescaled after each iterations","normalization = scaling of
    rows","centering method of rows","parameter for method").}
           \item{\code{n}:}{ number of rows, left dimension.}
           \item{\code{p1}:}{ right dimension of left matrix.}
           \item{\code{p2}:}{ left dimension of right matrix.}
           \item{\code{l}:}{ number of columns, right dimension.}
           \item{\code{center}:}{ vector of the centers.}
           \item{\code{scaleData}:}{ vector of the scaling factors.}
	   \item{\code{X}:}{ centered and scaled data matrix n x l.}
           \item{\code{L}:}{ left matrix n x p1.} 
           \item{\code{Z}:}{ right matrix p2 x l.}
           \item{\code{M}:}{ middle matrix p1 x p2.}
           \item{\code{LZ}:}{ matrix L x M x Z .}
           \item{\code{U}:}{ noise matrix.}
           \item{\code{avini}:}{  information of each bicluster, vector
           of length p2.} 
           \item{\code{xavini}:}{  information extracted from each sample, vector
           of length l.} 
           \item{\code{ini}:}{ information of each bicluster in each
           sample, matrix p2 x l.}
           \item{\code{Psi}:}{ noise variance per row, vector of length n.}
           \item{\code{lapla}:}{ prior information for each sample,
           vector of length l.}    
     }
  }
\section{Details}{
  
  This class contains the result of different matrix factorization
  methods. The methods may be generative or not.

  Methods my be "singular value decomposition" (M contains singular values
  as well as avini, L and Z are orthonormal matrices),
  "independent component analysis" (Z contains the projection/sources, L is the
  mixing matrix, M is unity), "factor analysis" (Z contains factors, L the loadings,
  M is unity, U the noise, Psi the noise covariance, lapla is a
  variational parameter for non-Gaussian factors, avini and ini are the
  information the factors convey about the observations).      
  
  
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
\keyword{classes}

