\name{extractBic}
\alias{extractBic}
\title{Extraction of Biclusters}
\description{

  \code{extractBic}: \R implementation of  \code{extractBic}.
}
\usage{

extractBic(fact,thresZ=0.5,thresL=NULL)

}
\arguments{
  \item{fact}{object of the class \code{Factorization}.}
  \item{thresZ}{threshold for sample belonging to bicluster; default 0.5.}
  \item{thresL}{threshold for loading belonging to bicluster (if not
    given it is estimated).}
}
\details{

  Essentially the model is the sum of
  outer products of vectors:
   \deqn{X  =  \sum_{i=1}^{p} \lambda_i   z_i^T  +  U}
  where the number of summands  \eqn{p}
  is the number of biclusters.
  The matrix factorization is
  \deqn{X  =  L  Z  +  U}

   Here \eqn{\lambda_i} are from \eqn{R^n}, \eqn{z_i} from
   \eqn{R^l}, \eqn{L} from \eqn{R^{n \times p}},
   \eqn{Z} from \eqn{R^{p \times l}}, and \eqn{X}, \eqn{U}
   from \eqn{R^{n \times l}}.


   \eqn{U } is the Gaussian noise with a diagonal covariance matrix
   which entries are given by \code{Psi}.
  
   The \eqn{Z } is locally approximated by a Gaussian with inverse
   variance given by  \code{lapla}.
   
   Using these values we can computer
   for each  \eqn{j} the variance \eqn{ z_i} given \eqn{ x_j}. Here

   \deqn{ x_j  =  L  z_j  +  u_j}

   This variance can be used to determine the information content
   of a bicluster.
   
   
   The \eqn{\lambda_i } and \eqn{z_i } are used to
   extract the bicluster \eqn{i}, where a threshold determines which
   observations and which samples belong the the bicluster.



   In \code{bic} the biclusters are extracted according to the
   largest absolute values of the component \eqn{i}, i.e.
   the largest values of \eqn{\lambda_i } and the
   largest values of \eqn{z_i }. The factors \eqn{ z_i}
   are normalized to variance 1.

   
   The components of \code{bic} are
   \code{binp}, \code{bixv},
   \code{bixn}, \code{biypv}, and \code{biypn}.
   

   \code{binp} give the size of the bicluster: number observations and
   number samples.
   \code{bixv} gives the values of the extracted
   observations that have absolute
   values above a threshold. They are sorted. 
   \code{bixn} gives  the extracted observation names (e.g. gene names).
   \code{biypv} gives the values of the extracted samples that have
   absolute values above a threshold.  They are sorted.
   \code{biypn} gives  the names of the extracted samples (e.g. sample names).

   In \code{bicopp} the opposite clusters to the biclusters are
   given. Opposite means that the negative pattern is present.
   
   The components of opposite clusters \code{bicopp} are
   \code{binn}, \code{bixv},
   \code{bixn}, \code{biypnv}, and \code{biynn}.
   
    \code{binp} give the size of the opposite bicluster: number observations and
   number samples.
   \code{bixv} gives the values of the extracted
   observations that have absolute
   values above a threshold. They are sorted. 
   \code{bixn} gives  the extracted observation names (e.g. gene names).
   \code{biynv} gives the values of the opposite extracted samples that have
   absolute values above a threshold.  They are sorted.
   \code{biynn} gives  the names of the opposite
   extracted samples (e.g. sample names).
   
   That means the samples are divided into two groups where one group
   shows large positive values and the other group has negative values
   with large absolute values. That means a observation pattern can be
   switched on or switched off relative to the average value.

   \code{numn} gives the indices of \code{bic} with components:
   \code{numng} = \code{bix} and \code{numnp}  =  \code{biypn}.
   
   \code{numn} gives the indices of \code{bicopp} with components:
   \code{numng} = \code{bix} and \code{numnn}  =  \code{biynn}.
   


  Implementation in \R.
}

\value{
  \item{bic}{extracted biclusters.}
  \item{numn}{indexes for the extracted biclusters.}
  \item{bicopp}{extracted opposite biclusters.}
  \item{numnopp}{indexes for the extracted opposite biclusters.}
  \item{X}{scaled and centered data matrix.}
  \item{np}{number of biclusters.}
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

rEx$bic[1,]
rEx$bic[2,]
rEx$bic[3,]


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

avini(resToy)

rToy$bic[1,]
rToy$bic[2,]
rToy$bic[3,]

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

avini(resBreast)

rBreast$bic[1,]
rBreast$bic[2,]
rBreast$bic[3,]
}

}
}
\concept{biclustering}
\concept{sparse coding}
\concept{sparse matrix factorization}
\concept{non-negative matrix factorization}
