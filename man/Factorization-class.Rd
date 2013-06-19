\name{Factorization-class}
\docType{class}

% Class
\alias{Factorization-class}
\alias{Factorization-method}

% Constructor-like function:
\alias{Factorization}
\alias{Factorization,ANY-method}
\alias{Factorization,list-method,numeric-method,vector-method,matrix-method,ANY-method}

% Coercion methods:

% Accessor methods:
\alias{parameters}
\alias{parameters,Factorization-method}
\alias{parameters<-}
\alias{parameters<-,Factorization,list-method}
\alias{n}
\alias{n,Factorization-method}
\alias{n<-}
\alias{n<-,Factorization,numeric-method}
\alias{p1}
\alias{p1,Factorization-method}
\alias{p1<-}
\alias{p1<-,Factorization,numeric-method}
\alias{p2}
\alias{p2,Factorization-method}
\alias{p2<-}
\alias{p2<-,Factorization,numeric-method}
\alias{l}
\alias{l,Factorization-method}
\alias{l<-}
\alias{l<-,Factorization,numeric-method}
\alias{center}
\alias{center,Factorization-method}
\alias{center<-}
\alias{center<-,Factorization,numeric-method}
\alias{center<-,Factorization,vector-method}
\alias{scaleData}
\alias{scaleData,Factorization-method}
\alias{scaleData<-}
\alias{scaleData<-,Factorization,numeric-method}
\alias{scaleData<-,Factorization,vector-method}
\alias{X}
\alias{X,Factorization-method}
\alias{X<-}
\alias{X<-,Factorization,matrix-method}
\alias{L}
\alias{L,Factorization-method}
\alias{L<-}
\alias{L<-,Factorization,matrix-method}
\alias{Z}
\alias{Z,Factorization-method}
\alias{Z<-}
\alias{Z<-,Factorization,matrix-method}
\alias{M}
\alias{M,Factorization-method}
\alias{M<-}
\alias{M<-,Factorization,matrix-method}
\alias{LZ}
\alias{LZ,Factorization-method}
\alias{LZ<-}
\alias{LZ<-,Factorization,matrix-method}
\alias{U}
\alias{U,Factorization-method}
\alias{U<-}
\alias{U<-,Factorization,matrix-method}
\alias{avini}
\alias{avini,Factorization-method}
\alias{avini<-}
\alias{avini<-,Factorization,numeric-method}
\alias{avini<-,Factorization,vector-method}
\alias{xavini}
\alias{xavini,Factorization-method}
\alias{xavini<-}
\alias{xavini<-,Factorization,numeric-method}
\alias{xavini<-,Factorization,vector-method}
\alias{ini}
\alias{ini,Factorization-method}
\alias{ini<-}
\alias{ini<-,Factorization,matrix-method}
\alias{Psi}
\alias{Psi,Factorization-method}
\alias{Psi<-}
\alias{Psi<-,Factorization,numeric-method}
\alias{Psi<-,Factorization,vector-method}
\alias{lapla}
\alias{lapla,Factorization-method}
\alias{lapla<-}
\alias{lapla<-,Factorization,matrix-method}




% show method:
\alias{show,Factorization-method}
% summary method:
\alias{summary,Factorization-method}
% showSelected method:
\alias{showSelected}
\alias{showSelected,Factorization-method}
\alias{showSelected,Factorization,numeric-method}
% plot method:
\alias{plot,Factorization-method}
\alias{plot,Factorization,missing-method}
% show method:
\alias{show,Factorization-method}

\title{Factorization instances}

\description{\code{Factorization} is a class to store results of matrix
  factorization algorithms. It has been designed for biclustering but
  can be used for  "principal component analysis",
  "singular value decomposition", "independent component analysis",
  "factor analysis", and "non-negative matrix factorization".
}

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
           \item{\code{LZ}:}{ matrix L x M x Z.}
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

\section{Constructor}{
 Constructor of
 class Factorization.
 
\code{Factorization(parameters=list(),n=1,p1=1,p2=1,l=1,center=as.vector(1),scaleData=as.vector(1),X=as.matrix(1),L=as.matrix(1),Z=as.matrix(1),M=as.matrix(1),LZ=as.matrix(1),U=as.matrix(1),avini=as.vector(1),xavini=as.vector(1),ini=as.matrix(1),Psi=as.vector(1),lapla=as.matrix(1))}
}

\section{Accessors}{
  In the following \code{x} denotes a Factorization object.
   
 
  \describe{
    \item{}{
      \code{parameters(x)}, \code{parameters(x) <- value}:
      Returns or sets \code{parameters}, where the return value and
      \code{value} are both an instance of \code{list}. Parameters of the factorization
    method are stored in a list: ("method","number of cycles","sparseness
    weight","sparseness prior for loadings","sparseness prior for
    factors","number biclusters","projection sparseness loadings",
    "projection sparseness factors","initialization range","are loadings
    rescaled after each iterations","normalization = scaling of
    rows","centering method of rows","parameter for method").
    }

    \item{}{
      \code{n(x)}, \code{n(x) <- value}:
      Returns or sets \code{n}, where the return value and
      \code{value} are both an instance of \code{numeric}. Number of rows, left dimension.
    }

    \item{}{
      \code{p1(x)}, \code{p1(x) <- value}:
      Returns or sets \code{p1}, where the return value and
      \code{value} are both an instance of \code{numeric}. Right dimension of left matrix
    }

    \item{}{
      \code{p2(x)}, \code{p2(x) <- value}:
      Returns or sets \code{p2}, where the return value and
      \code{value} are both an instance of \code{numeric}.
      Left dimension of right matrix.
    }

    \item{}{
      \code{l(x)}, \code{l(x) <- value}:
      Returns or sets \code{l}, where the return value and
      \code{value} are both an instance of \code{numeric}.
      Number of columns, right dimension.
    }

    \item{}{
      \code{center(x)}, \code{center(x) <- value}:
      Returns or sets \code{center}, where the return value and
      \code{value} are both an instance of \code{numeric}.
      Vector of the centers.
    }

    \item{}{
      \code{scaleData(x)}, \code{scaleData(x) <- value}:
      Returns or sets \code{scaleData}, where the return value and
      \code{value} are both an instance of \code{numeric}.
      Vector of the scaling factors.
    }

    \item{}{
      \code{X(x)}, \code{X(x) <- value}:
      Returns or sets \code{X}, where the return value and
      \code{value} are both an instance of \code{matrix}.
      Centered and scaled data matrix n x l.
    }

    \item{}{
      \code{L(x)}, \code{L(x) <- value}:
      Returns or sets \code{L}, where the return value and
      \code{value} are both an instance of \code{matrix}.
      Left matrix n x p1.
    }

    \item{}{
      \code{Z(x)}, \code{Z(x) <- value}:
      Returns or sets \code{Z}, where the return value and
      \code{value} are both an instance of \code{matrix}.
      Right matrix p2 x l.
    }

    \item{}{
      \code{M(x)}, \code{M(x) <- value}:
      Returns or sets \code{M}, where the return value and
      \code{value} are both an instance of \code{matrix}.
      Middle matrix p1 x p2.
    }

    \item{}{
      \code{LZ(x)}, \code{LZ(x) <- value}:
      Returns or sets \code{LZ}, where the return value and
      \code{value} are both an instance of \code{matrix}.
      Matrix L x M x Z.
    }

    \item{}{
      \code{U(x)}, \code{U(x) <- value}:
      Returns or sets \code{U}, where the return value and
      \code{value} are both an instance of \code{matrix}.
      Noise matrix.
    }

    \item{}{
      \code{avini(x)}, \code{avini(x) <- value}:
      Returns or sets \code{avini}, where the return value and
      \code{value} are both an instance of \code{numeric}.
      Information of each bicluster, vector
           of length p2.
    }

    \item{}{
      \code{xavini(x)}, \code{xavini(x) <- value}:
      Returns or sets \code{xavini}, where the return value and
      \code{value} are both an instance of \code{numeric}.
      Information extracted from each sample, vector
           of length l.
    }

    \item{}{
      \code{ini(x)}, \code{ini(x) <- value}:
      Returns or sets \code{ini}, where the return value and
      \code{value} are both an instance of \code{matrix}.
      Information of each bicluster in each
           sample, matrix p2 x l.
    }

    \item{}{
      \code{Psi(x)}, \code{Psi(x) <- value}:
      Returns or sets \code{Psi}, where the return value and
      \code{value} are both an instance of \code{numeric}.
      Noise variance per row, vector of length n.
    }
    
    \item{}{
      \code{lapla(x)}, \code{lapla(x) <- value}:
      Returns or sets \code{lapla}, where the return value and
      \code{value} are both an instance of \code{matrix}.
      Prior information for each sample,
           vector of length l.
    }


  }
}


\usage{


\S4method{plot}{Factorization}(x, Rm=NULL, Cm=NULL, dim = c(1, 2),
    zoom = rep(1, 2), col.group = NULL,
    colors = c("orange1", "red", rainbow(length(unique(col.group)),
               start=2/6, end=4/6)),
    col.areas = TRUE, col.symbols = c(1, rep(2, length(unique(col.group)))),
    sampleNames = TRUE, rot = rep(-1, length(dim)),
    labels = NULL, label.tol = 0.1, lab.size = 0.725, col.size = 10,
    row.size = 10, do.smoothScatter = FALSE, 
    do.plot = TRUE, ... )

\S4method{show}{Factorization}(object)

\S4method{showSelected}{Factorization}(object, which=c(1,2,3,4))

\S4method{summary}{Factorization}(object, ...) 

}


\section{Signatures}{
  \describe{


\item{plot}{\code{signature(x = "Factorization", y = "missing")}}{ Plot of a matrix factorization result}

\item{show}{\code{signature(object = "Factorization")}}{ Display statistics of a matrix factorization result}

\item{showSelected}{\code{signature(object = "Factorization", which =
  "numeric")}}{ Display
  particular statistics of a matrix factorization result}

\item{summary}{\code{signature(object = "Factorization")}}{ Summary of matrix factorization result}


}
}

\arguments{
  \item{}{\bold{PLOT:}}
  \item{x}{object of the class \code{Factorization}.}
  \item{Rm}{row weighting vector. If \code{NULL}, it defaults to
    \code{rep(1,nrow(L(x)))}.}
  \item{Cm}{column weighting vector. If \code{NULL}, it defaults to
    \code{rep(1,ncol(Z(x)))}.}
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
    (\code{label.tol<=1}) of rows or the number of rows
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
  \item{}{\bold{SHOW:}}
   \item{object}{An instance of \code{\link[fabia]{Factorization-class}}
   .}
  \item{}{\bold{SHOWSELECTED:}}
   \item{}{see \code{object} at \code{show}.}
   \item{which}{used to provide a
     list of which plots should be generated: 1=the
     information content of biclusters, 2=the information content of
     samples, 3=the loadings per bicluster, 4=the factors per bicluster,
     default c(1,2,3,4).}
  \item{}{\bold{SUMMARY:}}
   \item{}{see \code{object} at \code{show}.}
   \item{}{... further arguments.}

}


\value{
  \item{}{\bold{FACTORIZATION:}}
   \item{}{An instance of \code{\link[fabia]{Factorization-class}}
   .}
   \item{}{\bold{PLOT:}}
   \item{Rows}{a list with the X and Y coordinates of the rows and
    an indication \code{Select} of whether the row was selected
    according to \code{label.tol}.}
  \item{Columns}{a list with the X and Y coordinates of the columns.}
  \item{}{\bold{SHOW:}}
  \item{}{no value.}
  \item{}{\bold{SHOWSELECTED:}}
  \item{}{no value.}
  \item{}{\bold{SUMMARY:}}
  \item{}{no value.}
}



\section{Functions that return objects of this class}{
  Factorization objects are returned by \code{fabia}, \code{fabias}, \code{fabiap},
  \code{fabiasp},  \code{mfsc}, \code{nmfsc},
  \code{nmfdiv}, and \code{nmfeu}.} 

\section{Extension to store results of other methods}{
  
  The class \code{Factorization} may contain the result of different matrix factorization
  methods. The methods may be generative or not.

  Methods my be "singular value decomposition" (M contains singular values
  as well as avini, L and Z are orthonormal matrices),
  "independent component analysis" (Z contains the projection/sources, L is the
  mixing matrix, M is unity), "factor analysis" (Z contains factors, L the loadings,
  M is unity, U the noise, Psi the noise covariance, lapla is a
  variational parameter for non-Gaussian factors, avini and ini are the
  information the factors convey about the observations).      
  
 
  }

  \details{

  \describe{
\item{\bold{Plot}}{Produces a biplot of a matrix factorization result stored in an
  instance of the Factorization class.}

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
   

\item{\bold{Show}}{Statistics of a matrix factorization result stored in
  an instance of the Factorization class.}

  This function supplies statistics on a matrix factorization result
  which is stored as an instance of
  \code{\link[fabia]{Factorization-class}}.

  The following is plotted:
  \enumerate{
    \item{}{the information content of biclusters.}
    \item{}{the information content of samples.}
    \item{}{the loadings per bicluster.}
    \item{}{the factors per bicluster.}
    }

\item{\bold{ShowSelected}}{
  Lists selected statistics of a matrix factorization
  result stored in an instance of the Factorization class.}


  This function supplies selected statistics on a matrix factorization result
  which is stored as an instance of
  \code{\link[fabia]{Factorization-class}}.

  The following is plotted depending on the display selection
  variable \code{which}:
  \enumerate{
    \item{}{the information content of biclusters.}
    \item{}{the information content of samples.}
    \item{}{the loadings per bicluster.}
    \item{}{the factors per bicluster.}
    }
  



\item{\bold{Summary}}{ Summary of matrix factorization result stored in an instance
  of the Factorization class.}

  This function gives information on a matrix factorization result
  which is stored as an instance of
  \code{\link[fabia]{Factorization-class}}.


  The summary consists of following items:
  \enumerate{
    \item{}{the number or rows and columns of the original matrix.}

  \item{}{the number of clusters for rows and columns is given.}

  \item{}{for the row cluster the information content is given.}

  \item{}{for each column its information is given.}

  \item{}{for each column cluster a summary is given.}
  
  \item{}{for each row cluster a summary is given.}

  }

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

###################
# TEST
###################


#------------------
#   PLOT
#------------------



n=200
l=100
p=4

dat <- makeFabiaDataBlocks(n = n,l= l,p = p,f1 = 5,f2 = 5,
  of1 = 5,of2 = 10,sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,
  sd_z = 1.0,sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

X <- dat[[1]]
ZC <- dat[[3]]
LC <- dat[[4]]


resEx <- fabia(X,p,0.01,400)


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

colnames(X(resEx)) <- clab

rownames(X(resEx)) <- llab


plot(resEx,dim=c(1,2),label.tol=0.1,col.group = groups,lab.size=0.6)
plot(resEx,dim=c(1,3),label.tol=0.1,col.group = groups,lab.size=0.6)
plot(resEx,dim=c(2,3),label.tol=0.1,col.group = groups,lab.size=0.6)



#------------------
#   SHOW
#------------------


dat <- makeFabiaDataBlocks(n = 100,l= 50,p = 3,f1 = 5,f2 = 5,
  of1 = 5,of2 = 10,sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,
  sd_z = 1.0,sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

X <- dat[[1]]


resEx <- fabia(X,3,0.01,100)

show(resEx)



#------------------
# SHOWSELECTED
#------------------

dat <- makeFabiaDataBlocks(n = 100,l= 50,p = 3,f1 = 5,f2 = 5,
  of1 = 5,of2 = 10,sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,
  sd_z = 1.0,sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

X <- dat[[1]]


resEx <- fabia(X,3,0.01,100)

showSelected(resEx,which=1)
showSelected(resEx,which=2)



#------------------
# SUMMARY
#------------------

dat <- makeFabiaDataBlocks(n = 100,l= 50,p = 3,f1 = 5,f2 = 5,
  of1 = 5,of2 = 10,sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,
  sd_z = 1.0,sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

X <- dat[[1]]


resEx <- fabia(X,3,0.01,100)

summary(resEx)

}


\keyword{classes}
\keyword{manip}
\keyword{methods}
\keyword{multivariate}
\keyword{hplot}
\concept{biclustering}
