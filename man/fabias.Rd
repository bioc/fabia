\name{fabias}
\alias{fabias}
\title{Factor Analysis for Bicluster Acquisition: Sparseness Projection (FABIAS)}
\description{

  \code{fabias}: C implementation of  \code{fabias}.
}
\usage{

fabias(X,p=5,alpha=0.6,cyc=500,spz=0.5,non_negative=0,random=1.0,center=2,norm=1,lap=1.0,nL=0,lL=0,bL=0)

}
\arguments{
  \item{X}{the data matrix.}
  \item{p}{number of hidden factors = number of biclusters; default = 5.}
  \item{alpha}{sparseness loadings (0.1 - 1.0); default = 0.1.}
  \item{cyc}{number of iterations; default = 500.}
  \item{spz}{sparseness factors (0.5 - 2.0); default = 0.5 (Laplace).}
  \item{non_negative}{Non-negative factors and loadings if non_negative
    > 0; default = 0.}
  \item{random}{<=0: by SVD, >0: random initialization of loadings in
    [-random,random]; default = 1.0.}
  \item{center}{data centering: 1 (mean), 2
    (median), > 2 (mode), 0 (no); default = 2.}
  \item{norm}{data normalization:  1 (0.75-0.25 quantile), >1 (var=1), 0 (no); default = 1.}
  \item{lap}{minimal value of the variational parameter; default = 1.0.} 
  \item{nL}{maximal number of biclusters at which a row element can participate; default
    = 0 (no limit)} 
  \item{lL}{maximal number of row elements per bicluster; default
    = 0 (no limit)} 
  \item{bL}{cycle at which the nL or lL maximum starts; default
    = 0 (start at the beginning)} 
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


  The code is implemented in C.

}
\value{
  \item{}{object of the class \code{Factorization}. Containing
  \code{LZ} (estimated noise free data \eqn{L  Z}),
  \code{L} (loadings  \eqn{L}),
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
\code{\link{spfabia}},
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




resEx <- fabias(X,3,0.6,50)


\dontrun{
#-----------------
# DEMO1: Toy Data
#-----------------

n = 1000
l= 100
p = 10

dat <- makeFabiaDataBlocks(n = n,l= l,p = p,f1 = 5,f2 = 5,
  of1 = 5,of2 = 10,sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,
  sd_z = 1.0,sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

X <- dat[[1]]
Y <- dat[[2]]
ZC <- dat[[3]]
LC <- dat[[4]]

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

#### FABIAS

resToy2 <- fabias(X,13,0.6,400)

extractPlot(resToy2,ti="FABIAS",Y=Y)


raToy2 <- extractBic(resToy2)

if ((raToy2$bic[[1]][1]>1) && (raToy2$bic[[1]][2])>1) {
    plotBicluster(raToy2,1)
}
if ((raToy2$bic[[2]][1]>1) && (raToy2$bic[[2]][2])>1) {
    plotBicluster(raToy2,2)
}
if ((raToy2$bic[[3]][1]>1) && (raToy2$bic[[3]][2])>1) {
    plotBicluster(raToy2,3)
}
if ((raToy2$bic[[4]][1]>1) && (raToy2$bic[[4]][2])>1) {
    plotBicluster(raToy2,4)
}

colnames(X(resToy2)) <- clab

rownames(X(resToy2)) <- llab


plot(resToy2,dim=c(1,2),label.tol=0.1,col.group = groups,lab.size=0.6)
plot(resToy2,dim=c(1,3),label.tol=0.1,col.group = groups,lab.size=0.6)
plot(resToy2,dim=c(2,3),label.tol=0.1,col.group = groups,lab.size=0.6)


#------------------------------------------
# DEMO2: Laura van't Veer's gene expression  
#        data set for breast cancer 
#------------------------------------------

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

resBreast2 <- fabias(X,5,0.6,300)

extractPlot(resBreast2,ti="FABIAS Breast cancer(Veer)")


raBreast2 <- extractBic(resBreast2)

if ((raBreast2$bic[[1]][1]>1) && (raBreast2$bic[[1]][2])>1) {
    plotBicluster(raBreast2,1)
}
if ((raBreast2$bic[[2]][1]>1) && (raBreast2$bic[[2]][2])>1) {
    plotBicluster(raBreast2,2)
}
if ((raBreast2$bic[[3]][1]>1) && (raBreast2$bic[[3]][2])>1) {
    plotBicluster(raBreast2,3)
}
if ((raBreast2$bic[[4]][1]>1) && (raBreast2$bic[[4]][2])>1) {
    plotBicluster(raBreast2,4)
}

plot(resBreast2,dim=c(1,2),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast2,dim=c(1,3),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast2,dim=c(1,4),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast2,dim=c(1,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast2,dim=c(2,3),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast2,dim=c(2,4),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast2,dim=c(2,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast2,dim=c(3,4),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast2,dim=c(3,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast2,dim=c(4,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)

}

#-----------------------------------
# DEMO3: Su's multiple tissue types
#        gene expression data set 
#-----------------------------------

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

resMulti2 <- fabias(X,5,0.6,300)

extractPlot(resMulti2,ti="FABIAS Multiple tissues(Su)")

raMulti2 <- extractBic(resMulti2)

if ((raMulti2$bic[[1]][1]>1) && (raMulti2$bic[[1]][2])>1) {
    plotBicluster(raMulti2,1)
}
if ((raMulti2$bic[[2]][1]>1) && (raMulti2$bic[[2]][2])>1) {
    plotBicluster(raMulti2,2)
}
if ((raMulti2$bic[[3]][1]>1) && (raMulti2$bic[[3]][2])>1) {
    plotBicluster(raMulti2,3)
}
if ((raMulti2$bic[[4]][1]>1) && (raMulti2$bic[[4]][2])>1) {
    plotBicluster(raMulti2,4)
}

plot(resMulti2,dim=c(1,2),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti2,dim=c(1,3),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti2,dim=c(1,4),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti2,dim=c(1,5),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti2,dim=c(2,3),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti2,dim=c(2,4),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti2,dim=c(2,5),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti2,dim=c(3,4),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti2,dim=c(3,5),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti2,dim=c(4,5),label.tol=0.01,col.group=CMulti,lab.size=0.6)

}


#-----------------------------------------
# DEMO4: Rosenwald's diffuse large-B-cell
#        lymphoma gene expression data set 
#-----------------------------------------


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

resDLBCL2 <- fabias(X,5,0.6,300)

extractPlot(resDLBCL2,ti="FABIAS Lymphoma(Rosenwald)")

raDLBCL2 <- extractBic(resDLBCL2)

if ((raDLBCL2$bic[[1]][1]>1) && (raDLBCL2$bic[[1]][2])>1) {
    plotBicluster(raDLBCL2,1)
}
if ((raDLBCL2$bic[[2]][1]>1) && (raDLBCL2$bic[[2]][2])>1) {
    plotBicluster(raDLBCL2,2)
}
if ((raDLBCL2$bic[[3]][1]>1) && (raDLBCL2$bic[[3]][2])>1) {
    plotBicluster(raDLBCL2,3)
}
if ((raDLBCL2$bic[[4]][1]>1) && (raDLBCL2$bic[[4]][2])>1) {
    plotBicluster(raDLBCL2,4)
}

plot(resDLBCL2,dim=c(1,2),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL2,dim=c(1,3),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL2,dim=c(1,4),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL2,dim=c(1,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL2,dim=c(2,3),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL2,dim=c(2,4),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL2,dim=c(2,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL2,dim=c(3,4),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL2,dim=c(3,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL2,dim=c(4,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)

}



}
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
\concept{non-negative matrix factorization}
\concept{multivariate analysis}
\concept{latent variables}


