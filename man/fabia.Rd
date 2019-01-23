\name{fabia}
\alias{fabia}
\title{Factor Analysis for Bicluster Acquisition: Laplace Prior (FABIA)}
\description{

  \code{fabia}: C implementation of  \code{fabia}.
}
\usage{

fabia(X,p=13,alpha=0.01,cyc=500,spl=0,spz=0.5,non_negative=0,random=1.0,center=2,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0)

}
\arguments{
  \item{X}{the data matrix.}
  \item{p}{number of hidden factors = number of biclusters; default = 13.}
  \item{alpha}{sparseness loadings (0 - 1.0); default = 0.01.}
  \item{cyc}{number of iterations; default = 500.}
  \item{spl}{sparseness prior loadings (0 - 2.0); default = 0 (Laplace).}
  \item{spz}{sparseness factors (0.5 - 2.0); default = 0.5 (Laplace).}
  \item{non_negative}{Non-negative factors and loadings if non_negative
    > 0; default = 0.}
  \item{random}{<=0: by SVD, >0: random initialization of loadings in
    [-random,random]; default = 1.0.}
  \item{center}{data centering: 1 (mean), 2
    (median), > 2 (mode), 0 (no); default = 2.}
  \item{norm}{data normalization:   1 (0.75-0.25 quantile), >1 (var=1),
    0 (no); default = 1.}
  \item{scale}{loading vectors are scaled in each iteration to the given
    variance. 0.0 indicates non scaling; default = 0.0.}
  \item{lap}{minimal value of the variational parameter; default = 1.0} 
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
  \code{LZ} (estimated noise free data \eqn{L  Z}),
  \code{L} (loadings  \eqn{L}),
  \code{Z} (factors  \eqn{Z}),
  \code{U} (noise: \eqn{X-LZ}),
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
\code{\link{readSpfabiaResult}},
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


resEx <- fabia(X,3,0.01,50)


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



#### FABIA

resToy1 <- fabia(X,13,0.01,400)

extractPlot(resToy1,ti="FABIA",Y=Y)

raToy1 <- extractBic(resToy1)

if ((raToy1$bic[[1]][1]>1) && (raToy1$bic[[1]][2])>1) {
    plotBicluster(raToy1,1)
}
if ((raToy1$bic[[2]][1]>1) && (raToy1$bic[[2]][2])>1) {
plotBicluster(raToy1,2)
}
if ((raToy1$bic[[3]][1]>1) && (raToy1$bic[[3]][2])>1) {
plotBicluster(raToy1,3)
}
if ((raToy1$bic[[4]][1]>1) && (raToy1$bic[[4]][2])>1) {
plotBicluster(raToy1,4)
}



colnames(X(resToy1)) <- clab

rownames(X(resToy1)) <- llab


plot(resToy1,dim=c(1,2),label.tol=0.1,col.group = groups,lab.size=0.6)
plot(resToy1,dim=c(1,3),label.tol=0.1,col.group = groups,lab.size=0.6)
plot(resToy1,dim=c(2,3),label.tol=0.1,col.group = groups,lab.size=0.6)



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


resBreast1 <- fabia(X,5,0.1,400)

extractPlot(resBreast1,ti="FABIA Breast cancer(Veer)")


raBreast1 <- extractBic(resBreast1)

if ((raBreast1$bic[[1]][1]>1) && (raBreast1$bic[[1]][2])>1) {
    plotBicluster(raBreast1,1)
}
if ((raBreast1$bic[[2]][1]>1) && (raBreast1$bic[[2]][2])>1) {
    plotBicluster(raBreast1,2)
}
if ((raBreast1$bic[[3]][1]>1) && (raBreast1$bic[[3]][2])>1) {
    plotBicluster(raBreast1,3)
}
if ((raBreast1$bic[[4]][1]>1) && (raBreast1$bic[[4]][2])>1) {
    plotBicluster(raBreast1,4)
}


plot(resBreast1,dim=c(1,2),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast1,dim=c(1,3),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast1,dim=c(1,4),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast1,dim=c(1,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast1,dim=c(2,3),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast1,dim=c(2,4),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast1,dim=c(2,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast1,dim=c(3,4),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast1,dim=c(3,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast1,dim=c(4,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)

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

resMulti1 <- fabia(X,5,0.06,300,norm=2)

extractPlot(resMulti1,ti="FABIA Multiple tissues(Su)")

raMulti1 <- extractBic(resMulti1)

if ((raMulti1$bic[[1]][1]>1) && (raMulti1$bic[[1]][2])>1) {
    plotBicluster(raMulti1,1)
}
if ((raMulti1$bic[[2]][1]>1) && (raMulti1$bic[[2]][2])>1) {
    plotBicluster(raMulti1,2)
}
if ((raMulti1$bic[[3]][1]>1) && (raMulti1$bic[[3]][2])>1) {
    plotBicluster(raMulti1,3)
}
if ((raMulti1$bic[[4]][1]>1) && (raMulti1$bic[[4]][2])>1) {
    plotBicluster(raMulti1,4)
}

plot(resMulti1,dim=c(1,2),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti1,dim=c(1,3),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti1,dim=c(1,4),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti1,dim=c(1,5),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti1,dim=c(2,3),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti1,dim=c(2,4),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti1,dim=c(2,5),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti1,dim=c(3,4),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti1,dim=c(3,5),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti1,dim=c(4,5),label.tol=0.01,col.group=CMulti,lab.size=0.6)

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

resDLBCL1 <- fabia(X,5,0.1,400,norm=2)

extractPlot(resDLBCL1,ti="FABIA Lymphoma(Rosenwald)")

raDLBCL1 <- extractBic(resDLBCL1)

if ((raDLBCL1$bic[[1]][1]>1) && (raDLBCL1$bic[[1]][2])>1) {
    plotBicluster(raDLBCL1,1)
}
if ((raDLBCL1$bic[[2]][1]>1) && (raDLBCL1$bic[[2]][2])>1) {
    plotBicluster(raDLBCL1,2)
}
if ((raDLBCL1$bic[[3]][1]>1) && (raDLBCL1$bic[[3]][2])>1) {
    plotBicluster(raDLBCL1,3)
}
if ((raDLBCL1$bic[[4]][1]>1) && (raDLBCL1$bic[[4]][2])>1) {
    plotBicluster(raDLBCL1,4)
}

plot(resDLBCL1,dim=c(1,2),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL1,dim=c(1,3),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL1,dim=c(1,4),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL1,dim=c(1,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL1,dim=c(2,3),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL1,dim=c(2,4),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL1,dim=c(2,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL1,dim=c(3,4),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL1,dim=c(3,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL1,dim=c(4,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)

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


