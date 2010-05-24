\name{mfsc}
\alias{mfsc}
\title{Sparse Matrix Factorization for Bicluster Analysis (MFSC)}
\description{

  \code{mfsc}: \R implementation of  \code{mfsc}.
}
\usage{

mfsc(X,p=5,cyc=100,sL=0.6,sZ=0.6,center=2,norm=1)


}
\arguments{
  \item{X}{the data matrix.}
  \item{p}{number of hidden factors = number of biclusters; default = 5.}
  \item{cyc}{maximal number of iterations; default = 100.}
  \item{sL}{final sparseness loadings; default = 0.6.}
  \item{sZ}{final sparseness factors; default = 0.6.}
  \item{center}{data centering: 1 (mean), 2
    (median), > 2 (mode), 0 (no); default = 2.}
  \item{norm}{data normalization:   1 (0.75-0.25 quantile), >1 (var=1), 0 (no); default = 1.}
}
\details{

  Biclusters are found by sparse matrix factorization
  where \emph{both} factors are sparse.

  Essentially the model is the sum of
  outer products of vectors:
   \deqn{X  =  \sum_{i=1}^{p} \lambda_i  z_i^T}
  where the number of summands  \eqn{p}
  is the number of biclusters.
  The matrix factorization is
  \deqn{X  =  L  Z }
  

   Here \eqn{\lambda_i} are from \eqn{R^n}, \eqn{z_i} from
   \eqn{R^l}, \eqn{L} from \eqn{R^{n \times p}},
   \eqn{Z} from \eqn{R^{p \times l}}, and \eqn{X}
   from \eqn{R^{n \times l}}.
   
   \bold{No noise assumption:}
  In contrast to factor analysis there is no noise assumption.

  If the nonzero components of
  the sparse vectors are grouped together then the outer product results in
  a matrix with a nonzero block and zeros elsewhere.



	 The model selection is performed by a constraint optimization
	 according
  to  Hoyer, 2004. The Euclidean distance (the Frobenius norm) is
  minimized subject to sparseness constraints.

  Model selection is done by gradient descent on the Euclidean
  objective and thereafter projection of single vectors of
  \eqn{L} and single vectors of \eqn{Z}
  to fulfill the sparseness constraints.

  The projection minimize the Euclidean distance
  to the original vector  given an
  \eqn{l_1}-norm and an \eqn{l_2}-norm.

  The projection is a convex quadratic problem which is solved
  iteratively where at each iteration at least one component is set to
  zero. Instead of the  \eqn{l_1}-norm a sparseness measurement is used
  which relates the \eqn{l_1}-norm to the \eqn{l_2}-norm.

  The code is implemented in \R.

}
\value{
 \item{}{object of the class \code{Factorization}. Containing
  \code{LZ} (estimated noise free data \eqn{L  Z}),
  \code{L} (loadings \eqn{L}),
  \code{Z} (factors \eqn{Z}),
  \code{U} (noise \eqn{X-LZ}),
  \code{center} (centering vector),
  \code{scaleData} (scaling vector),
  \code{X} (centered and scaled data \eqn{X})
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


resEx <- mfsc(X,3,30,0.6,0.6)

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

#### MFSC

resToy4 <- mfsc(X,13,100,0.6,0.6)

rToy4 <- extractPlot(resToy4,ti="MFSC",Y=Y)

raToy4 <- extractBic(resToy4,thresZ=0.01,thresL=0.05)

if ((raToy4$bic[[1]][1]>1) && (raToy4$bic[[1]][2])>1) {
    plotBicluster(raToy4,1)
}
if ((raToy4$bic[[2]][1]>1) && (raToy4$bic[[2]][2])>1) {
    plotBicluster(raToy4,2)
}
if ((raToy4$bic[[3]][1]>1) && (raToy4$bic[[3]][2])>1) {
    plotBicluster(raToy4,3)
}
if ((raToy4$bic[[4]][1]>1) && (raToy4$bic[[4]][2])>1) {
    plotBicluster(raToy4,4)
}

colnames(resToy4@X) <- clab

rownames(resToy4@X) <- llab


plot(resToy4,dim=c(1,2),label.tol=0.1,col.group = groups,lab.size=0.6)
plot(resToy4,dim=c(1,3),label.tol=0.1,col.group = groups,lab.size=0.6)
plot(resToy4,dim=c(2,3),label.tol=0.1,col.group = groups,lab.size=0.6)


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

resBreast4 <- mfsc(X,5,100,0.6,0.6)

rBreast4 <- extractPlot(resBreast4,ti="MFSC Breast cancer(Veer)")


raBreast4 <- extractBic(resBreast4,thresZ=0.01,thresL=0.05)

if ((raBreast4$bic[[1]][1]>1) && (raBreast4$bic[[1]][2])>1) {
    plotBicluster(raBreast4,1)
}
if ((raBreast4$bic[[2]][1]>1) && (raBreast4$bic[[2]][2])>1) {
    plotBicluster(raBreast4,2)
}
if ((raBreast4$bic[[3]][1]>1) && (raBreast4$bic[[3]][2])>1) {
    plotBicluster(raBreast4,3)
}
if ((raBreast4$bic[[4]][1]>1) && (raBreast4$bic[[4]][2])>1) {
    plotBicluster(raBreast4,4)
}

plot(resBreast4,dim=c(1,2),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast4,dim=c(1,3),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast4,dim=c(1,4),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast4,dim=c(1,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast4,dim=c(2,3),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast4,dim=c(2,4),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast4,dim=c(2,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast4,dim=c(3,4),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast4,dim=c(3,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast4,dim=c(4,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)

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

resMulti4 <- mfsc(X,5,100,0.6,0.6)

rMulti4 <- extractPlot(resMulti4,ti="MFSC Multiple tissues(Su)")


raMulti4 <- extractBic(resMulti4,thresZ=0.01,thresL=0.05)

if ((raMulti4$bic[[1]][1]>1) && (raMulti4$bic[[1]][2])>1) {
    plotBicluster(raMulti4,1)
}
if ((raMulti4$bic[[2]][1]>1) && (raMulti4$bic[[2]][2])>1) {
    plotBicluster(raMulti4,2)
}
if ((raMulti4$bic[[3]][1]>1) && (raMulti4$bic[[3]][2])>1) {
    plotBicluster(raMulti4,3)
}
if ((raMulti4$bic[[4]][1]>1) && (raMulti4$bic[[4]][2])>1) {
    plotBicluster(raMulti4,4)
}

plot(resMulti4,dim=c(1,2),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti4,dim=c(1,3),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti4,dim=c(1,4),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti4,dim=c(1,5),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti4,dim=c(2,3),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti4,dim=c(2,4),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti4,dim=c(2,5),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti4,dim=c(3,4),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti4,dim=c(3,5),label.tol=0.01,col.group=CMulti,lab.size=0.6)
plot(resMulti4,dim=c(4,5),label.tol=0.01,col.group=CMulti,lab.size=0.6)

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


resDLBCL4 <- mfsc(X,5,100,0.6,0.6)

rDLBCL4 <- extractPlot(resDLBCL4,ti="MFSC Lymphoma(Rosenwald)")

raDLBCL4 <- extractBic(resDLBCL4,thresZ=0.01,thresL=0.05)

if ((raDLBCL4$bic[[1]][1]>1) && (raDLBCL4$bic[[1]][2])>1) {
    plotBicluster(raDLBCL4,1)
}
if ((raDLBCL4$bic[[2]][1]>1) && (raDLBCL4$bic[[2]][2])>1) {
    plotBicluster(raDLBCL4,2)
}
if ((raDLBCL4$bic[[3]][1]>1) && (raDLBCL4$bic[[3]][2])>1) {
    plotBicluster(raDLBCL4,3)
}
if ((raDLBCL4$bic[[4]][1]>1) && (raDLBCL4$bic[[4]][2])>1) {
    plotBicluster(raDLBCL4,4)
}

plot(resDLBCL4,dim=c(1,2),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL4,dim=c(1,3),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL4,dim=c(1,4),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL4,dim=c(1,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL4,dim=c(2,3),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL4,dim=c(2,4),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL4,dim=c(2,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL4,dim=c(3,4),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL4,dim=c(3,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL4,dim=c(4,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)

}



}
}
\references{
  S. Hochreiter et al.,
  \sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
  Bioinformatics, doi:10.1093/bioinformatics/btq227, 2010.
  http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btq227
  
  Patrik O. Hoyer,
  \sQuote{Non-negative Matrix Factorization with Sparseness Constraints},
  Journal of Machine Learning Research 5:1457-1469, 2004.
}
\keyword{methods}
\keyword{cluster}
\concept{biclustering}
\concept{sparse coding}
\concept{non-negative matrix factorization}

