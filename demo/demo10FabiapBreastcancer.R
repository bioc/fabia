#Data from Broad Institute "Cancer Program Data Sets"
#
#http://www.broadinstitute.org/cgi-bin/cancer/datasets.cgi
#
#Hoshida Y, Brunet J-P, Tamayo P, Golub TR, Mesirov JP (2007)
#Subclass Mapping: Identifying Common Subtypes in Independent Disease Data Sets
#PLoS ONE 2(11): e1195. doi:10.1371/journal.pone.0001195
#
#"Breast_A" without Array S54
#van't Veer LJ, Dai H, van de Vijver MJ, He YD, Hart AA, et al. (2002)
#Gene expression profiling predicts clinical outcome of breast cancer.
#Nature 415:530-536.
#
# Author: SEPP HOCHREITER
###############################################################################

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


resBreast3 <- fabiap(X,5,0.1,400,1.0,1.0,0.6,0.6)

rBreast3 <- extractPlot(resBreast3,ti="FABIAP Breast cancer(Veer)")

raBreast3 <- extractBic(resBreast3)

if ((raBreast3$bic[[1]][1]>1) && (raBreast3$bic[[1]][2])>1) {
    plotBicluster(raBreast3,1)
}
if ((raBreast3$bic[[2]][1]>1) && (raBreast3$bic[[2]][2])>1) {
    plotBicluster(raBreast3,2)
}
if ((raBreast3$bic[[3]][1]>1) && (raBreast3$bic[[3]][2])>1) {
    plotBicluster(raBreast3,3)
}
if ((raBreast3$bic[[4]][1]>1) && (raBreast3$bic[[4]][2])>1) {
    plotBicluster(raBreast3,4)
}

plot(resBreast3,dim=c(1,2),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast3,dim=c(1,3),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast3,dim=c(1,4),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast3,dim=c(1,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast3,dim=c(2,3),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast3,dim=c(2,4),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast3,dim=c(2,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast3,dim=c(3,4),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast3,dim=c(3,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)
plot(resBreast3,dim=c(4,5),label.tol=0.03,col.group=CBreast,lab.size=0.6)

}

