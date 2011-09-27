#Subclass Mapping: Identifying Common Subtypes in Independent Disease Data Sets
#PLoS ONE 2(11): e1195. doi:10.1371/journal.pone.0001195
#
#"Multi_A" without Array BR_U16
#
#Su AI, Cooke MP, Ching KA, Hakak Y, Walker JR, et al. (2002)
#Large-scale analysis of the human and mouse transcriptomes.
#Proc Natl Acad Sci U S A 99:4465-4470.
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

data(Multi_A)

X <- as.matrix(XMulti)

resMulti4 <- mfsc(X,5,100,0.6,0.6)

extractPlot(resMulti4,ti="MFSC Multiple tissues(Su)")


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

