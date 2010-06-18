#Data from Broad Institute "Cancer Program Data Sets"
#
#http://www.broadinstitute.org/cgi-bin/cancer/datasets.cgi
#
#Subclass Mapping: Identifying Common Subtypes in Independent Disease Data Sets
#PLoS ONE 2(11): e1195. doi:10.1371/journal.pone.0001195
#
#"DLBCL_B"
#
#Rosenwald A, Wright G, Chan WC, Connors JM, Campo E, et al. (2002)
#  The use of molecular profiling to predict survival after chemotherapy
#  for diffuse large-B-cell lymphoma.
#N Engl J Med 346: 1937-1947.
#
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


data(DLBCL_B)

X <- as.matrix(XDLBCL)


resDLBCL3 <- fabiap(X,5,0.1,400)

rDLBCL3 <- extractPlot(resDLBCL3,ti="FABIAP Lymphoma(Rosenwald)")
raDLBCL3 <- extractBic(resDLBCL3)

if ((raDLBCL3$bic[[1]][1]>1) && (raDLBCL3$bic[[1]][2])>1) {
    plotBicluster(raDLBCL3,1)
}
if ((raDLBCL3$bic[[2]][1]>1) && (raDLBCL3$bic[[2]][2])>1) {
    plotBicluster(raDLBCL3,2)
}
if ((raDLBCL3$bic[[3]][1]>1) && (raDLBCL3$bic[[3]][2])>1) {
    plotBicluster(raDLBCL3,3)
}
if ((raDLBCL3$bic[[4]][1]>1) && (raDLBCL3$bic[[4]][2])>1) {
    plotBicluster(raDLBCL3,4)
}

plot(resDLBCL3,dim=c(1,2),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL3,dim=c(1,3),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL3,dim=c(1,4),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL3,dim=c(1,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL3,dim=c(2,3),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL3,dim=c(2,4),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL3,dim=c(2,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL3,dim=c(3,4),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL3,dim=c(3,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
plot(resDLBCL3,dim=c(4,5),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)


}

