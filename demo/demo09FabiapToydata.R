#
#
# Author: SEPP HOCHREITER
###############################################################################

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

#### FABIAP

resToy3 <- fabiap(X,13,0.1,400,1.0,1.0,0.6,0.6)

rToy3 <- extractPlot(resToy3,ti="FABIAP",Y=Y)

raToy3 <- extractBic(resToy3)

if ((raToy3$bic[[1]][1]>1) && (raToy3$bic[[1]][2])>1) {
    plotBicluster(raToy3,1)
}
if ((raToy3$bic[[2]][1]>1) && (raToy3$bic[[2]][2])>1) {
    plotBicluster(raToy3,2)
}
if ((raToy3$bic[[3]][1]>1) && (raToy3$bic[[3]][2])>1) {
    plotBicluster(raToy3,3)
}
if ((raToy3$bic[[4]][1]>1) && (raToy3$bic[[4]][2])>1) {
    plotBicluster(raToy3,4)
}

colnames(resToy3@X) <- clab

rownames(resToy3@X) <- llab


plot(resToy3,dim=c(1,2),label.tol=0.1,col.group = groups,lab.size=0.6)
plot(resToy3,dim=c(1,3),label.tol=0.1,col.group = groups,lab.size=0.6)
plot(resToy3,dim=c(2,3),label.tol=0.1,col.group = groups,lab.size=0.6)

