#
#
# Author: SEPP HOCHREITER
###############################################################################


setMethod("showSelected", "Factorization", function(object, which=c(1,2,3,4))
{
    l<-object@l
    n<-object@n
    p1<-object@p1
    p2<-object@p2

    devAskNewPage(ask = FALSE)
    if (length(which) > 1 && dev.interactive()) {
        devAskNewPage(ask = TRUE)
    }

    showf <- c(FALSE, FALSE,FALSE, FALSE)
    showf[which] <- TRUE


    if (showf[1]){
        names <- paste("BC", 1:p2)
        title <- "Information Content of Biclusters"

        barplot(object@avini[1:p2],names.arg=names,main=title)
    }


    if (showf[2]){

        names1 <- paste("Sample", 1:l)
        title1 <- "Information Content of Samples"

        barplot(object@xavini[1:l],names.arg=names1,main=title1)
    }


    if (showf[3]){

        names <- paste("BC", 1:p2)

        title2 <- "Loadings of the Biclusters"

        boxplot(object@L,names=names,main=title2)
    }

    if (showf[4]){

        names <- paste("BC", 1:p2)

        title3 <- "Factors of the Biclusters"

        boxplot(t(object@Z),names=names,main=title3)
    }

    devAskNewPage(ask = FALSE)


})
