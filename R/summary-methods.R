#
#
# Author: SEPP HOCHREITER
###############################################################################


setMethod("summary", "Factorization",
function(object, ...)
{
    cat("\nAn object of class",class(object),"\n\n")
    cat("call:", deparse(object@parameters[[1]],0.75*getOption("width")),
        sep="\n\t")
    
    l<-object@l
    n<-object@n
    cat("\nNumber of rows: ",n, "\n")
    cat("\nNumber of columns: ",l, "\n")


    p1<-object@p1
    p2<-object@p2

    if (p1==p2) {
        cat("\nNumber of clusters: ",p2, "\n")
    } else {
        cat("\nNumber of row clusters: ",p1, "\n")
        cat("\nNumber of column clusters: ",p2, "\n")
    }


     cat("\n Information content of the clusters:\n")

    avini <- object@avini
    names(avini) <-paste("BC", c(1:p2,"sum"))


    print.default(round(avini,2), print.gap = 2)

     cat("\n\n")

     cat("\n Information content of the samples:\n")
    xavini <- object@xavini

    names(xavini) <-paste("Sample", c(1:l,"sum"))


    print.default(round(xavini,2), print.gap = 2)

     cat("\n Column clusters / Factors:\n")

    a <- summary.matrix(t(object@Z), digits = 3)
    colnames(a) <- paste("BC", c(1:p2))
    print(a)

     cat("\n")

     cat("\n Row clusters / Loadings:\n")
     a <- summary.matrix(object@L, digits = 3)
     colnames(a) <- paste("BC", c(1:p2))
   print(a)



})
