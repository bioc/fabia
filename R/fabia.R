#
#
# Author: SEPP HOCHREITER
###############################################################################


fabia <- function(X,p=5,alpha=0.1,cyc=500,spl=0,spz=0.5,non_negative=0,random=1.0,center=2,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0){
	## X - data matrix
	## cyc - maximum number of cycles
        ## alpha - sparseness
        ## p - factors

        if (missing(X)) {
            stop("Data matrix X missing. Stopped.")
        }

        if (!is.matrix(X)) {
            X <- as.matrix(X)
        }

        if (!is.numeric(X)) {
            stop("Data matrix X must be numeric. Stopped.")
        }

        l=ncol(X)
        n=nrow(X)

        if (p>min(l,n)) {
            stop("Too many biclusters. Stopped.")
        }

        rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "gene")
        colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
        rowna <- rownames(X)
        colna <- colnames(X)

       message("Running FABIA on a ",n,"x",l," matrix with")
       message("   Number of biclusters ---------------- p: ",p)
       message("   Sparseness factor --------------- alpha: ",alpha)
       message("   Number of iterations -------------- cyc: ",cyc)
       message("   Loading prior parameter ----------- spl: ",spl)
       message("   Factor prior parameter ------------ spz: ",spz)
       if (random>0) {
       message("   Initialization loadings--------- random: ",random," = interval")
       } else {
       message("   Initialization loadings--------- random: ",random," = SVD")
      }
       if (non_negative>0) {
       message("   Nonnegative Loadings and Factors ------: ",non_negative," = Yes")
       } else {
       message("   Nonnegative Loadings and Factors ------: ",non_negative," = No")
      }
       if (center>0) {
            if (center<2) {
       message("   Centering ---------------------- center: 1 = mean")
            } else {
                if (center<3) {
       message("   Centering ---------------------- center: 2 = median")
                } else {
       message("   Centering ---------------------- center: ", center," = mode")
                }
            }
        } else {
       message("   Centering ---------------------- center: ", center," = no centering")
        }
       if (norm>0) {
           if (norm > 1) {
       message("   Scaling to variance one: --------- norm: ",norm," = Yes")
           } else {
       message("   Quantile scaling (0.75-0.25): ---- norm: ",norm," = Yes")
           }} else {
       message("   Scaling -------------------------- norm: ",norm," = No")
       }
       if (scale>0) {
       message("   Scaling loadings per iteration -- scale: ",scale," = to ", scale)
       } else {
       message("   Scaling loadings per iteration -- scale: ", scale," = No")
       }

       message("   Constraint variational parameter -- lap: ", lap)
       if ((nL>0)&&(nL<p)) {
           message("   Max. number of biclusters per row -- nL: ", nL)
           if (bL>0) {
           message("         starting at ------------------ bL: ", bL)
           } else {
           message("         starting at ------------------ bL: ", bL, " = from start")
           }
       } else {
           message("   Max. number of biclusters per row -- nL: ", nL, " = no limit")
       }
       if ((lL>0)&&(lL<n)) {
           message("   Max. number of row elements / biclu. lL: ", lL)
           if (bL>0) {
           message("         starting at ------------------ bL: ", bL)
           } else {
           message("         starting at ------------------ bL: ", bL, " = from start")
           }
       } else {
           message("   Max. number of row elements / biclu. lL: ", lL, " = no limit")
       }






        eps <- as.double(1e-3)
        eps1 <- as.double(1e-10)
        init_lapla <- 1.0
        init_psi <- 0.2


        iin <-  1.0/l

        cent <- as.vector(rep(0.0,n))
        if (center>0) {
            if (center<2) {
                cent <- apply(X, 1, mean)
            } else {
                if (center<3) {
                    cent <- apply(X, 1, median)
                } else {
                    cent <- estimateMode(X)
                }
            }
        }

        X <- X - cent

        XX <- as.vector(rep(1,n))

        if (norm>0) {
            if (norm>1)
            {
                scaleData <-  1/sqrt(iin*apply(X,1,function(x) sum(x^2))+0.001*XX)
            }
            else
            {
                scaleData <-  1/sqrt(apply(X,1,function(x) {quantile(x,0.75) - quantile(x,0.25)})+0.001*XX)
            }
            X <- scaleData*X
        } else {
            scaleData <- as.vector(rep(1.0,n))
        }



         if (random>0) {
             L <- matrix(random*rnorm(n*p),nrow=n,ncol=p)
         } else {
             svd <- La.svd(X)
             L <- svd$u[,1:p]%*%diag(svd$d[1:p])
         }


        if (scale>0)
        {

            dL <- 1/sqrt( apply(L,1,function(x) sum(x^2))+0.001*as.vector(rep(1,n)))
            L <- (scale*dL)*L
        }


        lapla <-  init_lapla*matrix(1,nrow=l,ncol=p)

        Psi <- init_psi*XX



	cyc <- as.integer(cyc)
	nL <- as.integer(nL)
	lL <- as.integer(lL)
	bL <- as.integer(bL)
	alpha <- as.double(alpha)
	non_negative <- as.integer(non_negative)
	p <- as.integer(p)
	spz <- as.double(spz)
	scale <- as.double(scale)
	lap <- as.double(lap)



	res <- .Call("fabic", X,Psi,L,lapla,cyc ,alpha,eps,eps1,spl,spz,scale,lap,nL,lL,bL,non_negative,PACKAGE="fabia")


        if (is.null(res))
        {
            return(new('Factorization', parameters=list(),n=1,p1=1,p2=1,l=1,center=as.vector(1),scaleData=as.vector(1),X=as.matrix(1),L=noL,Z=nZ,M=as.matrix(1),LZ=as.matrix(1),U=as.matrix(1),avini=as.vector(1),xavini=as.vector(1),ini=as.matrix(1),Psi=as.vector(1),lapla=as.matrix(1)))
        }


        # INI call for biclusters


        vz <- iin*apply(res$E_SX_n,1,function(x) sum(x^2))

        vz <- sqrt(vz+1e-10)

        ivz <- 1/vz

        if(length(ivz)==1) {
            nZ <- ivz*res$E_SX_n

            noL <- vz*res$L
        }
        else {
            nZ <- ivz*res$E_SX_n

            noL <- t(vz*t(res$L))
        }


        ini <- matrix(0,l,(p+1))
        avini <- as.vector(rep(0.0,(p+1)))
        xavini <- as.vector(rep(0.0,(l+1)))

        idp <- diag(p)
        ppL <- crossprod(noL,(1/res$Psi)*noL)


        for (j in 1:l){
            mat <- idp + ppL/res$lapla[j,]
            ini[j,1:p] <- log(diag(mat))
            s <- log(det(mat))
            ini[j,p+1] <- s
            xavini[j] <-s
        }
        for (i in 1:p){
            avini[i] <- sum(ini[,i])
        }

        ss <- sum(ini[,p+1])
        xavini[l+1] <- ss
        avini[p+1] <- ss


        Lz=noL%*%nZ


        rownames(noL) <- rowna
        colnames(noL) <- colnames(noL, do.NULL = FALSE, prefix = "bicluster")
        clnames <- colnames(noL)

        rownames(nZ) <- clnames
        colnames(nZ) <- colna

        M <- diag(p)
        rownames(M) <- clnames
        colnames(M) <- clnames

        rownames(Lz) <- rowna
        colnames(Lz) <- colna

        U <- X-Lz

        rownames(U) <- rowna
        colnames(U) <- colna


        if ((avini[p+1]>1e-8)&&(p>1)) {

            soo <- sort(avini[1:p], decreasing = TRUE,index.return=TRUE)

            avini[1:p] <- avini[soo$ix]
            noL <- noL[,soo$ix]
            nZ <- nZ[soo$ix,]
            M <- M[soo$ix,soo$ix]

        }




        return(new('Factorization', parameters=list("fabia",cyc,alpha,spl,spz,p,NULL,NULL,random,scale,norm,center,lap,nL,lL,bL,non_negative),n=n,p1=p,p2=p,l=l,center=cent,scaleData=scaleData,X=X,L=noL,Z=nZ,M=M,LZ=Lz,U=U,avini=avini,xavini=xavini,ini=ini,Psi=res$Psi,lapla=res$lapla))


}

fabiap <- function(X,p=5,alpha=0.1,cyc=500,spl=0,spz=0.5,sL=0.6,sZ=0.6,non_negative=0,random=1.0,center=2,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0){
	## X - data matrix
	## cyc - maximum number of cycles
        ## alpha - sparseness
        ## p - factors
        ## sL - final sparseness L
        ## sZ - final sparseness z

        if (missing(X)) {
            stop("Data matrix X missing. Stopped.")
        }

        if (!is.matrix(X)) {
            X <- as.matrix(X)
        }

        if (!is.numeric(X)) {
            stop("Data matrix X must be numeric. Stopped.")
        }

        l=ncol(X)
        n=nrow(X)

        if (p>min(l,n)) {
            stop("Too many biclusters. Stopped.")
        }

        rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "gene")
        colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
        rowna <- rownames(X)
        colna <- colnames(X)



       message("Running FABIAP on a ",n,"x",l," matrix with")
       message("   Number of biclusters ---------------- p: ",p)
       message("   Sparseness factor --------------- alpha: ",alpha)
       message("   Number of iterations -------------- cyc: ",cyc)
       message("   Loading prior parameter ----------- spl: ",spl)
       message("   Factor prior parameter ------------ spz: ",spz)
       message("   Loading sparseness projection ------ sL: ",sL)
       message("   Factor sparseness projection ------- sZ: ",sZ)
       if (random>0) {
       message("   Initialization loadings--------- random: ",random," = interval")
       } else {
       message("   Initialization loadings--------- random: ",random," = SVD")
       }
       if (non_negative>0) {
       message("   Nonnegative Loadings and Factors ------: ",non_negative," = Yes")
       } else {
       message("   Nonnegative Loadings and Factors ------: ",non_negative," = No")
      }
       if (center>0) {
            if (center<2) {
       message("   Centering ---------------------- center: 1 = mean")
            } else {
                if (center<3) {
       message("   Centering ---------------------- center: 2 = median")
                } else {
       message("   Centering ---------------------- center: ", center," = mode")
                }
            }
        } else {
       message("   Centering ---------------------- center: ", center," = no centering")
        }
       if (norm>0) {
           if (norm > 1) {
       message("   Scaling to variance one: --------- norm: ",norm," = Yes")
           } else {
       message("   Quantile scaling (0.75-0.25): ---- norm: ",norm," = Yes")
           }} else {
       message("   Scaling -------------------------- norm: ",norm," = No")
       }
       if (scale>0) {
       message("   Scaling loadings per iteration -- scale: ",scale," = to ", scale)
       } else {
       message("   Scaling loadings per iteration -- scale: ", scale," = No")
       }

       message("   Constraint variational parameter -- lap: ", lap)
       if ((nL>0)&&(nL<p)) {
           message("   Max. number of biclusters per row -- nL: ", nL)
           if (bL>0) {
           message("         starting at ------------------ bL: ", bL)
           } else {
           message("         starting at ------------------ bL: ", bL, " = from start")
           }
       } else {
           message("   Max. number of biclusters per row -- nL: ", nL, " = no limit")
       }
       if ((lL>0)&&(lL<n)) {
           message("   Max. number of row elements / biclu. lL: ", lL)
           if (bL>0) {
           message("         starting at ------------------ bL: ", bL)
           } else {
           message("         starting at ------------------ bL: ", bL, " = from start")
           }
       } else {
           message("   Max. number of row elements / biclu. lL: ", lL, " = no limit")
       }






        eps <- as.double(1e-3)
        eps1 <- as.double(1e-10)
        init_lapla <- 1.0
        init_psi <- 0.2

        iin <-  1.0/l

        cent <- as.vector(rep(0.0,n))
        if (center>0) {
            if (center<2) {
                cent <- apply(X, 1, mean)
            } else {
                if (center<3) {
                    cent <- apply(X, 1, median)
                } else {
                    cent <- estimateMode(X)
                }
            }
        }

        X <- X - cent

        XX <- as.vector(rep(1,n))

        if (norm>0) {
            if (norm>1)
            {
                scaleData <-  1/sqrt(iin*apply(X,1,function(x) sum(x^2))+0.001*XX)
            }
            else
            {
                scaleData <-  1/sqrt(apply(X,1,function(x) {quantile(x,0.75) - quantile(x,0.25)})+0.001*XX)
            }
            X <- scaleData*X
        } else {
            scaleData <- as.vector(rep(1.0,n))
        }



        if (random>0) {
             L <- matrix(random*rnorm(n*p),nrow=n,ncol=p)
             if (non_negative>0)
             {
                 L <- abs(L)
             }
         } else {
             svd <- La.svd(X)
             L <- svd$u[,1:p]%*%diag(svd$d[1:p])
             if (non_negative>0)
             {
                 L <- abs(L)
             }
         }


        if (scale>0)
        {

            dL <- 1/sqrt( apply(L,1,function(x) sum(x^2))+0.001*as.vector(rep(1,n)))
            L <- (scale*dL)*L
        }


        lapla <-  init_lapla*matrix(1,nrow=l,ncol=p)

        Psi <- init_psi*XX


	cyc <- as.integer(cyc)
	nL <- as.integer(nL)
	lL <- as.integer(lL)
	bL <- as.integer(bL)
	alpha <- as.double(alpha)
	non_negative <- as.integer(non_negative)
	p <- as.integer(p)
	spz <- as.double(spz)
	scale <- as.double(scale)
	lap <- as.double(lap)


	res <- .Call("fabic", X,Psi,L,lapla,cyc,alpha,eps,eps1,spl,spz,scale,lap,nL,lL,bL,non_negative,PACKAGE="fabia")

        if (is.null(res))
        {
            return(new('Factorization', parameters=list(),n=1,p1=1,p2=1,l=1,center=as.vector(1),scaleData=as.vector(1),X=as.matrix(1),L=noL,Z=nZ,M=as.matrix(1),LZ=as.matrix(1),U=as.matrix(1),avini=as.vector(1),xavini=as.vector(1),ini=as.matrix(1),Psi=as.vector(1),lapla=as.matrix(1)))
        }

        rL <- res$L
        n <- nrow(rL)
        gh1 <- sqrt(1.0*n)-(sqrt(1.0*n)-1.0)*sL
        for (i in 1:p)
        {
            #le <- sum(res$L[,i]^2)
            le <- 1.0
            rL[,i] <- projFunc(s=as.vector(res$L[,i]),k1=gh1,k2=le)
        }

        rz <- res$E_SX_n
        l <- ncol(rz)
        gh2 <- sqrt(1.0*l)-(sqrt(1.0*l)-1.0)*sZ
        for (i in 1:p)
        {
            #le <- sum(res$z[,i]^2)
            le <- 1.0
            rz[i,] <- projFunc(s=as.vector(res$E_SX_n[i,]),k1=2,k2=le)
        }


        # INI call for biclusters

        vz <- iin*apply(rz,1,function(x) sum(x^2))

        vz <- sqrt(vz+1e-10)

        ivz <- 1/vz

        if(length(ivz)==1) {
            nZ <- ivz*res$E_SX_n

            noL <- vz*res$L
        }
        else {
            nZ <- ivz*res$E_SX_n

            noL <- t(vz*t(res$L))
        }




        ini <- matrix(0,l,(p+1))
        avini <- as.vector(rep(0.0,(p+1)))
        xavini <- as.vector(rep(0.0,(l+1)))

        idp <- diag(p)
        ppL <- crossprod(noL,(1/res$Psi)*noL)


        for (j in 1:l){
            mat <- idp + ppL/res$lapla[j,]
            ini[j,1:p] <- log(diag(mat))
            s <- log(det(mat))
            ini[j,p+1] <- s
            xavini[j] <-s
        }
        for (i in 1:p){
            avini[i] <- sum(ini[,i])
        }

        ss <- sum(ini[,p+1])
        xavini[l+1] <- ss
        avini[p+1] <- ss


        Lz=noL%*%nZ


        rownames(noL) <- rowna
        colnames(noL) <- colnames(noL, do.NULL = FALSE, prefix = "bicluster")
        clnames <- colnames(noL)

        rownames(nZ) <- clnames
        colnames(nZ) <- colna

        M <- diag(p)
        rownames(M) <- clnames
        colnames(M) <- clnames

        rownames(Lz) <- rowna
        colnames(Lz) <- colna

        U <- X-Lz

        rownames(U) <- rowna
        colnames(U) <- colna



        if ((avini[p+1]>1e-8)&&(p>1)) {

            soo <- sort(avini[1:p], decreasing = TRUE,index.return=TRUE)

            avini[1:p] <- avini[soo$ix]
            noL <- noL[,soo$ix]
            nZ <- nZ[soo$ix,]
            M <- M[soo$ix,soo$ix]

        }



    return(new('Factorization', parameters=list("fabiap",cyc,alpha,spl,spz,p,sL,sZ,random,scale,norm,center,lap,nL,lL,bL,non_negative),n=n,p1=p,p2=p,l=l,center=cent,scaleData=scaleData,X=X,L=noL,Z=nZ,M=M,LZ=Lz,U=U,avini=avini,xavini=xavini,ini=ini,Psi=res$Psi,lapla=res$lapla))

}

fabias <- function(X,p=5,alpha=0.6,cyc=500,spz=0.5,non_negative=0,random=1.0,center=2,norm=1,lap=1.0,nL=0,lL=0,bL=0){
	## X - data matrix
	## cyc - maximum number of cycles
        ## alpha - sparseness low value because enforced by projFunc
        ## p - factors

        if (missing(X)) {
            stop("Data matrix X missing. Stopped.")
        }

        if (!is.matrix(X)) {
            X <- as.matrix(X)
        }

        if (!is.numeric(X)) {
            stop("Data matrix X must be numeric. Stopped.")
        }


        l=ncol(X)
        n=nrow(X)

        if (p>min(l,n)) {
            stop("Too many biclusters. Stopped.")
        }

   rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "gene")
        colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
        rowna <- rownames(X)
        colna <- colnames(X)



       message("Running FABIAS on a ",n,"x",l," matrix with")
       message("   Number of biclusters ---------------- p: ",p)
       message("   Sparseness projection factor ---- alpha: ",alpha)
       message("   Number of iterations -------------- cyc: ",cyc)
       message("   Factor prior parameter ------------ spz: ",spz)
       if (random>0) {
       message("   Initialization loadings--------- random: ",random," = interval")
       } else {
       message("   Initialization loadings--------- random: ",random," = SVD")
       }
       if (non_negative>0) {
       message("   Nonnegative Loadings and Factors ------: ",non_negative," = Yes")
       } else {
       message("   Nonnegative Loadings and Factors ------: ",non_negative," = No")
      }
       if (center>0) {
            if (center<2) {
       message("   Centering ---------------------- center: 1 = mean")
            } else {
                if (center<3) {
       message("   Centering ---------------------- center: 2 = median")
                } else {
       message("   Centering ---------------------- center: ", center," = mode")
                }
            }
        } else {
       message("   Centering ---------------------- center: ", center," = no centering")
        }
       if (norm>0) {
           if (norm > 1) {
       message("   Scaling to variance one: --------- norm: ",norm," = Yes")
           } else {
       message("   Quantile scaling (0.75-0.25): ---- norm: ",norm," = Yes")
           }} else {
       message("   Scaling -------------------------- norm: ",norm," = No")
       }

       message("   Constraint variational parameter -- lap: ", lap)

       if ((nL>0)&&(nL<p)) {
           message("   Max. number of biclusters per row -- nL: ", nL)
           if (bL>0) {
           message("         starting at ------------------ bL: ", bL)
           } else {
           message("         starting at ------------------ bL: ", bL, " = from start")
           }
       } else {
           message("   Max. number of biclusters per row -- nL: ", nL, " = no limit")
       }
       if ((lL>0)&&(lL<n)) {
           message("   Max. number of row elements / biclu. lL: ", lL)
           if (bL>0) {
           message("         starting at ------------------ bL: ", bL)
           } else {
           message("         starting at ------------------ bL: ", bL, " = from start")
           }
       } else {
           message("   Max. number of row elements / biclu. lL: ", lL, " = no limit")
       }




       eps <- as.double(1e-3)
        eps1 <- as.double(1e-10)
        init_lapla <- 1.0
        init_psi <- 0.2

         iin <-  1.0/l

        cent <- as.vector(rep(0.0,n))
        if (center>0) {
            if (center<2) {
                cent <- apply(X, 1, mean)
            } else {
                if (center<3) {
                    cent <- apply(X, 1, median)
                } else {
                    cent <- estimateMode(X)
                }
            }
        }

        X <- X - cent

        XX <- as.vector(rep(1,n))

        if (norm>0) {
            if (norm>1)
            {
                scaleData <-  1/sqrt(iin*apply(X,1,function(x) sum(x^2))+0.001*XX)
            }
            else
            {
                scaleData <-  1/sqrt(apply(X,1,function(x) {quantile(x,0.75) - quantile(x,0.25)})+0.001*XX)
            }
            X <- scaleData*X
        } else {
            scaleData <- as.vector(rep(1.0,n))
        }



        if (random>0) {
             L <- matrix(random*rnorm(n*p),nrow=n,ncol=p)
             if (non_negative>0)
             {
                 L <- abs(L)
             }
         } else {
             svd <- La.svd(X)
             L <- svd$u[,1:p]%*%diag(svd$d[1:p])
             if (non_negative>0)
             {
                 L <- abs(L)
             }
         }




        lapla <-  init_lapla*matrix(1,nrow=l,ncol=p)

        Psi <- init_psi*XX


	cyc <- as.integer(cyc)
	nL <- as.integer(nL)
	lL <- as.integer(lL)
	bL <- as.integer(bL)
	alpha <- as.double(alpha)
	non_negative <- as.integer(non_negative)
	p <- as.integer(p)
	spz <- as.double(spz)
	lap <- as.double(lap)


	res <- .Call("fabics", X,Psi,L,lapla,cyc ,alpha,eps,spz,lap,nL,lL,bL,non_negative,PACKAGE="fabia")

        if (is.null(res))
        {
            return(new('Factorization', parameters=list(),n=1,p1=1,p2=1,l=1,center=as.vector(1),scaleData=as.vector(1),X=as.matrix(1),L=noL,Z=nZ,M=as.matrix(1),LZ=as.matrix(1),U=as.matrix(1),avini=as.vector(1),xavini=as.vector(1),ini=as.matrix(1),Psi=as.vector(1),lapla=as.matrix(1)))
        }

       # INI call for biclusters

        vz <- iin*apply(res$E_SX_n,1,function(x) sum(x^2))

        vz <- sqrt(vz+1e-10)

        ivz <- 1/vz

        if(length(ivz)==1) {
            nZ <- ivz*res$E_SX_n

            noL <- vz*res$L
        }
        else {
            nZ <- ivz*res$E_SX_n

            noL <- t(vz*t(res$L))
        }





        ini <- matrix(0,l,(p+1))
        avini <- as.vector(rep(0.0,(p+1)))
        xavini <- as.vector(rep(0.0,(l+1)))

        idp <- diag(p)
        ppL <- crossprod(noL,(1/res$Psi)*noL)


        for (j in 1:l){
            mat <- idp + ppL/res$lapla[j,]
            ini[j,1:p] <- log(diag(mat))
            s <- log(det(mat))
            ini[j,p+1] <- s
            xavini[j] <-s
        }
        for (i in 1:p){
            avini[i] <- sum(ini[,i])
        }

        ss <- sum(ini[,p+1])
        xavini[l+1] <- ss
        avini[p+1] <- ss

        Lz=noL%*%nZ


        rownames(noL) <- rowna
        colnames(noL) <- colnames(noL, do.NULL = FALSE, prefix = "bicluster")
        clnames <- colnames(noL)

        rownames(nZ) <- clnames
        colnames(nZ) <- colna

        M <- diag(p)
        rownames(M) <- clnames
        colnames(M) <- clnames

        rownames(Lz) <- rowna
        colnames(Lz) <- colna

        U <- X-Lz

        rownames(U) <- rowna
        colnames(U) <- colna

        if ((avini[p+1]>1e-8)&&(p>1)) {

            soo <- sort(avini[1:p], decreasing = TRUE,index.return=TRUE)

            avini[1:p] <- avini[soo$ix]
            noL <- noL[,soo$ix]
            nZ <- nZ[soo$ix,]
            M <- M[soo$ix,soo$ix]

        }

    return(new('Factorization', parameters=list("fabias",cyc,alpha,NULL,spz,p,NULL,NULL,random,scale,norm,center,lap,nL,lL,bL,non_negative),n=n,p1=p,p2=p,l=l,center=cent,scaleData=scaleData,X=X,L=noL,Z=nZ,M=M,LZ=Lz,U=U,avini=avini,xavini=xavini,ini=ini,Psi=res$Psi,lapla=res$lapla))


}


fabi <- function(X,p=5,alpha=0.1,cyc=500,spl=0,spz=0.5,center=2,norm=1,lap=1.0){

        if (missing(X)) {
            stop("Data matrix X missing. Stopped.")
        }

        if (!is.matrix(X)) {
            X <- as.matrix(X)
        }

        if (!is.numeric(X)) {
            stop("Data matrix X must be numeric. Stopped.")
        }


        l=ncol(X)
        n=nrow(X)

        if (p>min(l,n)) {
            stop("Too many biclusters. Stopped.")
        }

        rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "gene")
        colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
        rowna <- rownames(X)
        colna <- colnames(X)



       message("Running FABI on a ",n,"x",l," matrix with")
       message("   Number of biclusters ---------------- p: ",p)
       message("   Sparseness factor --------------- alpha: ",alpha)
       message("   Number of iterations -------------- cyc: ",cyc)
       message("   Loading prior parameter ----------- spl: ",spl)
       message("   Factor prior parameter ------------ spz: ",spz)
       if (center>0) {
            if (center<2) {
       message("   Centering ---------------------- center: 1 = mean")
            } else {
                if (center<3) {
       message("   Centering ---------------------- center: 2 = median")
                } else {
       message("   Centering ---------------------- center: ", center,"  = mode")
                }
            }
        } else {
       message("   Centering ---------------------- center: ", center," = no centering")
        }
       if (norm>0) {
           if (norm > 1) {
       message("   Scaling to variance one: --------- norm: ",norm," = Yes")
           } else {
       message("   Quantile scaling (0.75-0.25): ---- norm: ",norm," = Yes")
           }} else {
       message("   Scaling -------------------------- norm: ",norm," = No")
       }

        eps <- as.double(1e-3)
        eps1 <- as.double(1e-10)
        init_lapla <- 1.0
        init_psi <- 0.2


        iin <-  1.0/l


        cent <- as.vector(rep(0.0,n))
        if (center>0) {
            if (center<2) {
                cent <- apply(X, 1, mean)
            } else {
                if (center<3) {
                    cent <- apply(X, 1, median)
                } else {
                    cent <- estimateMode(X)
                }
            }
        }

        X <- X - cent

        XX <- as.vector(rep(1,n))

        if (norm>0) {
            if (norm>1)
            {
                scaleData <-  1/sqrt(iin*apply(X,1,function(x) sum(x^2))+0.001*XX)
            }
            else
            {
                scaleData <-  1/sqrt(apply(X,1,function(x) {quantile(x,0.75) - quantile(x,0.25)})+0.001*XX)
            }
            X <- scaleData*X
        } else {
            scaleData <- as.vector(rep(1.0,n))
        }




       lapla <-  init_lapla*matrix(1,nrow=l,ncol=p)

        Psi <- init_psi*XX


        eps2 <- 1e-1



        kvect <- as.vector(rep(1,p))
        nvect <- as.vector(rep(1,n))
        nk_one <- matrix(1,n,p)
        nk_zero <- matrix(0,n,p)
        kk_zero <- matrix(0,p,p)
        kk_one <- diag(p)
	epsv<-eps1*kvect
	epsn<-eps*nvect
        E_SX_n <- matrix(0,p,l)
        E_SSXX_n <- list()

        L <- (0.5*XX^(.5))*nk_one
	for (i in 1:cyc){
		LPsi<-diag(1/Psi)%*%L
		LPsiL<-crossprod(L,LPsi)
 		sum1<- nk_zero
		sum2<- eps2*kk_one
 		for (j in 1:l){
                        laj <- lapla[j,]
			tmp <- chol2inv(chol(LPsiL+diag(laj)))
			x_j <- as.vector(X[,j])
			e_sx_n <- as.vector(tcrossprod(tmp,LPsi)%*%x_j)
			#e_sx_n[which(e_sx_n<0)] <- 0
			e_ssxx_n <-  tmp + tcrossprod(e_sx_n,e_sx_n)
			sum1 <- sum1 +  tcrossprod(x_j,e_sx_n)
			sum2 <- sum2 + e_ssxx_n
			laj <- (epsv+diag(e_ssxx_n))^(-spz)
                        laj[which(laj<lap)] <- lap
                        lapla[j,] <- laj
		}
#                L <- (sum1 - alpha*sign(sum1))
#                L[which(abs(sum1)<alpha)] <- 0
                sll <- chol2inv(chol(sum2))
                L <- sum1%*%sll
                ddL <- alpha*Psi*sign(L)*abs(nk_one*eps+L)^{-spl}
                L <- L - ddL
                L[which(abs(L)<abs(ddL))] <- 0

                Psi <- epsn+abs(XX - diag(tcrossprod(L,sum1)/n))
                if (i %% 20==0) { print(i)}
	}

        LPsi<-diag(1/Psi)%*%L
        LPsiL<-crossprod(L,LPsi)
        for (j in 1:l){
            tmp <- chol2inv(chol(LPsiL+diag(lapla[j,])))
            x_j <- as.vector(X[,j])
            e_sx_n <- as.vector(tcrossprod(tmp,LPsi)%*%x_j)
            E_SX_n[,j] <- e_sx_n
            E_SSXX_n[[j]] <-  tmp + tcrossprod(e_sx_n,e_sx_n)
        }

        # INI call for biclusters

        vz <- iin*apply(E_SX_n,1,function(x) sum(x^2))

        vz <- sqrt(vz+1e-10)

        ivz <- 1/vz

        if(length(ivz)==1) {
            nZ <- ivz*E_SX_n

            noL <- vz*L
        }
        else {
            nZ <- ivz*E_SX_n

            noL <- t(vz*t(L))
        }



        ini <- matrix(0,l,(p+1))
        avini <- as.vector(rep(0.0,(p+1)))
        xavini <- as.vector(rep(0.0,(l+1)))

        idp <- diag(p)
        ppL <- crossprod(noL,(1/Psi)*noL)

        for (j in 1:l){
            mat <- idp + ppL/lapla[j,]
            ini[j,1:p] <- log(diag(mat))
            s <- log(det(mat))
            ini[j,p+1] <- s
            xavini[j] <-s
        }
        for (i in 1:p){
            avini[i] <- sum(ini[,i])
        }

        ss <- sum(ini[,p+1])
        xavini[l+1] <- ss
        avini[p+1] <- ss


        Lz=noL%*%nZ


        rownames(noL) <- rowna
        colnames(noL) <- colnames(noL, do.NULL = FALSE, prefix = "bicluster")
        clnames <- colnames(noL)

        rownames(nZ) <- clnames
        colnames(nZ) <- colna

        M <- diag(p)
        rownames(M) <- clnames
        colnames(M) <- clnames

        rownames(Lz) <- rowna
        colnames(Lz) <- colna

        U <- X-Lz

        rownames(U) <- rowna
        colnames(U) <- colna

        if ((avini[p+1]>1e-8)&&(p>1)) {

            soo <- sort(avini[1:p], decreasing = TRUE,index.return=TRUE)

            avini[1:p] <- avini[soo$ix]
            noL <- noL[,soo$ix]
            nZ <- nZ[soo$ix,]
            M <- M[soo$ix,soo$ix]

        }

        return(new('Factorization', parameters=list("fabi",cyc,alpha,spl,spz,p,NULL,NULL,NULL,scale,norm,center,NULL),n=n,p1=p,p2=p,l=l,center=cent,scaleData=scaleData,X=X,L=noL,Z=nZ,M=M,LZ=Lz,U=U,avini=avini,xavini=xavini,ini=ini,Psi=Psi,lapla=lapla))

}


fabiasp <- function(X,p=5,alpha=0.6,cyc=500,spz=0.5,center=2,norm=1,lap=1.0){


         if (missing(X)) {
            stop("Data matrix X missing. Stopped.")
        }

        if (!is.matrix(X)) {
            X <- as.matrix(X)
        }

        if (!is.numeric(X)) {
            stop("Data matrix X must be numeric. Stopped.")
        }


        l=ncol(X)
        n=nrow(X)

         if (p>min(l,n)) {
             stop("Too many biclusters. Stopped.")
         }

         rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "gene")
        colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
        rowna <- rownames(X)
        colna <- colnames(X)



       message("Running FABIASP on a ",n,"x",l," matrix with")
       message("   Number of biclusters ---------------- p: ",p)
       message("   Sparseness projection factor ---- alpha: ",alpha)
       message("   Number of iterations -------------- cyc: ",cyc)
       message("   Factor prior parameter ------------ spz: ",spz)
       if (center>0) {
            if (center<2) {
       message("   Centering ---------------------- center: 1 = mean")
            } else {
                if (center<3) {
       message("   Centering ---------------------- center: 2 = median")
                } else {
       message("   Centering ---------------------- center: ", center,"  = mode")
                }
            }
        } else {
       message("   Centering ---------------------- center: ", center," = no centering")
        }
       if (norm>0) {
           if (norm > 1) {
       message("   Scaling to variance one: --------- norm: ",norm," = Yes")
           } else {
       message("   Quantile scaling (0.75-0.25): ---- norm: ",norm," = Yes")
           }} else {
       message("   Scaling -------------------------- norm: ",norm," = No")
       }





        eps <- as.double(1e-3)
        eps1 <- as.double(1e-10)
        init_lapla <- 1.0
        init_psi <- 0.2


         iin <-  1.0/l

        cent <- as.vector(rep(0.0,n))
        if (center>0) {
            if (center<2) {
                cent <- apply(X, 1, mean)
            } else {
                if (center<3) {
                    cent <- apply(X, 1, median)
                } else {
                    cent <- estimateMode(X)
                }
            }
        }

        X <- X - cent

         XX <- as.vector(rep(1,n))

        if (norm>0) {
            if (norm>1)
            {
                scaleData <-  1/sqrt(iin*apply(X,1,function(x) sum(x^2))+0.001*XX)
            }
            else
            {
                scaleData <-  1/sqrt(apply(X,1,function(x) {quantile(x,0.75) - quantile(x,0.25)})+0.001*XX)
            }
            X <- scaleData*X
        } else {
            scaleData <- as.vector(rep(1.0,n))
        }



        lapla <-  init_lapla*matrix(1,nrow=l,ncol=p)

        Psi <- init_psi*XX


        eps2 <- 1e-1



        L1a <- sqrt(n)-(sqrt(n)-1)*alpha


        kvect <- as.vector(rep(1,p))
        nvect <- as.vector(rep(1,n))
        nk_one <- matrix(1,n,p)
        nk_zero <- matrix(0,n,p)
        kk_zero <- matrix(0,p,p)
        kk_one <- diag(p)
	epsv<-eps1*kvect
	epsn<-eps*nvect
        E_SX_n <- matrix(0,p,l)
        E_SSXX_n <- list()

        L <- (0.5*XX^(.5))*nk_one
        L <- L/sqrt(sum(L*L))

	for (i in 1:cyc){
		LPsi<-diag(1/Psi)%*%L
		LPsiL<-crossprod(L,LPsi)
 		sum1<- nk_zero
		sum2<- eps2*kk_one
 		for (j in 1:l){
                        laj <- lapla[j,]
			tmp <- chol2inv(chol(LPsiL+diag(laj)))
			x_j <- as.vector(X[,j])
			e_sx_n <- as.vector(tcrossprod(tmp,LPsi)%*%x_j)
			#e_sx_n[which(e_sx_n<0)] <- 0
			e_ssxx_n <-  tmp + tcrossprod(e_sx_n,e_sx_n)
			sum1 <- sum1 +  tcrossprod(x_j,e_sx_n)
			sum2 <- sum2 + e_ssxx_n
			laj <- (epsv+diag(e_ssxx_n))^(-spz)
                        laj[which(laj<lap)] <- lap
                        lapla[j,] <- laj
		}
                L <- sum1%*%chol2inv(chol(sum2))

                Psi <- epsn+abs(XX - diag(tcrossprod(L,sum1)/n))

                for (j in 1:p)
                {
                    L[,j] <- projFunc(s=as.vector(L[,j]),k1=L1a,k2=1.0)
                }

                if (i %% 20==0) { print(i)}
	}

        LPsi<-diag(1/Psi)%*%L
        LPsiL<-crossprod(L,LPsi)
        for (j in 1:l){
            tmp <- chol2inv(chol(LPsiL+diag(lapla[j,])))
            x_j <- as.vector(X[,j])
            e_sx_n <- as.vector(tcrossprod(tmp,LPsi)%*%x_j)
            E_SX_n[,j] <- e_sx_n
            E_SSXX_n[[j]] <-  tmp + tcrossprod(e_sx_n,e_sx_n)
        }


        # INI call for biclusters

        vz <- iin*apply(E_SX_n,1,function(x) sum(x^2))


        vz <- sqrt(vz+1e-10)

        ivz <- 1/vz

        if(length(ivz)==1) {
            nZ <- ivz*E_SX_n

            noL <- vz*L
        }
        else {
            nZ <- ivz*E_SX_n

            noL <- t(vz*t(L))
        }







        ini <- matrix(0,l,(p+1))
        avini <- as.vector(rep(0.0,(p+1)))
        xavini <- as.vector(rep(0.0,(l+1)))

        idp <- diag(p)
        ppL <- crossprod(noL,(1/Psi)*noL)

        for (j in 1:l){
            mat <- idp + ppL/lapla[j,]
            ini[j,1:p] <- log(diag(mat))
            s <- log(det(mat))
            ini[j,p+1] <- s
            xavini[j] <- s
        }
        for (i in 1:p){
            avini[i] <- sum(ini[,i])
        }


        ss <- sum(ini[,p+1])
        xavini[l+1] <- ss
        avini[p+1] <- ss



        Lz=noL%*%nZ

        rownames(noL) <- rowna
        colnames(noL) <- colnames(noL, do.NULL = FALSE, prefix = "bicluster")
        clnames <- colnames(noL)

        rownames(nZ) <- clnames
        colnames(nZ) <- colna

        M <- diag(p)
        rownames(M) <- clnames
        colnames(M) <- clnames

        rownames(Lz) <- rowna
        colnames(Lz) <- colna

        U <- X-Lz

        rownames(U) <- rowna
        colnames(U) <- colna

        if ((avini[p+1]>1e-8)&&(p>1)) {

            soo <- sort(avini[1:p], decreasing = TRUE,index.return=TRUE)

            avini[1:p] <- avini[soo$ix]
            noL <- noL[,soo$ix]
            nZ <- nZ[soo$ix,]
            M <- M[soo$ix,soo$ix]

        }


         return(new('Factorization', parameters=list("fabiasp",cyc,alpha,NULL,spz,p,NULL,NULL,NULL,NULL,norm,center,NULL),n=n,p1=p,p2=p,l=l,center=cent,scaleData=scaleData,X=X,L=noL,Z=nZ,M=M,LZ=Lz,U=U,avini=avini,xavini=xavini,ini=ini,Psi=Psi,lapla=lapla))

}



fabiaVersion <- function() {
  version <- packageDescription("fabia",fields="Version")
  message('\nFABIA Package Version ', version, '\n')
  message('\nCopyright (c) 2010 by Sepp Hochreiter\n\n')
}


estimateMode <- function(X,maxiter=50,tol=0.001,alpha=0.1,a1=4.0,G1=FALSE) {

        if (missing(X)) {
            stop("Data matrix X missing for mode. Stopped.")
        }

        if (!is.matrix(X)) {
            X <- as.matrix(X)
        }

        if (!is.numeric(X)) {
            stop("Data matrix X must be numeric for mode. Stopped.")
        }




    l=ncol(X)
    n=nrow(X)


    iin <-  1.0/l
    u <- apply(X, 1, median)
    xu <- X - u
    XX <-  sqrt(iin*apply(xu,1,function(x) sum(x^2))+0.001*as.vector(rep(1,n)))
    dXX <- 1/XX
    X <- dXX*X

    u <- apply(X, 1, median)
    gxu <- 10.0*tol
    iter <- 1
    alpha <- alpha/l

    xu <- X - u

    if (G1) {
        while (abs(gxu) > tol && iter < maxiter) {
            gxu <- apply(xu, 1, function(x) alpha*sum(tanh(a1 * x)))
            u <- u + gxu
            xu <- X - u
            iter <- iter + 1
        }
    } else {
        while (abs(gxu) > tol && iter < maxiter) {
            gxu <- apply(xu, 1, function(x) alpha*sum(x * exp(-a1*(x^2)/2)))
            u <- u + gxu
            xu <- X - u
            iter <- iter + 1
        }
    }

    u <- XX*u
    xu <- XX*xu

    return(list(u=u,xu=xu))
}

matrixImagePlot <- function(x,xLabels=NULL, yLabels=NULL, zlim=NULL, title=NULL){

        if (missing(x)) {
            stop("Data matrix x missing. Stopped.")
        }

        if (!is.matrix(x)) {
            x <- as.matrix(x)
        }

        if (!is.numeric(x)) {
            stop("Data matrix x must be numeric. Stopped.")
        }



     min <- min(x)
     max <- max(x)

      if( !is.null(zlim) ){
       min <- zlim[1]
       max <- zlim[2]
    }

    if( !is.null(yLabels) ){
       yLabels <- c(yLabels)
    } else {
       yLabels <- rownames(x)
   }

    if( !is.null(xLabels) ){
      xLabels <- c(xLabels)
    } else {
     xLabels <- colnames(x)
    }

     if( is.null(title) ){
     title <- c()
    }

# check for null values
if( is.null(xLabels) ){
   xLabels <- c(1:ncol(x))
}
if( is.null(yLabels) ){
   yLabels <- c(1:nrow(x))
}

#layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(7,1), heights=c(1,1))

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
 ColorRamp <- rgb( seq(0,1,length=256),  # Red
                   seq(0,1,length=256),  # Green
                   seq(1,0,length=256))  # Blue
 ColorLevels <- seq(min, max, length=length(ColorRamp))

 # Reverse Y axis
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

 # Data Map
# par(mar = c(3,5,2.5,2))
 par(mar = c(3,5,2.5,1))
 image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
 ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
 axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
 cex.axis=0.7)

 # Color Scale
# par(mar = c(3,2.5,2.5,2))
 par(mar = c(3,2,2.5,1))
 image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")

 layout(1)
}



plotBicluster <- function(r,p,opp=FALSE,zlim=NULL,title=NULL,which=c(1, 2)){



        if (missing(r)) {
            stop("r (= the result of extractBic) is missing. Stopped.")
        }
        if (missing(p)) {
            stop("The bicluster to plot is missing. Stopped.")
        }


        x <- r$X
        if (!is.matrix(x)) {
            x <- as.matrix(x)
        }

        if (!is.numeric(x)) {
            stop("Data matrix X must be numeric. Stopped.")
        }
        if (!is.numeric(r$np)) {
            stop("Number of biclusters must be numeric. Stopped.")
        }


        if (p>r$np) {
            stop("Bicluster number is too large. Stopped.")
        }

        if ((nrow(x)<2)||(ncol(x)<2)) {
            stop("Data matrix X too small or missing. Stopped.")
        }


     devAskNewPage(ask = FALSE)

   if (length(which) > 1 && dev.interactive()) {
        devAskNewPage(ask = TRUE)
    }

   showf <- c(FALSE, FALSE)
    showf[which] <- TRUE




     n <- nrow(x)

     np <- r$np

     l <- ncol(x)

     if (!opp) {
         if (r$bic[[p]][1]<1) {
            stop("No genes in bicluster. Stopped.")
         }
         if (r$bic[[p]][2]<1) {
            stop("No samples in bicluster. Stopped.")
         }

        observations <- unlist(r$bic[p,3])
         samples <- unlist(r$bic[p,5])
     } else {
         if (r$bicopp[[p]][1]<1) {
            stop("No genes in bicluster. Stopped.")
         }
         if (r$bicopp[[p]][2]<1) {
            stop("No samples in bicluster. Stopped.")
         }
         observations <- unlist(r$bicopp[p,3])
         samples <- unlist(r$bicopp[p,5])
     }

     if (length(observations)<1) {
         stop("No genes in bicluster. Stopped.")
     }

     if (length(samples)<1) {
         stop("No samples in bicluster. Stopped.")
     }

     min <- min(x)
     max <- max(x)

   if( !is.null(zlim) ){
       min <- zlim[1]
       max <- zlim[2]
    }

    yLabels <- rownames(x)
    xLabels <- colnames(x)


     if( is.null(title) ){
     title <- c()
    }



 # Reverse Y axis
 reverse <- n : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

     obs <- match(observations,yLabels)
     samp <- match(samples,xLabels)

    ybl <- length(obs)
    if (ybl == 0) {
      obs <- match(observations,as.character(n:1))
      ybl <- length(obs)

    }
    if (ybl == 0) {
      obs <- match(observations,c(n:1))
      ybl <- length(obs)

    }

    xbl <- length(samp)
    if (xbl == 0) {
      samp <- match(samples,as.character(1:l))
      xbl <- length(samp)

    }
    if (xbl == 0) {
      samp <- match(samples,c(1:l))
      xbl <- length(samp)

    }



    pmL <- diag(n)
    for (i in 1:ybl){
        ii <- n-i+1
        pmL[ii,ii] <- 0
        pmL[ii,obs[i]] <- 1
        pmL[obs[i],obs[i]] <- 0
        pmL[obs[i],ii] <- 1
        tmp <- yLabels[ii]
        yLabels[ii] <- yLabels[obs[i]]
        yLabels[obs[i]] <- tmp
    }

    pmZ <- diag(l)
    for (i in 1:xbl){
        pmZ[i,i] <- 0
        pmZ[i,samp[i]] <- 1
        pmZ[samp[i],samp[i]] <- 0
        pmZ[samp[i],i] <- 1
        tmp <- xLabels[i]
        xLabels[i] <- xLabels[samp[i]]
        xLabels[samp[i]] <- tmp
   }

    x <- pmL%*%x%*%pmZ

    rownames(x) <-  yLabels
    colnames(x) <-  xLabels


    if (n>100)
    {
        lll <- 12

        yl <- c(1,round(n*(1:lll)/lll))

        yll <- rep("",n)

        yll[yl] <- rownames(x)[yl]

        yLabels <- yll
    }
    for (i in 1:ybl){
        ii <- n-i+1
        yLabels[ii] <- rownames(x)[ii]
    }

if (showf[1]){

 layout(matrix(data=c(2,1), nrow=1, ncol=2), widths=c(7,1), heights=c(1,1))

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
 ColorRamp <- rgb( seq(0,1,length=256),  # Red
                   seq(0,1,length=256),  # Green
                   seq(1,0,length=256))  # Blue
 ColorLevels <- seq(min, max, length=length(ColorRamp))

 # Color Scale
 par(mar = c(3,2,2.5,1))
 image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")




 # Data Map
 par(mar = c(3,5,2.5,1))
 image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
 ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
 axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
 cex.axis=0.7)

#rect(xleft, ybottom, xright, ytop, density = NULL, angle = 45,
#     col = NA, border = NULL, lty = par("lty"), lwd = par("lwd"),
#     ...)

rect(xleft=0.6, ybottom= n-ybl, xright=xbl+1, ytop=n+0.3,lwd=3,border="red")

layout(1)

}


if (showf[2]){

layout(matrix(data=c(2,1), nrow=1, ncol=2), widths=c(7,1), heights=c(1,1))

 ColorRamp <- rgb( seq(0,1,length=256),  # Red
                   seq(0,1,length=256),  # Green
                   seq(1,0,length=256))  # Blue
 ColorLevels <- seq(min, max, length=length(ColorRamp))

# Color Scale
 par(mar = c(3,2,2.5,1))
 image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")


 # Data Map
 par(mar = c(3,5,2.5,1))
 image(1:xbl, 1:ybl, t(x[(n-ybl+1):n,1:xbl]), col=ColorRamp, xlab="",
 ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
axis(BELOW<-1, at=1:xbl, labels=xLabels[1:xbl], cex.axis=0.7)
 axis(LEFT <-2, at=1:ybl, labels=yLabels[(n-ybl+1):n], las= HORIZONTAL<-1,
 cex.axis=0.7)

 layout(1)

}
devAskNewPage(ask = FALSE)
}




#########################
##Data
makeFabiaData <- function(n,l,p,f1,f2,of1,of2,sd_noise,sd_z_noise,mean_z,sd_z,sd_l_noise,mean_l,sd_l){

    if(!is.numeric(n)) {
        stop("n must be numeric")
    }
    if(!is.numeric(l)) {
        stop("l must be numeric")
    }

    if(!is.numeric(p)) {
        stop("p must be numeric")
    }
    if(!is.numeric(f1)) {
        stop("f1 must be numeric")
    }



    za <- floor(runif(p)*(l/f1)+of1)
    la <- floor(runif(p)*(n/f2)+of2)
    ZC <- list()
    LC <- list()

    X <- matrix(rnorm(l*n,mean = 0, sd = sd_noise),n,l)
    Y <- matrix(0,n,l)
    rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "gene")
    colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
    rownames(Y) <- rownames(X)
    colnames(Y) <- colnames(X)

for (i in 1:p){

        zj <- rnorm(l,mean = 0, sd = sd_z_noise)
        z0 <- rep(0,l)
        zf <- floor(runif(za[i])*l)
        ZC[[i]] <- zf
        z1 <- rnorm(za[i],mean = mean_z, sd = sd_z)
        zj[zf] <- z1
        z0[zf] <- z1

        lj <- rnorm(n,mean = 0, sd = sd_l_noise)
        l0 <- rep(0,n)
        lf <- floor(runif(la[i])*n)
        LC[[i]] <- lf
        l1 <- rnorm(la[i],mean = mean_l, sd = sd_l)
        lsign <- 2*round(runif(la[i]))-1
        l1 <- l1*lsign
        lj[lf] <- l1
        l0[lf] <- l1

        X <- X + lj%*%t(zj)
        Y <- Y + l0%*%t(z0)



    }


return(list(X=X,Y=Y,ZC=ZC,LC=LC))
}

#########################
##Blocks
makeFabiaDataBlocks <- function(n,l,p,f1,f2,of1,of2,sd_noise,sd_z_noise,mean_z,sd_z,sd_l_noise,mean_l,sd_l){

    if(!is.numeric(n)) {
        stop("n must be numeric")
    }
    if(!is.numeric(l)) {
        stop("l must be numeric")
    }

    if(!is.numeric(p)) {
         stop("p must be numeric")
    }
    if(!is.numeric(f1)) {
        stop("f1 must be numeric")
    }


    za <- floor(runif(p)*(l/f1)+of1)
    la <- floor(runif(p)*(n/f2)+of2)
    ZC <- list()
    LC <- list()


    X <- matrix(rnorm(l*n,mean = 0, sd = sd_noise),n,l)
    Y <- matrix(0,n,l)

    rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "gene")
    colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
    rownames(Y) <- rownames(X)
    colnames(Y) <- colnames(X)

   for (i in 1:p){


        zj <- rnorm(l,mean = 0, sd = sd_z_noise)
        z0 <- as.vector(rep(0,l))

        tf <- floor(runif(1)*(l-za[i]))
        zf <- tf:(tf+za[i]-1)
        ZC[[i]] <- zf
        z1 <- rnorm(za[i],mean = mean_z, sd = sd_z)
        zj[zf] <- z1
        z0[zf] <- z1

        lj <- rnorm(n,mean = 0, sd = sd_l_noise)
        l0 <- as.vector(rep(0,n))
        uf <- floor(runif(1)*(n-la[i]))
        lf <- uf:(uf+la[i]-1)
        LC[[i]] <- lf
        l1 <- rnorm(la[i],mean = mean_l, sd = sd_l)
        lsign <- 2*round(runif(la[i]))-1
        l1 <- l1*lsign
        lj[lf] <- l1
        l0[lf] <- l1

        X <- X + lj%*%t(zj)
        Y <- Y + l0%*%t(z0)


    }
    return(list(X=X,Y=Y,ZC=ZC,LC=LC))
}

#########################
##Data
makeFabiaDataPos <- function(n,l,p,f1,f2,of1,of2,sd_noise,sd_z_noise,mean_z,sd_z,sd_l_noise,mean_l,sd_l){

    if(!is.numeric(n)) {
        stop("n must be numeric")
    }
    if(!is.numeric(l)) {
        stop("l must be numeric")
    }

    if(!is.numeric(p)) {
        stop("p must be numeric")
    }
    if(!is.numeric(f1)) {
        stop("f1 must be numeric")
    }



    za <- floor(runif(p)*(l/f1)+of1)
    la <- floor(runif(p)*(n/f2)+of2)
    ZC <- list()
    LC <- list()

    X <- matrix(rnorm(l*n,mean = 0, sd = sd_noise),n,l)
    Y <- matrix(0,n,l)

    rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "gene")
    colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
    rownames(Y) <- rownames(X)
    colnames(Y) <- colnames(X)

   for (i in 1:p){

        zj <- rnorm(l,mean = 0, sd = sd_z_noise)
        z0 <- rep(0,l)
        zf <- floor(runif(za[i])*l)
        ZC[[i]] <- zf
        z1 <- rnorm(za[i],mean = mean_z, sd = sd_z)
        zj[zf] <- z1
        z0[zf] <- z1

        lj <- rnorm(n,mean = 0, sd = sd_l_noise)
        l0 <- rep(0,n)
        lf <- floor(runif(la[i])*n)
        LC[[i]] <- lf
        l1 <- rnorm(la[i],mean = mean_l, sd = sd_l)
        lj[lf] <- l1
        l0[lf] <- l1

        X <- X + lj%*%t(zj)
        Y <- Y + l0%*%t(z0)


    }
    return(list(X=X,Y=Y,ZC=ZC,LC=LC))
}

#########################
##Blocks
makeFabiaDataBlocksPos<- function(n,l,p,f1,f2,of1,of2,sd_noise,sd_z_noise,mean_z,sd_z,sd_l_noise,mean_l,sd_l){

    if(!is.numeric(n)) {
        stop("n must be numeric")
    }
    if(!is.numeric(l)) {
        stop("l must be numeric")
    }

    if(!is.numeric(p)) {
         stop("p must be numeric")
    }
    if(!is.numeric(f1)) {
        stop("f1 must be numeric")
    }


    za <- floor(runif(p)*(l/f1)+of1)
    la <- floor(runif(p)*(n/f2)+of2)
    ZC <- list()
    LC <- list()


    X <- matrix(rnorm(l*n,mean = 0, sd = sd_noise),n,l)
    Y <- matrix(0,n,l)


    rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "gene")
    colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
    rownames(Y) <- rownames(X)
    colnames(Y) <- colnames(X)

    for (i in 1:p){


        zj <- rnorm(l,mean = 0, sd = sd_z_noise)
        z0 <- as.vector(rep(0,l))

        tf <- floor(runif(1)*(l-za[i]))
        zf <- tf:(tf+za[i]-1)
        ZC[[i]] <- zf
        z1 <- rnorm(za[i],mean = mean_z, sd = sd_z)
        zj[zf] <- z1
        z0[zf] <- z1

        lj <- rnorm(n,mean = 0, sd = sd_l_noise)
        l0 <- as.vector(rep(0,n))
        uf <- floor(runif(1)*(n-la[i]))
        lf <- uf:(uf+la[i]-1)
        LC[[i]] <- lf
        l1 <- rnorm(la[i],mean = mean_l, sd = sd_l)
        lj[lf] <- l1
        l0[lf] <- l1

        X <- X + lj%*%t(zj)
        Y <- Y + l0%*%t(z0)


    }
    return(list(X=X,Y=Y,ZC=ZC,LC=LC))
}


# mfsc - matrix factorization with sparseness constraints
#
# Written by Sepp Hochreiter according to
# Patrik O. Hoyer
# 'Non-negative matrix factorization with sparseness constraints'
# Journal of Machine Learning Research  5:1457-1469, 2004.
#
#
#
# SYNTAX:
# res <- mfsc(X,p=5,cyc=100,sL=0.6,sZ=0.6,center=2,norm=1)
# res$L
# res$Z
#
# INPUTS:
# X          - data matrix
# p       - number of components (inner dimension of factorization)
# sL         - sparseness of L, in [0,1]
# sZ         - sparseness of Z, in [0,1]
# cyc    - maximal number of iterations
#
# Note: Sparseness is measured on the scale [0,1] where 0 means
# completely distributed and 1 means ultimate sparseness.
#
# NOTE: There is NO CONVERGENCE CRITERION.
#


mfsc <- function(X,p=5,cyc=100,sL=0.6,sZ=0.6,center=2,norm=1) {

        if (missing(X)) {
            stop("Data matrix X missing. Stopped.")
        }

        if (!is.matrix(X)) {
            X <- as.matrix(X)
        }

        if (!is.numeric(X)) {
            stop("Data matrix X must be numeric. Stopped.")
        }

        l=ncol(X)
        n=nrow(X)

        if (p>min(l,n)) {
            stop("Too many biclusters. Stopped.")
        }


        rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "gene")
        colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
        rowna <- rownames(X)
        colna <- colnames(X)


       message("Running MFSC on a ",n,"x",l," matrix with")
       message("   Number of biclusters ---------------- p: ",p)
       message("   Number of maximal iterations ------ cyc: ",cyc)
       message("   Loading sparseness projection ------ sL: ",sL)
       message("   Factor sparseness projection ------- sZ: ",sZ)
       if (center>0) {
            if (center<2) {
       message("   Centering ---------------------- center: 1 = mean")
            } else {
                if (center<3) {
       message("   Centering ---------------------- center: 2 = median")
                } else {
       message("   Centering ---------------------- center: ", center,"  = mode")
                }
            }
        } else {
       message("   Centering ---------------------- center: ", center," = no centering")
        }
       if (norm>0) {
           if (norm > 1) {
       message("   Scaling to variance one: --------- norm: ",norm," = Yes")
           } else {
       message("   Quantile scaling (0.75-0.25): ---- norm: ",norm," = Yes")
           }} else {
       message("   Scaling -------------------------- norm: ",norm," = No")
       }

        iin <-  1.0/l

        ini <- matrix(1,l,(p+1))
        avini <- as.vector(rep(1.0,(p+1)))
        xavini <- as.vector(rep(1.0,(l+1)))

        cent <- as.vector(rep(0.0,n))
        if (center>0) {
            if (center<2) {
                cent <- apply(X, 1, mean)
            } else {
                if (center<3) {
                    cent <- apply(X, 1, median)
                } else {
                    cent <- estimateMode(X)
                }
            }
        }

        X <- X - cent

        XX <- as.vector(rep(1,n))

        if (norm>0) {
            if (norm>1)
            {
                scaleData <-  1/sqrt(iin*apply(X,1,function(x) sum(x^2))+0.001*XX)
            }
            else
            {
                scaleData <-  1/sqrt(apply(X,1,function(x) {quantile(x,0.75) - quantile(x,0.25)})+0.001*XX)
            }
            X <- scaleData*X
        } else {
            scaleData <- as.vector(rep(1.0,n))
        }


# Create initial matrices
L <- matrix (rnorm(n*p),nrow=n,ncol=p)
Z <- matrix (rnorm(l*p),nrow=p,ncol=l)



# Make initial matrices have correct sparseness
L1a <- sqrt(n)-(sqrt(n)-1)*sL
L1s <- sqrt(l)-(sqrt(l)-1)*sZ

for (i in 1:p)
  {
    L[,i] <- projFunc(s=as.vector(L[,i]),k1=L1a,k2=1.0)
    Z[i,] <- projFunc(s=as.vector(Z[i,]),k1=L1s,k2=1.0)
  }


# Calculate initial objective
objhistory <- sum((X-L%*%Z)^2)



# Initial stepsizes
stepsizeL <- 1.0
stepsizeZ <- 1.0



# Start iteration
for (i in 1:cyc)
  {


    # Save old values
    Lold <- L
    Zold <- Z

    # ----- Update Z ---------------------------------------


	# Gradient for Z
	dZ <- t(L)%*%(L%*%Z-X)

        begobj <- objhistory

        newobj <- 2.0*begobj




     # Make sure we decrease the objective!
	while (newobj>begobj) {


         # Take step in direction of negative gradient, and project
          Znew <- Z - stepsizeZ*dZ
	# Make initial matrices have correct sparseness
          for (j in 1:p)
            {
              Znew[j,] <- projFunc(s=as.vector(Znew[j,]),k1=L1s,k2=1.0)
            }


         # Calculate new objective
	    newobj <- sum((X-L%*%Znew)^2)

	    # If the objective decreased, we can continue...
            stepsizeZ <- stepsizeZ/2.0
            if (stepsizeZ<1e-10) {

                rownames(L) <- rowna
                colnames(L) <- colnames(L, do.NULL = FALSE, prefix = "bicluster")
                clnames <- colnames(L)

                rownames(Z) <- clnames
                colnames(Z) <- colna

                M <- diag(p)
                rownames(M) <- clnames
                colnames(M) <- clnames
                LZ <- L%*%Z

                rownames(LZ) <- rowna
                colnames(LZ) <- colna

                U <- X-LZ

                rownames(U) <- rowna
                colnames(U) <- colna

                for (i in 1:p){
                    avini[i] <- sum(L[,i]^2)*sum(Z[i,]^2)
                }
                for (i in 1:l){
                    xavini[i] <-  sum(L[i,]^2)
                }

                avini[p+1] <- sum(avini[1:p])
                xavini[l+1] <- sum(xavini[1:l])


                if (avini[p+1]>1e-8) {

                    soo <- sort(avini[1:p], decreasing = TRUE,index.return=TRUE)

                    avini[1:p] <- avini[soo$ix]
                    L <- L[,soo$ix]
                    Z <- Z[soo$ix,]
                    M <- M[soo$ix,soo$ix]

                }

                return(new('Factorization', parameters=list("mfsc",cyc,NULL,NULL,NULL,p,sL,sZ,NULL,NULL,norm,center,NULL),n=n,p1=p,p2=p,l=l,center=cent,scaleData=scaleData,X=X,L=L,Z=Z,M=M,LZ=LZ,U=U,avini=avini,xavini=xavini,ini=ini,Psi=as.vector(1),lapla=as.matrix(1)))
            }
        }

	# Slightly increase the stepsize
	stepsizeZ <- stepsizeZ*1.2
	Z <- Znew



    # ----- Update L ---------------------------------------


	# Gradient for L
	dL <- (L%*%Z-X)%*%t(Z)
	begobj <- sum((X-L%*%Z)^2)

        newobj <- 2*begobj

	# Make sure we decrease the objective!
	while (newobj>begobj) {



	    # Take step in direction of negative gradient, and project
	    Lnew <- L - stepsizeL*dL
	    norms <- sqrt(sum(Lnew^2))
            for (j in 1:p)
              {
                Lnew[,j] <- projFunc(s=as.vector(Lnew[,j]),k1=L1a,k2=1.0)
              }

	    # Calculate new objective
	    newobj <- sum((X-Lnew%*%Z)^2)

	    stepsizeL <- stepsizeL/2

            if (stepsizeZ<1e-10) {

                rownames(L) <- rowna
                colnames(L) <- colnames(L, do.NULL = FALSE, prefix = "bicluster")
                clnames <- colnames(L)

                rownames(Z) <- clnames
                colnames(Z) <- colna

                M <- diag(p)
                rownames(M) <- clnames
                colnames(M) <- clnames
                LZ <- L%*%Z

                rownames(LZ) <- rowna
                colnames(LZ) <- colna

                U <- X-LZ

                rownames(U) <- rowna
                colnames(U) <- colna

                for (i in 1:p){
                    avini[i] <- sum(L[,i]^2)*sum(Z[i,]^2)
                }
                for (i in 1:l){
                    xavini[i] <-  sum(L[i,]^2)
                }

                avini[p+1] <- sum(avini[1:p])
                xavini[l+1] <- sum(xavini[1:l])


                if (avini[p+1]>1e-8) {

                    soo <- sort(avini[1:p], decreasing = TRUE,index.return=TRUE)

                    avini[1:p] <- avini[soo$ix]
                    L <- L[,soo$ix]
                    Z <- Z[soo$ix,]
                    M <- M[soo$ix,soo$ix]

                }


     return(new('Factorization', parameters=list("mfsc",cyc,NULL,NULL,NULL,p,sL,sZ,NULL,NULL,norm,center,NULL),n=n,p1=p,p2=p,l=l,center=cent,scaleData=scaleData,X=X,L=L,Z=Z,M=M,LZ=LZ,U=U,avini=avini,xavini=xavini,ini=ini,Psi=as.vector(1),lapla=as.matrix(1)))
            }


          }

	# Slightly increase the stepsize
	stepsizeL <- stepsizeL*1.2
	L <- Lnew


    # Calculate objective
    newobj <- 0.5*sum(sum((X-L%*%Z)^2))
    objhistory <- newobj



  }
}




# Solves the following problem:
# Given a vector s, find the vector v having sum(abs(v))=k1
# and sum(v^2)=k2 which is closest to s in the euclidian sense.


projFunc <- function(s, k1, k2) {


  N <- length(s)

  isneg <- as.vector(rep(0,N))
  isneg[which(s<0)] <- 1
  ones <- as.vector(rep(1,N))

  s <- abs(s)


# Start by projecting the point to the sum constraint hyperplane
  v <- s + (k1-sum(s))/N

# Initialize zerocoeff (initially, no elements are assumed zero)
  # zerocoeff <- which(v<=0)
  # v[zerocoeff] <- 0

  zerocoeff <- vector(mode = "integer", length = 0)


  j <- 0
  ende <- 0
  while (ende<1) {

    # This does the proposed projection operator
    midpoint <- ones*k1/(N-length(zerocoeff))
    midpoint[zerocoeff] <- 0
    w <- v-midpoint
    a <- sum(w^2)
    b <- 2*crossprod(w,v)
    c <- sum(v^2)-k2
    t <- b^2-4*a*c
    if (t<0) {t <- 0}
    if (a<1e-10) {a <- 1e-10}

    alphap <- (-b+sqrt(t))/(2*a)
    v <- alphap*w + v

    if ((length(which(v<0))==0)&&(j>1)) {
      ende <- 2
    }
    else {

      j <- j+1

    # Set negs to zero, subtract appropriate amount from rest
      zerocoeff <- which(v<=0)
      v[zerocoeff] <- 0
      tempsum <- sum(v)
      tta <- N - length(zerocoeff)
      if (tta>0)  {
       v <- v + (k1-tempsum)/tta
       }
      v[zerocoeff] <- 0
     }
  }

  v <- (-2*isneg + ones)*v

  return(v)

}


##########################################################
extractPlot <- function(fact,thresZ=0.5,ti="FABIA",thresL=NULL,Y=NULL,which=c(1,2,3,4,5,6,7,8)){

        if (missing(fact)) {
            stop("Object fact of class Factorization is missing. Stopped.")
        }

        if (is(fact) != "Factorization") {
            stop("Object fact is not of the class Factorization. Stopped.")
        }


    devAskNewPage(ask = FALSE)
    if (length(which) > 1 && dev.interactive()) {
      devAskNewPage(ask = TRUE)
      if ((length(which) == 2) && (is.null(Y)) && (which[1]==1)) {
          devAskNewPage(ask = FALSE)
      }
    }



    showf <- c(FALSE, FALSE,FALSE, FALSE,FALSE, FALSE,FALSE, FALSE)
    showf[which] <- TRUE





    X <- as.matrix(fact@X)

    misX <- FALSE
    if (!is.matrix(X)||(nrow(X)<2)||(ncol(X)<2)) {
     showf[2] <- FALSE
     showf[4] <- FALSE
     showf[8] <- FALSE
     misX <- TRUE
    }

    noL <- as.matrix(fact@L)
    nZ <- as.matrix(fact@Z)

    l=ncol(nZ)
    n=nrow(noL)


    p <- ncol(noL)




    lll <- 12

    yl <- c(1,round(n*(1:lll)/lll))

    yll <- rep("",n)

    if (misX) {
        yll[yl] <- as.character(yl)
    } else {
        yll[yl] <- rownames(X)[yl]
    }
    LZ <- noL%*%nZ

    tt <- paste("(",as.character(n)," genes, ",as.character(l)," samples, ",as.character(p)," biclusters )")

    if (!is.null(Y) && showf[1]){
      matrixImagePlot(Y,title=paste(ti,": noise free data\n",tt,sep=""),yLabels= yll)
    }

    if (showf[2]){
        matrixImagePlot(X,title=paste(ti,": data\n",tt,sep=""),yLabels= yll)
    }

    if (showf[3]){
        matrixImagePlot(LZ,title=paste(ti,": reconstructed data\n",tt,sep=""),yLabels= yll)
    }
    if (showf[4]){
        matrixImagePlot(LZ-X,title=paste(ti,": error\n",tt,sep=""),yLabels= yll)
    }



    if (showf[5]){
        matrixImagePlot(abs(noL),title=paste(ti,": absolute loadings\n",tt,sep=""),yLabels= yll)
    }
    if (showf[6]){
        matrixImagePlot(abs(nZ),title=paste(ti,": absolute factors\n",tt,sep=""))
    }




    hnzc <- kmeans(t(nZ),p, iter.max = 20, nstart = 10)

    hnzs <- sort(hnzc$cluster,index.return = TRUE)

    if (misX) {
        xLabels <- as.character(1:l)
    } else {
        xLabels <- colnames(X)
    }

    pmZ <- matrix(0,l,l)
    for (i in 1:l){
         pmZ[hnzs$ix[i],i] <- 1
         tmp <- xLabels[i]
         xLabels[i] <- xLabels[hnzs$ix[i]]
         xLabels[hnzs$ix[i]] <- tmp
    }




    hnLc <- kmeans(noL,p, iter.max = 20, nstart = 10)

    hnLs <- sort(hnLc$cluster,index.return = TRUE)

    biclustx <- list()
    biclusty <- list()
    for (i in 1:p){

        biclustx[[i]] <- hnLs$ix[which(hnLs$x==i)]
        biclusty[[i]] <- hnzs$ix[which(hnzs$x==i)]

    }

    biclust <- cbind(biclustx,biclusty)

    pmL <- matrix(0,n,n)
    if (misX) {
        yLabels <- as.character(1:n)
    } else {
        yLabels <- rownames(X)
    }
    for (i in 1:n){
         pmL[i,hnLs$ix[i]] <- 1
         tmp <- yLabels[i]
         yLabels[i] <- yLabels[hnLs$ix[i]]
         yLabels[hnLs$ix[i]] <- tmp
     }


    lll <- 12

    yl <- c(1,round(n*(1:lll)/lll))

    yll <- rep("",n)

    yll[yl] <- yLabels[yl]

    reconstr <- pmL%*%noL%*%nZ%*%pmZ

    rownames(reconstr) <- yLabels
    colnames(reconstr) <- xLabels

    if (showf[7]){
        matrixImagePlot(reconstr,title=paste(ti,": reconstructed matrix sorted\n",tt,sep=""), yLabels= yll)
    }

   if (!misX) {
       reord <- pmL%*%as.matrix(X)%*%pmZ
       rownames(reord) <- yLabels
       colnames(reord) <- xLabels

       if (showf[8]){
           matrixImagePlot(reord,title=paste(ti,": original matrix sorted\n",tt,sep=""), yLabels= yll)
       }
   } else {
       reord <- NULL
   }

   devAskNewPage(ask = FALSE)

    return;

}


##########################################################
extractBic <- function(fact,thresZ=0.5,thresL=NULL)
{

    if (missing(fact)) {
        stop("Object fact of class Factorization is missing. Stopped.")
    }

    if (is(fact) != "Factorization") {
        stop("Object fact is not of the class Factorization. Stopped.")
    }



    noL <- as.matrix(fact@L)
    nZ <- as.matrix(fact@Z)
    X <- as.matrix(fact@X)

    n <- nrow(noL)

    p <- ncol(noL)

    l <- ncol(nZ)



    binn <- list()
    binp <- list()
    bixv <- list()
    bixn <- list()
    numng <- list()

    numnp <- list()
    biyp <- list()
    biypv <- list()
    biypn <- list()

    numnn <- list()
    biyn <- list()
    biynv <- list()
    biynn <- list()
    if (is.null(rownames(X))||(length(rownames(X))<2))
    {
        gene_names <- as.character(1:n)
    } else {
        gene_names <- rownames(X)
    }
    rownames(noL) <- gene_names

    if (is.null(colnames(X))||(length(colnames(X))<2))
    {
       sample_names <- as.character(1:l)
    } else {
        sample_names <- colnames(X)
    }
    colnames(nZ) <- sample_names

    if (is.null(thresL)) {
        mom <- 0
        for (i in 1:p) {
            tmom <- as.numeric(tcrossprod(noL[,i],nZ[i,]))
            mom <- mom + sum(tmom^2)
        }
        mom <- mom/(length(tmom)*p)

        thresL <- sqrt(mom)/thresZ
    }


    for (i in 1:p){

        snl <- sort(abs(noL[,i]),decreasing = TRUE,index.return = TRUE)
        gene_namest <- gene_names[snl$ix]
        sL <- noL[snl$ix,i]


        numng[[i]] <- as.vector(which(abs(noL[,i])>thresL))
        cho <- which(abs(sL)>thresL)
        bixv[[i]] <- sL[cho]
        bixn[[i]] <- gene_namest[cho]


        snz <- sort(abs(nZ[i,]),decreasing = TRUE,index.return = TRUE)
        sample_namest <- sample_names[snz$ix]
        sz <- nZ[i,snz$ix]

        chop <- which(sz>thresZ)
        ss2 <- sum(abs(sz[chop]))
        chon <- which(sz< -thresZ)
        ss3 <- sum(abs(sz[chon]))

        if (ss2>=ss3) {
            binp[[i]] <- c(length(cho),length(chop))
            numnp[[i]] <-  as.vector(which(nZ[i,]>thresZ))
            biypv[[i]] <- sz[chop]
            biypn[[i]] <-  sample_namest[chop]

            binn[[i]] <- c(length(cho),length(chon))
            numnn[[i]] <-  as.vector(which(nZ[i,]< -thresZ))
            biynv[[i]] <- sz[chon]
            biynn[[i]] <-  sample_namest[chon]
        } else {
            binn[[i]]    <- c(length(cho),length(chop))
            numnn[[i]]   <-  as.vector(which(nZ[i,]>thresZ))
            biynv[[i]]   <- sz[chop]
            biynn[[i]]   <-  sample_namest[chop]

            binp[[i]]    <- c(length(cho),length(chon))
            numnp[[i]]   <-  as.vector(which(nZ[i,]< -thresZ))
            biypv[[i]]   <- sz[chon]
            biypn[[i]]   <-  sample_namest[chon]

        }





    }

    bic <- cbind(binp,bixv,bixn,biypv,biypn)
    numn <- cbind(numng,numnp)

    bicopp <- cbind(binn,bixv,bixn,biynv,biynn)
    numnopp <- cbind(numng,numnn)




    return(list(bic=bic,numn=numn,bicopp=bicopp,numnopp=numnopp,X=X,np=p))

}





nmfdiv <- function(X,p=5,cyc=100) {


        if (missing(X)) {
            stop("Data matrix X missing. Stopped.")
        }

        if (!is.matrix(X)) {
            X <- as.matrix(X)
        }

        if (!is.numeric(X)) {
            stop("Data matrix X must be numeric. Stopped.")
        }

  # Check that we have non-negative data
if (min(X)<0) {
    stop("Negative values in matrix X. Stopped.")
}

# Dimensions
n <- nrow(X)
l <- ncol(X)

        rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "gene")
        colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
        rowna <- rownames(X)
        colna <- colnames(X)


       message("Running NMFDIV on a ",n,"x",l," matrix with")
       message("   Number of biclusters ---------------- p: ",p)
       message("   Number of maximal iterations ------ cyc: ",cyc)



# Globally rescale data to avoid potential overflow/underflow
X <- X/max(X)

cent <- as.vector(rep(0.0,n))
scaleData <- as.vector(rep(max(X),n))

# Create initial matrices
L <- matrix (abs(rnorm(n*p)),nrow=n,ncol=p)
Z <- matrix (abs(rnorm(l*p)),nrow=p,ncol=l)
#Z <- Z/(sqrt(rowsum(Z*Z,1:n))*as.vector(rep(1,l)))

# Calculate initial objective
objhistory <- 0.5*sum(sum((X-L%*%Z)^2))


ini <- matrix(1,l,(p+1))
avini <- as.vector(rep(1.0,(p+1)))
xavini <- as.vector(rep(1.0,(l+1)))



# Start iteration
for (i in 1:cyc)
  {

    # Show progress
    print(i)
    print(objhistory)

    # Save old values
    Lold <- L
    Zold <- Z


    # Compute new L and Z (Lee and Seung; NIPS*2000)

    #L <- L*((X/(L%*%Z + 1e-9))%*%t(Z))/rowSums(Z)
    #Z <- Z*(t(L)%*%(X/(L%*%Z + 1e-9)))/colSums(L)

    L <- L*((X/(L%*%Z + 1e-9))%*%t(Z))/tcrossprod(rep(1,n),rowSums(Z))
    Z <- Z*(t(L)%*%(X/(L%*%Z + 1e-9)))/tcrossprod(colSums(L),rep(1,l))



    # Renormalize so rows of Z have constant energy
    norms <- sqrt(rowSums(Z^2))
#    Z <- Z/tcrossprod(norms,rep(1,l))
#    L <- L*tcrossprod(rep(1,n),norms)

    Z <- Z/norms
    L <- L*norms

    # Calculate objective
    newobj <- sum(X*log((X+1e-10)/(L%*%Z)) - X + L%*%Z)
    objhistory <- newobj


  }

        rownames(L) <- rowna
        colnames(L) <- colnames(L, do.NULL = FALSE, prefix = "bicluster")
        clnames <- colnames(L)

        rownames(Z) <- clnames
        colnames(Z) <- colna

        M <- diag(p)
        rownames(M) <- clnames
        colnames(M) <- clnames

        LZ <- L%*%Z

        rownames(LZ) <- rowna
        colnames(LZ) <- colna

        U <- X-LZ

        rownames(U) <- rowna
        colnames(U) <- colna

  return(new('Factorization', parameters=list("nmfdiv",cyc,NULL,NULL,NULL,p,NULL,NULL,NULL,NULL,NULL,NULL,NULL),n=n,p1=p,p2=p,l=l,center=cent,scaleData=scaleData,X=X,L=L,Z=Z,M=M,LZ=LZ,U=U,avini=avini,xavini=xavini,ini=ini,Psi=as.vector(1),lapla=as.matrix(1)))
}


nmfeu <- function(X,p=5,cyc=100) {

        if (missing(X)) {
            stop("Data matrix X missing. Stopped.")
        }

        if (!is.matrix(X)) {
            X <- as.matrix(X)
        }

        if (!is.numeric(X)) {
            stop("Data matrix X must be numeric. Stopped.")
        }


# Check that we have non-negative data
if (min(X)<0) {
    stop("Negative values in matrix X. Stopped.")
}

# Dimensions
n <- nrow(X)
l <- ncol(X)

        rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "gene")
        colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
        rowna <- rownames(X)
        colna <- colnames(X)


       message("Running NMFEU on a ",n,"x",l," matrix with")
       message("   Number of biclusters ---------------- p: ",p)
       message("   Number of maximal iterations ------ cyc: ",cyc)




# Globally rescale data to avoid potential overflow/underflow
X <- X/max(X)


cent <- as.vector(rep(0.0,n))
scaleData <- as.vector(rep(max(X),n))

ini <- matrix(1,l,(p+1))
avini <- as.vector(rep(1.0,(p+1)))
xavini <- as.vector(rep(1.0,(l+1)))


# Create initial matrices
L <- matrix (abs(rnorm(n*p)),nrow=n,ncol=p)
Z <- matrix (abs(rnorm(l*p)),nrow=p,ncol=l)
#Z <- Z/(sqrt(rowsum(Z*Z,1:n))*as.vector(rep(1,l)))

# Calculate initial objective
objhistory <- 0.5*sum(sum((X-L%*%Z)^2))


# Start iteration
for (i in 1:cyc)
  {

    # Show progress
    print(i)
    print(objhistory)

    # Save old values
    Lold <- L
    Zold <- Z


    # Compute new L and Z (Lee and Seung; NIPS*2000)
    Z <- Z*(t(L)%*%X)/(t(L)%*%L%*%Z + 1e-9)
    L <- L*(X%*%t(Z))/(L%*%Z%*%t(Z) + 1e-9)

    # Renormalize so rows of Z have constant energy
    norms <- sqrt(rowSums(Z^2))
#    Z <- Z/tcrossprod(norms,rep(1,l))
#    L <- L*tcrossprod(rep(1,n),norms)

    Z <- Z/norms
    L <- L*norms



    # Calculate objective
    newobj <- sum(X*log((X+1e-10)/(L%*%Z)) - X + L%*%Z)
    objhistory <- newobj


  }

        rownames(L) <- rowna
        colnames(L) <- colnames(L, do.NULL = FALSE, prefix = "bicluster")
        clnames <- colnames(L)

        rownames(Z) <- clnames
        colnames(Z) <- colna

        M <- diag(p)
        rownames(M) <- clnames
        colnames(M) <- clnames

        LZ <- L%*%Z

        rownames(LZ) <- rowna
        colnames(LZ) <- colna

        U <- X-LZ

        rownames(U) <- rowna
        colnames(U) <- colna

        return(new('Factorization', parameters=list("nmfeu",cyc,NULL,NULL,NULL,p,NULL,NULL,NULL,NULL,NULL,NULL,NULL),n=n,p1=p,p2=p,l=l,center=cent,scaleData=scaleData,X=X,L=L,Z=Z,M=M,LZ=LZ,U=U,avini=avini,xavini=xavini,ini=ini,Psi=as.vector(1),lapla=as.matrix(1)))
}


projFuncPos <- function(s, k1, k2) {


  N <- length(s)


  ones <- as.vector(rep(1,N))



# Start by projecting the point to the sum constraint hyperplane
  v <- s + (k1-sum(s))/N

# Initialize zerocoeff (initially, no elements are assumed zero)

  # zerocoeff <- which(v<=0)
  # v[zerocoeff] <- 0

  zerocoeff <- vector(mode = "integer", length = 0)



  j <- 0
  ende <- 0
  while (ende<1) {

    # This does the proposed projection operator
    midpoint <- ones*k1/(N-length(zerocoeff))
    midpoint[zerocoeff] <- 0
    w <- v-midpoint
    a <- sum(w^2)
    b <- 2*crossprod(w,v)
    c <- sum(v^2)-k2
    t <- b^2-4*a*c
    if (t<0) {t <- 0}
    if (a<1e-10) {a <- 1e-10}

    alphap <- (-b+sqrt(t))/(2*a)
    v <- alphap*w + v

    if ((length(which(v<0))==0)&&(j>1)) {
      ende <- 2
    }
    else {

      j <- j+1

    # Set negs to zero, subtract appropriate amount from rest
      zerocoeff <- which(v<=0)
      v[zerocoeff] <- 0
      tempsum <- sum(v)
      tta <- N - length(zerocoeff)
      if (tta>0)  {
       v <- v + (k1-tempsum)/tta
       }
      v[zerocoeff] <- 0
    }
  }


  return(v)

}




# nmfsc - nonnegative matrix factorization with sparseness constraints
#
# Written by Sepp Hochreiter according to
# Patrik O. Hoyer
# 'Non-negative matrix factorization with sparseness constraints'
# Journal of Machine Learning Research  5:1457-1469, 2004.
#
#
# SYNTAX:
# res <- nmfsc( X, p, sL, sZ,cyc=100)
# res$L
# res$Z
#
# INPUTS:
# X          - data matrix
# p       - number of components (inner dimension of factorization)
# sL         - sparseness of L, in [0,1]
# sZ         - sparseness of Z, in [0,1]
# cyc    - maximal number of iterations
#
# Note: Sparseness is measured on the scale [0,1] where 0 means
# completely distributed and 1 means ultimate sparseness.
#
# NOTE: There is NO CONVERGENCE CRITERION.
#

nmfsc <- function(X,p=5,cyc=100,sL=0.6,sZ=0.6) {

        if (missing(X)) {
            stop("Data matrix X missing. Stopped.")
        }

        if (!is.matrix(X)) {
            X <- as.matrix(X)
        }

        if (!is.numeric(X)) {
            stop("Data matrix X must be numeric. Stopped.")
        }


# Check that we have non-negative data
if (min(X)<0) {
    stop("Negative values in matrix X. Stopped.")
}


# Data dimensions
n <- nrow(X)
l <- ncol(X)

        rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "gene")
        colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
        rowna <- rownames(X)
        colna <- colnames(X)

       message("Running NMFSC on a ",n,"x",l," matrix with")
       message("   Number of biclusters ---------------- p: ",p)
       message("   Number of maximal iterations ------ cyc: ",cyc)
       message("   Loading sparseness projection ------ sL: ",sL)
       message("   Factor sparseness projection ------- sZ: ",sZ)


cent <- as.vector(rep(0.0,n))
scaleData <- as.vector(rep(1.0,n))

# Create initial matrices
L <- matrix (rnorm(n*p),nrow=n,ncol=p)
Z <- matrix (rnorm(l*p),nrow=p,ncol=l)


ini <- matrix(1,l,(p+1))
avini <- as.vector(rep(1.0,(p+1)))
xavini <- as.vector(rep(1.0,(l+1)))


# Make initial matrices have correct sparseness
L1a <- sqrt(n)-(sqrt(n)-1)*sL
L1s <- sqrt(l)-(sqrt(l)-1)*sZ

for (i in 1:p)
  {
    L[,i] <- projFuncPos(s=as.vector(L[,i]),k1=L1a,k2=1.0)
    Z[i,] <- projFuncPos(s=as.vector(Z[i,]),k1=L1s,k2=1.0)
  }


# Calculate initial objective
objhistory <- sum((X-L%*%Z)^2)



# Initial stepsizes
stepsizeL <- 1.0
stepsizeZ <- 1.0



# Start iteration
for (i in 1:cyc)
  {


    # Save old values
    Lold <- L
    Zold <- Z

    # ----- Update Z ---------------------------------------


	# Gradient for Z
	dZ <- t(L)%*%(L%*%Z-X)

        begobj <- objhistory

        newobj <- 2.0*begobj




     # Make sure we decrease the objective!
	while (newobj>begobj) {


         # Take step in direction of negative gradient, and project
          Znew <- Z - stepsizeZ*dZ
	# Make initial matrices have correct sparseness
          for (j in 1:p)
            {
              Znew[j,] <- projFuncPos(s=as.vector(Znew[j,]),k1=L1s,k2=1.0)
            }


         # Calculate new objective
	    newobj <- sum((X-L%*%Znew)^2)

	    # If the objective decreased, we can continue...
            stepsizeZ <- stepsizeZ/2.0
            if (stepsizeZ<1e-10) {

                rownames(L) <- rowna
                colnames(L) <- colnames(L, do.NULL = FALSE, prefix = "bicluster")
                clnames <- colnames(L)

                rownames(Z) <- clnames
                colnames(Z) <- colna

                M <- diag(p)
                rownames(M) <- clnames
                colnames(M) <- clnames

                LZ <- L%*%Z

                rownames(LZ) <- rowna
                colnames(LZ) <- colna

                U <- X-LZ

                rownames(U) <- rowna
                colnames(U) <- colna

         return(new('Factorization', parameters=list("nmfsc",cyc,NULL,NULL,NULL,p,sL,sZ,NULL,NULL,NULL,NULL,NULL),n=n,p1=p,p2=p,l=l,center=cent,scaleData=scaleData,X=X,L=L,Z=Z,M=M,LZ=LZ,U=U,avini=avini,xavini=xavini,ini=ini,Psi=as.vector(1),lapla=as.matrix(1)))
            }
        }

	# Slightly increase the stepsize
	stepsizeZ <- stepsizeZ*1.2
	Z <- Znew



    # ----- Update L ---------------------------------------


	# Gradient for L
	dL <- (L%*%Z-X)%*%t(Z)
	begobj <- sum((X-L%*%Z)^2)

        newobj <- 2*begobj

	# Make sure we decrease the objective!
	while (newobj>begobj) {



	    # Take step in direction of negative gradient, and project
	    Lnew <- L - stepsizeL*dL
	    norms <- sqrt(sum(Lnew^2))
            for (j in 1:p)
              {
                Lnew[,j] <- projFuncPos(s=as.vector(Lnew[,j]),k1=L1a,k2=1.0)
              }

	    # Calculate new objective
	    newobj <- sum((X-Lnew%*%Z)^2)

	    stepsizeL <- stepsizeL/2

            if (stepsizeZ<1e-10) {
                rownames(L) <- rowna
                colnames(L) <- colnames(L, do.NULL = FALSE, prefix = "bicluster")
                clnames <- colnames(L)

                rownames(Z) <- clnames
                colnames(Z) <- colna

                M <- diag(p)
                rownames(M) <- clnames
                colnames(M) <- clnames

                LZ <- L%*%Z

                rownames(LZ) <- rowna
                colnames(LZ) <- colna

                U <- X-LZ

                rownames(U) <- rowna
                colnames(U) <- colna

              return(new('Factorization', parameters=list("nmfsc",cyc,NULL,NULL,NULL,p,sL,sZ,NULL,NULL,NULL,NULL,NULL),n=n,p1=p,p2=p,l=l,center=cent,scaleData=scaleData,X=X,L=L,Z=Z,M=M,LZ=LZ,U=U,avini=avini,xavini=xavini,ini=ini,Psi=as.vector(1),lapla=as.matrix(1)))
            }


          }

	# Slightly increase the stepsize
	stepsizeL <- stepsizeL*1.2
	L <- Lnew


    # Calculate objective
    newobj <- 0.5*sum(sum((X-L%*%Z)^2))
    objhistory <- newobj



  }
}

# modified eqscplot from the package MASS
#
# Plots with Geometrically Equal Scales
#
# From Version: 7.3-3, Date: 2009-10-15
# Maintainer: Brian Ripley <ripley@stats.ox.ac.uk>
#
#References:
#
#     Venables, W. N. and Ripley, B. D. (2002) _Modern Applied
#     Statistics with S._ Fourth edition.  Springer.

plotEqScale <- function (x, y, ratio = 1, tol = 0.04, uin, ...)
{
    dots <- list(...)
    nmdots <- names(dots)
    Call <- match.call()
    Call$ratio <- Call$tol <- Call$uin <- NULL
    if (is.matrix(x)) {
        y <- x[, 2]
        x <- x[, 1]
        if (!is.null(dn <- colnames(x))) {
            xlab0 <- dn[1L]
            ylab0 <- dn[2L]
        }
        else {
            xlab0 <- ""
            ylab0 <- ""
        }
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
        xlab0 <- "x"
        ylab0 <- "y"
    }
    else {
        xlab0 <- deparse(substitute(x))
        ylab0 <- deparse(substitute(y))
    }
    Call$x <- x
    Call$y <- y
    Call$xlab <- if ("xlab" %in% nmdots)
        dots$xlab
    else xlab0
    Call$ylab <- if ("ylab" %in% nmdots)
        dots$ylab
    else ylab0
    xlim <- if ("xlim" %in% nmdots)
        dots$xlim
    else range(x[is.finite(x)])
    ylim <- if ("ylim" %in% nmdots)
        dots$ylim
    else range(y[is.finite(y)])
    midx <- 0.5 * (xlim[2L] + xlim[1L])
    xlim <- midx + (1 + tol) * 0.5 * c(-1, 1) * (xlim[2L] - xlim[1L])
    midy <- 0.5 * (ylim[2L] + ylim[1L])
    ylim <- midy + (1 + tol) * 0.5 * c(-1, 1) * (ylim[2L] - ylim[1L])
    oldpin <- par("pin")
    xuin <- oxuin <- oldpin[1L]/abs(diff(xlim))
    yuin <- oyuin <- oldpin[2L]/abs(diff(ylim))
    if (missing(uin)) {
        if (yuin > xuin * ratio)
            yuin <- xuin * ratio
        else xuin <- yuin/ratio
    }
    else {
        if (length(uin) == 1L)
            uin <- uin * c(1, ratio)
        if (any(c(xuin, yuin) < uin))
            stop("'uin' is too large to fit plot in")
        xuin <- uin[1L]
        yuin <- uin[2L]
    }
    xlim <- midx + oxuin/xuin * c(-1, 1) * diff(xlim) * 0.5
    ylim <- midy + oyuin/yuin * c(-1, 1) * diff(ylim) * 0.5
    Call$xlim <- xlim
    Call$ylim <- ylim
    Call$xaxs <- Call$yaxs <- "i"
    Call[[1L]] <- as.name("plot")
    eval.parent(Call)
}






spfabia <- function(X,p=5,alpha=0.1,cyc=500,spl=0,spz=0.5,non_negative=0,random=1.0,write_file=1,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0,samples=0,initL=0,iter=1,quant=0.001,lowerB=0.0,upperB=1000.0){
	## X - name of the data file
	## cyc - maximum number of cycles
        ## alpha - sparseness
        ## p - factors

        if (missing(X)) {
            stop("Data file name missing. Stopped.")
        }


       message("Running sparse FABIA on a sparse matrix in file >",X,"<:")
       message("   Number of biclusters ---------------- p: ",p)
       message("   Sparseness factor --------------- alpha: ",alpha)
       message("   Number of iterations -------------- cyc: ",cyc)
       message("   Loading prior parameter ----------- spl: ",spl)
       message("   Factor prior parameter ------------ spz: ",spz)
        if (random>0) {
       message("   Initialization loadings--------- random: ",random," = positive")
       } else {
       message("   Initialization loadings--------- random: ",random," = all values")
      }
        if (non_negative>0) {
       message("   Nonnegative Loadings and Factors ------: ",non_negative," = Yes")
       } else {
       message("   Nonnegative Loadings and Factors ------: ",non_negative," = No")
      }
       if (write_file>0) {
       message("   Write results files -------------------: ",write_file, " = Yes")
       } else {
       message("   Write results files -------------------: ",write_file, " = No")
      }
        if (norm>0) {
       message("   Scaling to variance one: --------- norm: ",norm," = Yes")
           } else {
       message("   Scaling -------------------------- norm: ",norm," = No")
       }
       if (scale>0) {
       message("   Scaling loadings per iteration -- scale: ",scale," = to ", scale)
       } else {
       message("   Scaling loadings per iteration -- scale: ", scale," = No")
       }

       message("   Constraint variational parameter -- lap: ", lap)
       if ((nL>0)&&(nL<p)) {
           message("   Max. number of biclusters per row -- nL: ", nL)
           if (bL>0) {
           message("         starting at ------------------ bL: ", bL)
           } else {
           message("         starting at ------------------ bL: ", bL, " = from start")
           }
       } else {
           message("   Max. number of biclusters per row -- nL: ", nL, " = no limit")
       }
       if (lL>0) {
           message("   Max. number of row elements / biclu. lL: ", lL)
           if (bL>0) {
           message("         starting at ------------------ bL: ", bL)
           } else {
           message("         starting at ------------------ bL: ", bL, " = from start")
           }
       } else {
           message("   Max. number of row elements / biclu. lL: ", lL, " = no limit")
       }
       if (samples[1]==0) {
           message("   Number of samples ------------- samples: ", samples, " = all samples")
       } else {
           message("   Number of samples  ------------ samples: ", length(samples))
       }
       if (initL[1]==0) {
           message("   Random initialization of L ------ initL: ", initL)
       } else {
           message("   Biclusters init. by samples ----- initL: ", length(initL))
       }
           message("   Iterations ------------------------iter: ",iter)
        if (iter>1) {
           message("   Quantile of largest L removed --- quant: ",quant)
        }
           message("   Lower column sum bound --------- lowerB: ",lowerB)
           message("   Upper column sum bound --------- upperB: ",upperB)






        eps <- as.double(1e-3)
        eps1 <- as.double(1e-10)
        init_lapla <-  as.double(1.0)
        init_psi <-  as.double(0.2)

	samples <- as.integer(sort.int(as.integer(unique(samples))))
	initL <- as.integer(initL)
	iter <- as.integer(iter)
	norm <- as.integer(norm)
	cyc <- as.integer(cyc)
	nL <- as.integer(nL)
	lL <- as.integer(lL)
	bL <- as.integer(bL)
	quant <- as.double(quant)
	lowerB <- as.double(lowerB)
	upperB <- as.double(upperB)
	alpha <- as.double(alpha)
	non_negative <- as.integer(non_negative)
        write_file<- as.integer(write_file)
	p <- as.integer(p)
	spz <- as.double(spz)
	scale <- as.double(scale)
	lap <- as.double(lap)



	res <- .Call("spfabic",X,p,alpha,cyc,spl,spz,non_negative,random,write_file,init_psi,init_lapla,norm,scale,lap,nL,lL,bL,eps,eps1,samples,initL,iter,quant,lowerB,upperB,PACKAGE="fabia")

        if (is.null(res))
        {
            return(new('Factorization', parameters=list(),n=1,p1=1,p2=1,l=1,center=as.vector(1),scaleData=as.vector(1),X=as.matrix(1),L=noL,Z=nZ,M=as.matrix(1),LZ=as.matrix(1),U=as.matrix(1),avini=as.vector(1),xavini=as.vector(1),ini=as.matrix(1),Psi=as.vector(1),lapla=as.matrix(1)))
        }

        l=ncol(res$E_SX_n)
        n=nrow(res$L)

        pi=p*iter

        eps <- as.double(1e-3)
	nvect <- as.vector(rep(1,n))
        epsn<-eps*nvect

        iin <-  1.0/l

        # INI call for biclusters


        vz <- iin*apply(res$E_SX_n,1,function(x) sum(x^2))

        vz <- sqrt(vz+1e-10)

        ivz <- 1/vz

        if(length(ivz)==1) {
            nZ <- ivz*res$E_SX_n

            noL <- vz*res$L
        }
        else {
            nZ <- ivz*res$E_SX_n

            noL <- t(vz*t(res$L))
        }


        ini <- matrix(0,l,(pi+1))
        avini <- as.vector(rep(0.0,(pi+1)))
        xavini <- as.vector(rep(0.0,(l+1)))

        idp <- diag(pi)
        ppL <- matrix(0,pi,pi)

        for (ite in 0:(iter-1))
        {
            ppL[((ite*p)+1):((ite+1)*p),((ite*p)+1):((ite+1)*p)] <- crossprod(noL[,((ite*p)+1):((ite+1)*p)],(1/(res$Psi[((ite*n)+1):((ite+1)*n)]+epsn))*noL[,((ite*p)+1):((ite+1)*p)])
        }

        for (j in 1:l){
            mat <- idp + ppL/res$lapla[j,]
            ini[j,1:pi] <- log(diag(mat))
            s <- log(det(mat))
            ini[j,pi+1] <- s
            xavini[j] <- s
        }
        for (i in 1:pi){
            avini[i] <- sum(ini[,i])
        }

        ss <- sum(ini[,pi+1])
        xavini[l+1] <- ss
        avini[pi+1] <- ss



        if ((avini[pi+1]>1e-8)&&(pi>1)) {

            soo <- sort(avini[1:pi], decreasing = TRUE,index.return=TRUE)

            avini[1:pi] <- avini[soo$ix]
            noL <- noL[,soo$ix]
            nZ <- nZ[soo$ix,]

        }



        return(new('Factorization', parameters=list("spfabia",cyc,alpha,spl,spz,p,NULL,NULL,random,scale,norm,NULL,lap,nL,lL,bL,non_negative,write_file,init_lapla,init_psi,samples,initL,iter,quant,lowerB,upperB),n=n,p1=pi,p2=pi,l=l,center=as.vector(1),scaleData=as.vector(1),X=as.matrix(1),L=noL,Z=nZ,M=as.matrix(1),LZ=as.matrix(1),U=as.matrix(1),avini=avini,xavini=xavini,ini=ini,Psi=res$Psi,lapla=res$lapla))


}

readSamplesSpfabia <- function(X,samples=0,lowerB=0.0,upperB=1000.0){

        if (missing(X)) {
            stop("Data file name missing. Stopped.")
        }


	samples <- as.integer(sort.int(as.integer(unique(samples))))
	lowerB <- as.double(lowerB)
	upperB <- as.double(upperB)

	res <- .Call("readSamplesSpfabic",X,samples,lowerB,upperB,PACKAGE="fabia")

        if (is.null(res))
        {
            return(X=as.matrix(1))

        }

        return(res$X)


}


readSpfabiaResult <- function(X){
	## X - name of the data files

        if (missing(X)) {
            stop("Data file name missing. Stopped.")
        }

	res <- .Call("readSpfabicResult",X,PACKAGE="fabia")

        l=ncol(res$E_SX_n)
        n=nrow(res$L)

        p=ncol(res$L)

        iin <-  1.0/l

        # INI call for biclusters


        eps <- as.double(1e-3)
	nvect <- as.vector(rep(1,n))
        epsn<-eps*nvect

        vz <- iin*apply(res$E_SX_n,1,function(x) sum(x^2))

        vz <- sqrt(vz+1e-10)

        ivz <- 1/vz

        if(length(ivz)==1) {
            nZ <- ivz*res$E_SX_n

            noL <- vz*res$L
        }
        else {
            nZ <- ivz*res$E_SX_n

            noL <- t(vz*t(res$L))
        }


        ini <- matrix(0,l,(p+1))
        avini <- as.vector(rep(0.0,(p+1)))
        xavini <- as.vector(rep(0.0,(l+1)))

        idp <- diag(p)
        ppL <- crossprod(noL,(1/(res$Psi+epsn))*noL)


        for (j in 1:l){
            mat <- idp + ppL/res$lapla[j,]
            ini[j,1:p] <- log(diag(mat))
            s <- log(det(mat))
            ini[j,p+1] <- s
            xavini[j] <- s
        }
        for (i in 1:p){
            avini[i] <- sum(ini[,i])
        }

        ss <- sum(ini[,p+1])
        xavini[l+1] <- ss
        avini[p+1] <- ss



        if ((avini[p+1]>1e-8)&&(p>1)) {

            soo <- sort(avini[1:p], decreasing = TRUE,index.return=TRUE)

            avini[1:p] <- avini[soo$ix]
            noL <- noL[,soo$ix]
            nZ <- nZ[soo$ix,]

        }




        return(new('Factorization', parameters=list(),n=n,p1=p,p2=p,l=l,center=as.vector(1),scaleData=as.vector(1),X=as.matrix(1),L=noL,Z=nZ,M=as.matrix(1),LZ=as.matrix(1),U=as.matrix(1),avini=avini,xavini=xavini,ini=ini,Psi=res$Psi,lapla=res$lapla))


}

