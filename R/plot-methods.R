#
#
# Author: SEPP HOCHREITER
###############################################################################



# The function plot.fabia is based on
# the function plot.mpm in the R package mpm
#
# Version: 1.0-16, Date: 2009-08-26
# Title: Multivariate Projection Methods
# Maintainer: Tobias Verbeke <tobias.verbeke@openanalytics.be>
# Author: Luc Wouters <wouters_luc@telenet.be>


setMethod("plot", signature(x="Factorization", y="missing"),
 function(x,
    Rm=NULL,
    Cm=NULL,
    dim = c(1, 2), # Factors to plot.
    zoom = rep(1, 2), # Zoom factor for factor scores and loadings
    col.group = NULL,
    colors = c("orange1", "red", rainbow(length(unique(col.group)), start=2/6, end=4/6)),
    col.areas = TRUE,
    col.symbols = c(1, rep(2, length(unique(col.group)))),
    sampleNames = TRUE,
    rot = rep(-1, length(dim)), # Mirror all axes
    labels = NULL, # character vector of labels (to allow labels to differ from row.names)
    label.tol = 0.1,
    lab.size = 0.725,
    col.size = 10,
    row.size = 10,
    do.smoothScatter = FALSE, # Plot individual points or density maps
    do.plot = TRUE, # This routine can also be used to calculate the
                    # coordinates, without plotting
                    # legend = FALSE
          ... )
{

  if (is.null(Rm)) {
    Rm=rep(1,nrow(x@L))
   }
  if (is.null(Cm)) {
    Cm=rep(1,ncol(x@Z))
  }
  if (is.null(col.group)) {
    col.group = rep(1, ncol(x@Z))
  }


  LI <- x@L
  ZI <- x@Z


  p <- nrow(ZI)

  if (p<2) {
        stop("At least two biclusters are needed for a 2dim plot. Stopped.")
    }

    if (dev.interactive()) {
      x11()
    }

  n <- nrow(LI)
  l <- ncol(ZI)

  X <- x@X
  if (ncol(X)==l) {
      colnames(ZI) <- colnames(X)
  }
  if (nrow(X)==n) {
      rownames(LI) <- rownames(X)

  }


  if (length(x@avini) < p) {
      avini <- rep(1,p)
  } else {
      avini <- x@avini[1:p]
  }

  if (do.smoothScatter && label.tol == 1){
    # if we require all points in the plot to be labelled, no points are left for the density map
    warning("All points selected for labelling, continuing without density map.\n")
    do.smoothScatter <- FALSE
  }
  if(do.smoothScatter && !do.plot){
    warning("Density map plotting requested but plotting not selected. Continuing with plot.\n")
    do.plot <- TRUE # If smoothScatter is set, then also do a plot
  }

  dLL <- 1/(apply(LI,2,function(x) max(abs(x)))+0.001*as.vector(rep(1,p)))
  #  dLL <- 1/sqrt((1/n)* apply(LI,2,function(x) sum(x^2))+0.001*as.vector(rep(1,p)))
  LI <- t(dLL*t(LI))

  dZZ <- 1/(apply(ZI,1,function(x) max(abs(x)))+0.001*as.vector(rep(1,p)))
  #  dZZ <- 1/sqrt( (1/l)*apply(ZI,1,function(x) sum(x^2))+0.001*as.vector(rep(1,p)))
  ZI <- dZZ*ZI



  # labeling of samples on the plot
  if (is.logical(sampleNames)){
    if (length(sampleNames) > 1){
      stop("'sampleNames' should either be a logical of length one or a character of length equal to the number of observations")
    } else {
      sampleLabels <- if (sampleNames) colnames(ZI) else rep("", l)
    }
  } else {
    if (length(sampleNames) != l){
      stop("If 'sampleNames' is a character vector, the length should be equal to the number of observations")
    } else {
      sampleLabels <- sampleNames
    }
  }


  if (is.data.frame(col.group))
    col.group <- unlist(col.group)
  if (length(col.group) != l)
    stop("Length of 'col.group' not equal to number of columns in data.")
  col.group <- as.numeric(as.factor(col.group))


  sud <- sum(avini)
  if (sud>1e-8) {

      soo <- sort(avini, decreasing = TRUE,index.return=TRUE)

      d <- avini[soo$ix]
      U <- LI[,soo$ix]
      V <- t(ZI[soo$ix,])
      contrib <- d/sud
  } else {
      U <- LI
      V <- t(ZI)
      contrib <- rep((1/p),p)
  }

  #
  # Scaling, scores and loadings
  #
  ### scaling: alpha and beta values for different scaling options
  # Calculate weighted factor scores (for rows)
  Z <- V[, dim]

  # Calculate weighted factor loadings (for columns)
  L <- U[,dim]

  # Rotation
  # Default: flips both the X and Y axes
  Z <- Z * matrix(rot, ncol = ncol(Z), nrow = nrow(Z), byrow = TRUE)
  L <- L * matrix(rot, ncol = ncol(L), nrow = nrow(L), byrow = TRUE)

  #
  # zooming
  zz <- Z * zoom[1]
  ll <- L * zoom[2]


  #
  # Calculate which rows are far from the origin (most specific)
  #

  iz <- rep(TRUE, nrow(Z))

  il <- rep(TRUE, nrow(L))

  if (length(dim) == 2){
      DL <- L[, 1]^2 + L[, 2]^2  # distance from center
  # Threshold for labels
      if (label.tol > 1) {# label.tol most distant rows are plotted as circles and labelled
          thres <- sort(DL, decreasing = TRUE)[min(length(DL), floor(label.tol))]
      }
      else {# label.tol percent most distant rows are plotted as circles and labelled
          thres <- if (label.tol == 0) Inf else quantile(DL, probs = 1-label.tol)}

  #
  # Positioning
  #
      iselt <- il & (DL >= thres)
  } else {


      DL <- L[, 1]^2 + L[, 2]^2 + L[, 3]^2  # distance from center
  # Threshold for labels
      if (label.tol > 1) {# label.tol most distant rows are plotted as circles and labelled
          thres <- sort(DL, decreasing = TRUE)[min(length(DL), floor(label.tol))]
      }
      else {# label.tol percent most distant rows are plotted as circles and labelled
          thres <- if (label.tol == 0) Inf else quantile(DL, probs = 1-label.tol)
      }

  #
  # Positioning
  #
      iselt <- il & (DL >= thres)

  }


  isel <- which(iselt)
  # Only draw a plot if the projection is 2D and the user requested a plot (default)
 #
  if (do.plot && (length(dim) == 2)){

   # Compute range of plot
    #

   xrange <- range(ll[, 1], zz[, 1], 0)
    yrange <- range(ll[, 2], zz[, 2], 0)
    xrange <- ifelse(rep(diff(xrange) == 0, 2),
        c(-1, 1), # Default range [-1,1] if calculated range is 0
        xrange + diff(xrange) * c(-0.1, 0.1)) # Else expand range by 10%
    yrange <- ifelse(rep(diff(yrange) == 0, 2),
        c(-1, 1),
        yrange + diff(yrange) * c(-0.1, 0.1))

    zz[,1] <- pmin(pmax(zz[,1], xrange[1]), xrange[2])
    zz[,2] <- pmin(pmax(zz[,2], yrange[1]), yrange[2])
    ll[,1] <- pmin(pmax(ll[,1], xrange[1]), xrange[2])
    ll[,2] <- pmin(pmax(ll[,2], yrange[1]), yrange[2])

   #
    # Set-up plot
    #
    opar <- par(pty = "m") # preserve configuration
    # Create a window with the maximal plotting region
    # Equal scale plot function from MASS library
    dotList <- list(...)
    if (is.null(dotList$sub)){
      dotList$sub <- paste("FABIA", sep="")
      if (is.null(dotList$cex.sub))
        dotList$cex.sub <- 0.85
    }
    dotList$cex.sub <- if (is.null(dotList$cex.sub)) 0.85 else dotList$cex.sub
    if (is.null(dotList$xlab))
      dotList$xlab <- paste("BC", dim[1], ": ", 100 * round(contrib[dim[1]], 2), "%,  ", floor(d[dim[1]]),sep = "")
    if (is.null(dotList$ylab))
      dotList$ylab <- paste("BC", dim[2], ": ", 100 * round(contrib[dim[2]], 2), "%,  ", floor(d[dim[2]]), sep = "")

    dotList$x <- xrange
    dotList$y <- yrange
    dotList$ratio <- 1
    dotList$tol <- 0
    dotList$type <- "n"
    dotList$axes <- FALSE
    dotList$cex.lab = 0.85


    do.call("plotEqScale", dotList)


    #
    # Scales
    usr <- par("usr") # Retrieve extremes of coordinates in the plotting region
    sx <- diff(usr[c(1,2)]) / (25.4 * par("pin")[1, drop = TRUE]) # scale in mm
    sy <- diff(usr[c(3,4)]) / (25.4 * par("pin")[2, drop = TRUE])

  #
    # Plot rows close to 0,0 as unlabelled dots or as density maps
    #
    # Select rows to be plotted
    idlt <- il & (DL < thres)

    idl <- which(idlt)

    if (!do.smoothScatter) # plot inner points as unlabelled dots
      points(ll[idl,1], ll[idl,2], col = colors[1], cex = 0.825, lwd = 2)
    else # plot as density maps
      smoothScatter(x = ll[idl,1], y = ll[idl,2], nbin = 256, nrpoints = 0,
          add = TRUE, colramp = colorRampPalette(c("white", "burlywood")))


    ### plot distant rows as circles with areas proportional to Rm
    sqs <- 0.5 * sx * pmax(0.02, row.size * sqrt((Rm) / (max(Rm))))
    yoffset <- sy * (2 + sqs / sx)






    if (length(isel) > 0){ # if there is at least 1 point to plot
      symbols(ll[isel, 1], ll[isel, 2], circle = sqs[isel],
        inches = FALSE, lwd = 3, add = TRUE, fg = colors[2])
      if (is.null(labels)) labels <- rownames(LI)
      text(ll[isel, 1], ll[isel, 2] - yoffset[isel], adj = c(0.5, 1),
        cex = lab.size, labels = labels[isel], # x$row.names[isel]
        col = colors[2])
    }

    ### plot columns with indication of column-grouping
    iGroup <- unique(col.group[iz]) # unique groups in columns to be plotted
    for (i in 1:length(iGroup)){
      ii <- iz & (col.group == iGroup[i]) # Select columns in group i selected for plotting
      if (col.areas)
      { # use squares with size (ie. area) proportional to x$Cm
        sqs <- 0.5 * sx * pmax(0.02, col.size * sqrt((Cm[ii]) / (max(Cm))))
        yoffset <- sy * (5 + sqs / (2 * sx))
        symbols(zz[ii, 1], zz[ii, 2],
            square = sqs, inches = FALSE, lwd = 3, add = TRUE, fg = colors[2+iGroup[i]])
        # if (sampleNames){
          text(zz[ii,1], zz[ii,2] + yoffset,
              adj=c(0.5, 1), cex=lab.size, labels=sampleLabels[ii], col=colors[2+iGroup[i]]) # x$col.names[ii]
        # }
      }
      else # Use different symbols, ignore size
      {
        yoffset <- sy * (5 + col.size / 5)
        points(zz[ii, 1], zz[ii, 2],
            pch = col.symbols[iGroup[i]], col = colors[2+iGroup[i]],
            cex = col.size / (25.4 * par("csi")), lwd = 3)
        text(zz[ii, 1], zz[ii, 2] + yoffset,
            adj = c(0.5, 1), cex = lab.size, labels = sampleLabels[ii], col = colors[2+iGroup[i]]) # x$col.names[ii]
      }
    }

    ### finish plot, put crozz on 0,0, box, and legend
    lines(c(-2.5,2.5) * sx, c(0,0), lwd=3)
    lines(c(0,0), c(-2.5,2.5) * sy, lwd=3)
    box()

    # legend
#    legendArg <- match.arg(legend)
#    if (legendArg != "none"){
#      legend(legendArg,
#        legend = levels(pData(ALL)$BT),
#        text.col = col.group, bty='n')
#    }

    par(opar) # Restore plotting configuration
}



                                        #
  # Return value: list of coordinates and indication of most distant (most specific) points
  #

  if (length(dim) == 1){
      rows <- cbind(ll, iselt)
      dimnames(rows) <- list(rownames(LI), c("X", "Select"))
    dimnames(zz) <- list(colnames(ZI), c("X"))
  } else if (length(dim) == 2){
      rows <- cbind(ll, iselt)
      dimnames(rows) <- list(rownames(LI), c("X", "Y", "Select"))
      dimnames(zz) <- list(colnames(ZI), c("X", "Y"))
  } else if (length(dim) == 3){
    dimnames(rows) <- list(rownames(LI), c(paste("Prf", dim, sep=""), "Select"))
    dimnames(zz) <- list(colnames(ZI), paste("Pcf", dim, sep=""))
  } else { # More than 3 dimensions requested
    dimnames(rows) <- list(rownames(LI), c(paste("Prf", dim, sep=""), "Select"))
    dimnames(zz) <- list(colnames(ZI), paste("Pcf", dim, sep=""))
  }
  r <- list(Rows = rows, Columns = zz)
})
