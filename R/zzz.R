# Init file for package fabia
.onLoad <- function(lib, pkg)
{
    library.dynam("fabia", pkg, lib)

    if ((.Platform$OS.type == "windows") && (.Platform$GUI ==
        "Rgui") && interactive()) {
        vigFile = system.file("Meta", "vignette.rds", package = "fabia")
        if (!file.exists(vigFile)) {
            warning(sprintf("fabia vignette is missing, nothing is added to the menu bar"))
        }
        else {
            vigMtrx = readRDS(vigFile)
            vigs = file.path(chartr("\\", "/", .find.package("fabia")), "doc", vigMtrx[,
                "PDF"])
            names(vigs) = vigMtrx[, "Title"]
            if (!"Vignettes" %in% winMenuNames())
                winMenuAdd("Vignettes")
            pkgMenu = paste("Vignettes", "fabia", sep = "/")
            winMenuAdd(pkgMenu)
            for (i in seq(along = vigs)) winMenuAddItem(pkgMenu,
                names(vigs)[i], paste("shell.exec(\"", vigs[i],
                  "\")", sep = ""))
        }
    }

  packageStartupMessage("+----------------------------+                                          \n",
      "|............................|                                          \n",
      "|............................|                                          \n",
      "|..............########......|  #######    #    ######    ###      #    \n",
      "|..............########......|  #         # #   #     #    #      # #   \n",
      "|.....####.....########......|  #        #   #  #     #    #     #   #  \n",
      "|.....####.....########......|  #####   #     # ######     #    #     # \n",
      "|.....####...................|  #       ####### #     #    #    ####### \n",
      "|.....####...........###.....|  #       #     # #     #    #    #     # \n",
      "|....................###.....|  #       #     # ######    ###   #     # \n",
      "|....................###.....|                                          \n",
      "|............................|                                          \n",
      "+----------------------------+                                          \n")

    version <- packageDescription("fabia",fields="Version")
    packageStartupMessage( "Citation: S. Hochreiter et al.,","\n",
      "FABIA: Factor Analysis for Bicluster Acquisition,","\n",
      "Bioinformatics 26(12):1520-1527, 2010.","\n",
      "BibTex: enter 'toBibtex(citation(\"fabia\"))'","\n\n",
      "Homepage: http://www.bioinf.jku.at/software/fabia/fabia.html","\n\n",
      "FABIA Package Version ", version, "\n")
}

.onUnload <- function(libpath)
{
    library.dynam.unload("fabia", libpath)
}


# Runs fabia demos interactively.
fabiaDemo <- function() {
  demoCode <- demo(package="fabia")$results[,3]
  if(length(demoCode) == 0) {
    return('No demos found in package fabia')
  }
  while(TRUE) {
    for(i in seq(1,length(demoCode))) {
      cat(i, ': ', demoCode[i], '\n')
    }
    suppressWarnings(demoNum <- as.integer(readline(prompt="Enter demo number to run [0 to view descriptions, RETURN to exit]: ")))

    if(is.na(demoNum)) {
      cat("Didn't get a number, exiting.\n")
      break
    }
    else if(demoNum == 0) {
      print(demo(package="fabia"))
    }
    else {
      if(demoNum < 0 || demoNum > length(demoCode)) {
        cat('Number out of range, try again...\n')
        next
      }
      demo(demoCode[demoNum], character.only=TRUE)
    }
  }
}

