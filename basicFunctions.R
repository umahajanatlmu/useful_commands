#' A ImputeTransformScale 
#'
#' This function is fro install and load packages
#' @param x list of packages to load
#' @keywords install, load cran and bioconductor packages 
#' @export
#' @examples
#' installScriptLibs(x)

installScriptLibs <- function(packages) {
  options(warn=-1) ## supress all warning
  for (i in packages) {
    if (i %in% .packages(all.available = TRUE)) {
      library(i, character.only = TRUE)
    } else if (i %in% names(utils::available.packages()[,1])) {
      print(paste("installing:",i))
      utils::install.packages(i, character.only = TRUE)
    } else if (i %in% BiocManager::available(i)) {
      print(paste("installing:",i))
      BiocManager::install(i)
    } else
      print(paste("package", i, "not found in CRAN or Bioconductor"))
  }
  options(warn=0) ## reset all warnings
}


#' print banner
#'
#' This function is to print header
#' @param txt text to print
#' @param txt symbols for brackets eg. "="
#' @keywords banner
#' @export
#' @examples
#' installScriptLibs(x)

banner <-function(txt, Char ="-") {
  nchar <- 64
  ## head tail
  headTail <- strrep(Char, nchar)
  hash <- paste0("##")
  charEnd <- strrep(Char,2)
  headTailSent <- paste0(hash, headTail)
  ## text
  textChar <- nchar(txt)
  if (textChar > nchar) {
    txt <- substring(txt, first = 1, last = 60)
  }
  space <- " "
  centering <- (nchar - 2) - textChar
  spacesBeforeAfter <- strrep(space, centering/2)
  
  textSent <- paste0(hash,
                     spacesBeforeAfter,
                     txt,
                     spacesBeforeAfter, 
                     charEnd)
  
  return(cat(paste0(headTailSent ,
                    "\n",textSent ,
                    "\n",headTailSent,
                    "\n")))
  
}

#' plot confusion matrix
#'
#' This function is for plotting confusion matriz
#' @param cm confusion matrix
#' @keywords install, load cran and bioconductor packages 
#' @export
#' @examples
#' plotCM(cm)

require(alluvial)

plotCM <- function(cm){
  cmdf <- as.data.frame(cm[["table"]])
  cmdf[["color"]] <- ifelse(cmdf[[1]] == cmdf[[2]], "#377EB8", "#E41A1C")
  
  alluvial::alluvial(cmdf[,1:2]
                     , freq = cmdf$Freq
                     , col = cmdf[["color"]]
                     , alpha = 0.5
                     , hide  = cmdf$Freq == 0
  )
}


# basic table theme  ------------------------------------------------------
datTable <- function(Dat) {
  if (!require(DT)) {
    print("installing DT-------------")
    install.packages("DT")
    library("DT")
  } else {
    library("DT")
  }
  
  table <- datatable(Dat, 
                     class = "cell-border stripe", 
                     rownames = TRUE, 
                     filter = "top", 
                     extensions = "Buttons", 
                     options = list(dom = "Bfrtip", 
                                    buttons = c("excel", "pdf", "print")))
  return(table)
  
}

# Function to plot Hotelling's T-squared ellipse
# Adapted from https://github.com/tyrannomark/bldR/blob/master/R/L2017.R
# GPL-3 license
gg_circle <- function(rx, ry, xc, yc, color = "black", fill = NA, ...) {
  x <- xc + rx * cos(seq(0, pi, length.out = 100))
  ymax <- yc + ry * sin(seq(0, pi, length.out = 100))
  ymin <- yc + ry * sin(seq(0, -pi, length.out = 100))
  annotate(
    "ribbon",
    x = x, ymin = ymin, ymax = ymax,
    color = color, fill = fill, ...
  )
}
