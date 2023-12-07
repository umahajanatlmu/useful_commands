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
      library(i, character.only = TRUE)
    } else if (i %in% BiocManager::available(i)) {
      print(paste("installing:",i))
      BiocManager::install(i)
      library(i, character.only = TRUE)
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
