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

### ggplot publication theme
theme_publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = rel(1.2)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background =  element_rect(colour = NA),
            panel.border = element_rect(colour="black", linewidth= rel(1)),
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            #axis.line = element_line(colour="black", size = rel(0.5)),
            axis.ticks = element_line(colour="black", linewidth = rel(0.5)),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            #legend.direction = "horizontal",
            legend.key.size= unit(rel(0.2), "cm"),
            legend.margin = margin(rel(0.01),rel(0.01),rel(0.01),rel(0.01), "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(rel(10),rel(5),rel(5),rel(5)),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


