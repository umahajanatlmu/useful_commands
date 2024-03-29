#' A ImputeTransformScale
#'
#' This function median Imputation, transformation and scaling by different
#' @param Data data-frame with column names
#' @param Impute whether to perform median imputation, default is TRUE
#' @param Transform whether to perform log 10 transformation, default is TRUE
#' @param Scaling whether to perform scaling, default is TRUE
#' @param ScaleType Scaling method, either of "Centering", "Auto", "Range","Pareto", "Vast", "Level".
#' @keywords Imputations, transformation, Scaling
#' @export
#' @examples
#' ImputeTransformScale (Data = Data, exclude = TRUE, dropList = NULL, Impute = TRUE, Transform = TRUE, Scaling = TRUE,
#' ScaleType = c("Centering", "Auto", "Range","Pareto", "Vast", "Level"))

ImputeTransformScale <- function(Data,
                                 Impute = FALSE,
                                 Transform = FALSE,
                                 Scaling = FALSE,
                                 ScaleType,
                                 drop.variables) {
  ## subset data
  nCol <- ncol(Data)
  ## empty data-frame
  scaleData <-
    data.frame(matrix(NA, nrow = nrow(Data), ncol = nCol))
  ## column names to empty dataframe
  colnames(scaleData) <- colnames(Data)
  
  for (i in 1:nCol) {
    ## select metabolite
    met <- colnames(Data)[[i]]
    if (!met %in% drop.variables) {
      if (class(Data[[i]]) %in% c("numeric", "integer")) {
        if (Impute == TRUE) {
          ## median of data
          impute <- median(Data[[i]], na.rm = TRUE)
          ## median imputation
          Data[is.na(Data[[i]]), i] <- impute
        }
        if (Transform == TRUE) {
          ## transformed data
          Data[[i]] <- log10(Data[[i]])
        }
        if (Scaling == TRUE) {
          ## mean of data
          mean <- mean(Data[[i]], na.rm = TRUE)
          ## centred data
          centreData <- Data[[i]] - mean
          ## sd data
          sd <- sd(centreData, na.rm = TRUE)
          ## min
          min <- min(Data[[i]])
          ## max
          max <- max(Data[[i]])
          ## centering
          if (ScaleType == "Centering") {
            scaleData[[met]] <- centreData
          }
          ## Auto-scaling
          else if (ScaleType == "Auto") {
            scaleData[[met]] <- centreData / sd
          }
          ## range scaling
          else if (ScaleType == "Range") {
            scaleData[[met]] <- (Data[[i]] - min) / (max - min)
          }
          ## pareto scaling
          else if (ScaleType == "Pareto") {
            scaleData[[met]] <- centreData / sqrt(sd)
          }
          ## vast scaling
          else if (ScaleType == "Vast") {
            scaleData[[met]] <- (centreData / mean) * (mean / sd)
          }
          ## level scaling
          else if (ScaleType == "Level") {
            scaleData[[met]] <- (centreData) / mean
          }
        }
      }
    } else
      scaleData[[met]] <- Data[[i]]
  }
  return(scaleData)
}