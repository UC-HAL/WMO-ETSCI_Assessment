#this function generates spatial descriptive statistics from gridded data
    #(i.e., gridded mean, median, min, max)

generateSpatialStats_func <- function(inputDat, computDim, valRound) {
  
  #descriptive statistics to calculate: mean, median, min, max
  spatialMean <- round(apply(inputDat, computDim, mean, na.rm=TRUE), valRound)
  spatialMedian <- round(apply(inputDat, computDim, median, na.rm=TRUE), valRound)
  spatialMin <- round(apply(inputDat, computDim, min, na.rm=TRUE), valRound)
  spatialMax <- round(apply(inputDat, computDim, max, na.rm=TRUE), valRound)
  
  #replacing infinite values (Inf or -Inf)
  spatialMean[sapply(spatialMean, is.infinite)] <- NA
  spatialMedian[sapply(spatialMedian, is.infinite)] <- NA
  spatialMin[sapply(spatialMin, is.infinite)] <- NA
  spatialMax[sapply(spatialMax, is.infinite)] <- NA
  
  #return the outputs of interest
  return(list("mean"=spatialMean, "median"=spatialMedian, 
              "min"=spatialMin, "max"=spatialMax))
  
} #end of function