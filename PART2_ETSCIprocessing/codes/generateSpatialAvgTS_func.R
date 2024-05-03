#this function generates spatially averaged time series from gridded data


generateSpatialAvgTS_func <- function(inputDat, computDim, valRound) {
  
  #spatially averaged time series
  spatialAvgTS <- round(apply(inputDat, computDim, mean, na.rm=TRUE), valRound)
  spatialAvgTS[is.infinite(spatialAvgTS)] <- NA #replacing infinite values (Inf or -Inf)
  
  #return the output of interest
  return(spatialAvgTS)
  
} #end of function