#this function generates spatially averaged time series from original cropped (but non-masked) gridded data

    #these time series are used to define spatial maps colorbar limits common to all data sources, target regions, and time periods evaluated.


generateSpatialDataRange_func <- function(inputDat, computDim, valRound) {
  
  #region-wide single annual min. time series
  spatialMinTS <- round(apply(inputDat, computDim, min, na.rm=TRUE), valRound)
  spatialMinTS[is.infinite(spatialMinTS)] <- NA #replacing infinite values (Inf or -Inf)
  
  #region-wide single annual max. time series
  spatialMaxTS <- round(apply(inputDat, computDim, max, na.rm=TRUE), valRound)
  spatialMaxTS[is.infinite(spatialMaxTS)] <- NA #replacing infinite values (Inf or -Inf)
  
  #region-wide single annual mean time series
  spatialMeanTS <- round(apply(inputDat, computDim, mean, na.rm=TRUE), valRound)
  spatialMeanTS[is.infinite(spatialMeanTS)] <- NA #replacing infinite values (Inf or -Inf)
  
  #region-wide single annual median time series
  spatialMedianTS <- round(apply(inputDat, computDim, median, na.rm=TRUE), valRound)
  spatialMedianTS[is.infinite(spatialMedianTS)] <- NA #replacing infinite values (Inf or -Inf)
  
  #return the output of interest
  return(list("minTS"=spatialMinTS, "maxTS"=spatialMaxTS, 
              "meanTS"=spatialMeanTS, "medianTS"=spatialMedianTS))
  
} #end of function