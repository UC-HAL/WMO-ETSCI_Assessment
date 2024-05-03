#this function generates long-term descriptive statistics from temporal annual time series
#(i.e., mean, sd, min, max)

generateTemporalStats_func <- function(inputDat, valRound) {
  
  #descriptive statistics to calculate: mean, median, min, max
  TemporalMean <-  round(mean(inputDat, na.rm = TRUE), valRound)
  TemporalSD <- round(sd(inputDat, na.rm = TRUE), valRound)
  TemporalMin <- round(min(inputDat, na.rm = TRUE), valRound)
  TemporalMax <- round(max(inputDat, na.rm = TRUE), valRound)
  
  #return the outputs of interest
  return(list("mean"=TemporalMean, "sd"=TemporalSD, 
              "min"=TemporalMin, "max"=TemporalMax))
  
} #end of function