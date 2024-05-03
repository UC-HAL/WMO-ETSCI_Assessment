general_sen_slope <- function ( timeSeries ) {
  
  # timeseries   vector    timeSeries of data (can contain NAs, but needs spacing preserved)
  
  n <- length(timeSeries) # find number of real data points
  pairSlope <- NA # initialize slopes with dummy entry
  
  for (i in 1:(n - 1)) { # FOR all entries up to second last
    for (ii in (i + 1):n) { # FOR all entries after cureent first
      pairSlope <- c(pairSlope,
                     (timeSeries[ii] - timeSeries[i]) / (ii - i)) # compute all pairwise slopes
    } # ENDFOR second entry to end
  } # ENDFOR all entries up to second last
  
  sensSlope <- median(pairSlope, 
                      na.rm = TRUE) # median of paired slopes from sequential points
  sensInter <- median(timeSeries - sensSlope * seq(1,
                                                   n),
                      na.rm = TRUE) # median intercept from slope
  
  return(list("slope" = sensSlope,
              "inter" = sensInter)) # RETURN Sen-Theil characteristics 
}