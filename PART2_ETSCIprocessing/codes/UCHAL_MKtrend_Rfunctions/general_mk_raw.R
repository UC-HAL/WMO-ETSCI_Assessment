general_mk_raw <- function ( timeSeries ) {
  
  # timeseries   vector    timeSeries of data (can contain NAs, but needs spacing preserved)
  
  n <- sum(!is.na(timeSeries)) # find number of data points
  mkZ <- 0 # initialize z-score to zero
  
  mkS <- outer(timeSeries,
               timeSeries,
               "-") # sign-change term
  mkS <- sum(sign(mkS) * lower.tri(mkS),
             na.rm = TRUE) # sign-change summation
  mkSvar <- n * (n - 1) * (2 * n + 5) / 18 # variability term
  
  if (mkS > 0) { # IF increasing
    mkZ <- (mkS - 1) / sqrt(mkSvar) # compute Z distribution term
    
  } else if (mkS < 0) { # ELSEIF decreasing
    mkZ <- (mkS + 1) / sqrt(mkSvar) # compute Z distribution term
    
  } # ENDIF increasing
  
  pValue <- 2 * pnorm(-1 * abs(mkZ)) # convert to two-sided p-value
  
 return(pValue) # RETURN z-score p-value for normal distribution
}