general_ac_lag <- function ( timeSeries,
                             acLag ) {
  
  # timeseries   vector    timeSeries of data (can contain NAs, but needs spacing preserved)
  # acLag        integer   number of timesteps used for lag
  
  n <- sum(!is.na(timeSeries)) # find number of real data points
  mu <- mean(timeSeries,
             na.rm = TRUE) # mean of timeseries
  
  rawTerm <- (timeSeries[(1 + acLag):n] - mu) # variance of un-lagged data
  lagTerm <- (timeSeries[1:(n - acLag)] - mu) # variance of lagged data
  covSum <- sum(rawTerm * lagTerm,
                na.rm = TRUE) # sum of covariance vectors of lagged, un-lagged
  autocovar <- sum((timeSeries - mu)^2,
               na.rm = TRUE) # autocovariance of un-lagged data 
  
  ack <- (covSum / (n - acLag)) / (autocovar / n) # autocorrelation coeff for given lag
  
  return(ack) # RETURN autocorrelation coefficient
}