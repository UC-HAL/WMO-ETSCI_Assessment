general_test_mk_prewhiten <- function ( timeSeries, 
                                        acLag, 
                                        alpha ) {
  
  # Compute Mann-Kendall significance, Sen's slope, Sen's intercept of pre-whtiened timeseries 
  # Andrew Tefs, January 2021, University of Calgary Geography, andrew.tefs@ucalgary.ca
  
  # timeseries   vector    timeSeries of data (can contain NAs, but needs spacing preserved)
  # acLag        integer   number of timesteps used for lag
  # alpha        double    significance of Mann-Kendall test as fraction (i.e. 0.05 for 5%)
  
  # Method for correcting trend tests summarized in pg. 1821 - 1822:
  # Yue et al., 2002, Hydrological Processes, DOI: 10.1002/hyp.1095,
  # "The influence of autocorrelation on the ability to detect trend..."
  
  # !! STILL NEEDS SECTION ADDED TO ACCOUNT FOR TIES TO MATCH YUE ET AL. (2002) EXACTLY !!
  
  n <- length(timeSeries) # entries in timesries
  slpThr <- 0.01 * mean(timeSeries,
                        na.rm = TRUE) / n # slopes shallower than this will be excluded (0.01%)
  acThr <- abs(qnorm(alpha / 2)) / sqrt(n - acLag) # AC coefficient above which pre-whitening is necessary
  
  acRaw <-  general_ac_lag(timeSeries,
                         acLag) # compute AC coefficient for k steps of lag (specified in input)
  pValue <- general_mk_raw(timeSeries) # compute Mann-Kendall p-value of non-transformed data
  senData <- general_sen_slope(timeSeries) # compute the raw Sen Thiel slope and intercept
  senSlp <- senData$slope # Sen's slope
  senInt <- senData$inter # intercept relative to Sen's slope
  
  if (abs(senSlp) <= slpThr) { # IF Sen's slope is too small to detect significance (corrupts trend MK calculation)
    pValue <- NA # MK p-value is not detectable due to tiny slope, return no-data value
    
  } else { # ELSE Sen's slope is too small to detect significance (corrupts trend MK calculation)
    if (abs(acRaw) >= acThr) { # IF autocorrelation coefficient exceeds critical value
      timeSeriesDT <- timeSeries - senSlp * seq(1, 
                                                n) # de-trend data by removing linear equivalent of Thiel-Sen slope
      acDTr <- general_ac_lag(timeSeriesDT,
                              acLag) # compute autocorreatlion component of timeseries for lag k
      
      acComponent <- acDTr * timeSeriesDT[1:(n - acLag)] # k-timestep lagged autocorrelation component
      timeSeriesPW <- c(timeSeriesDT[1:acLag],
                        timeSeriesDT[(1 + acLag):n] - acComponent) # remove k-timestep lagged autocorrelation component
      
      timeSeriesBl <- timeSeriesPW + senSlp * seq(1,
                                                  n) # blend timeseries by re-applying trend to the pre-whitened data
      pValue <- general_mk_raw(timeSeriesBl) # re-compute Mann-Kendall test on blended timeseries
      
    } # ENDIF autocorrelation coefficient exceeds critical value
  } # ENDIF Sen's slope is too small to detect significance (corrupts trend MK calculation)
  
  return(list("pVal" = pValue, 
              "slope" = senSlp, 
              "interc" = senInt)) # RETURN p-value, Sen slope, and intercept of Sen slope
} 