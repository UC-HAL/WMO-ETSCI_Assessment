#this function run trend analysis for single time series

trendAnalysis_func <- function (dataTS, acLagVal, alphaVal, valDigit) {
  
  #call and run the trend analysis function
  mkTest <- general_test_mk_prewhiten(dataTS, acLagVal, alphaVal)
  MK_pValue <- mkTest[[1]]
  MK_slope <- round(mkTest[[2]], valDigit)
  MK_intercept <- round(mkTest[[3]], valDigit)
  
  #return the outputs of interest
  return(list("pvalue"=MK_pValue, "slope"=MK_slope, "intercept"=MK_intercept))
  
} #end of function