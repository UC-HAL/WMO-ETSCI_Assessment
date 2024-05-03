#this function loads required libraries and sources user-defined functions called throughout the workflow main script
libLoad_funcSource_func <- function() {
  
  # Load necessary libraries (install them 1st if not yet done)
  library(ncdf4) #manipulate data in netcdf format
  library(ncdf4.helpers) #manage time dimension in dates format
  
  options("sp_evolution_status"=2) #adding the option to prevent sp to call most code in rgdal or rgeos retiring packages (as of october 2023); visit https://r-spatial.org/r/2023/05/15/evolution4.html
  library(sp) #required for 'raster' package; 
  
  library(raster) #manipulate raster objects
  library(sf) #work with geospatial data (CAUTION: IMPORTANT TO CALL IT AFER THE RASTER PACKAGE!!!)
  
  library(reshape2) #convert data frame to long format (melt func)
  library(ggplot2) #create graphs (https://stackoverflow.com/tags/ggplot2)
  
  library(RColorBrewer) # for color palettes
  library(viridisLite) #required for "viridis" package
  library(viridis) #Colorblind-Friendly Color Maps for R
  
  library(greekLetters) #for greek letters like Deltas [https://cran.r-project.org/web/packages/greekLetters/greekLetters.pdf]
  library(stringr) #operations with strings (https://stringr.tidyverse.org/)
  
  library(gridExtra) #for arranging plots
  library(grid) #for text_grob function use
  library(ggpubr) #to use ggarrange function to add a common Legend for combined ggplots
  
  
  #source all the functions (created for data processing) to be called: CAUTION: VERY IMPORTANT TO DO THAT!!!
  source(paste0(wkDir, slash, "ncfile_cropping_func.R")) #crop netCDF to the extent of the target region shapefile (spatial analysis), and mask area outside cropped raster (temporal analysis)
  source(paste0(wkDir, slash, "generateSpatialStats_func.R")) #calculate spatial statistics (mean, median, min and max) over a given time period, from gridded annual time series
  source(paste0(wkDir, slash, "generateSpatialAvgTS_func.R")) #spatially average annual time series over target region for a given time period, from gridded annual time series
  source(paste0(wkDir, slash, "generateSpatialDataRange_func.R")) #get region-wide min., max., mean, and median annual time series from original cropped (but not masked) gridded annual time series for the study full time period
  source(paste0(wkDir, slash, "getIndicesDataRange_forSpatialMaps_func.R")) #spatial maps colorbar limit setting: get for each index, data range common to all data sources, time periods and target regions.
  source(paste0(wkDir, slash, "generateSpatialMaps_func.R")) #create spatial maps of descriptive statistics (mean, median, min and max) calculated over a given time period, from gridded annual time series (original and changes)
  source(paste0(wkDir, slash, "getSpatialClimatologieMaps_func.R")) #get spatial maps from long-term gridded descriptive statistics
  source(paste0(wkDir, slash, "trendAnalysis_func.R")) #run trend analysis for single time series
  source(paste0(wkDir, slash, "aggTSsubplots_func.R")) #aggregate annual time series subplots into a single plot for each index and target region
  source(paste0(wkDir, slash, "generateTemporalStats_func.R")) #calculate long-term statistics (mean, sd, min and max) from temporal annual time series
  source(paste0(wkDir, slash, "getTempDeltasAggBoxplots_func.R")) #aggregate temporal change boxplots for different groups of indices
  
  #source functions developed for trend analysis by Andrew Tefs (UC-HAL: https://ucalgary.ca/labs/hydrological-analysis/team)
  	#note that these are within another directory
  source(paste0(trendFunc_dir, slash, "general_test_mk_prewhiten.R"))
  source(paste0(trendFunc_dir, slash, "general_ac_lag.R"))
  source(paste0(trendFunc_dir, slash, "general_mk_raw.R"))
  source(paste0(trendFunc_dir, slash, "general_sen_slope.R"))
  
} #end of function

