# R version 4.3.1 (2023-06-16 ucrt) 
# Development period: Fall 2023 throughout Winter 2024 | by Alida Thiombiano (alida.thiombiano@ucalgary.ca)

# This script enable spatial and temporal analyses of the ET-SCI annual time series
# resulting in a climate data-driven assessment of hydroclimatic annual conditions for any target region globally.

###############################
  #NOTE1: To enable analysis at monthly and/or seasonal scales, 
          #the user may select the indices with monthly time series,
          #use the CDO software (https://code.mpimet.mpg.de/projects/cdo/embedded/index.html) to aggregate them seasonally,
          #and adapt the current script and functions accordingly.

  #NOTE2: SECTIONS 1, 2 and 3 contain the DEFAULT SETTINGS for the gridded indices (netCDF files) pre-processings and analyses in the subsequent sections.
          #Hence, lot of hard coding there, and MOST OF THE necessary modifications in this workflow should be brought in these top 3 sections ONLY!!!
          #However, the User's choices made in Sections1-3, may require some adjustment across the subsequent sections where hard coding are made based on these beforehand settings.

  #NOTE3: THIS WORKFLOW SHOULD BE USED WITH DATA SOURCES SHARING SAME TIME COVERAGE LIKE CMIP SIMULATIONS (https://www.wcrp-climate.org/wgcm-cmip; e.g., CMIP6).
              #THE ANALYSES FOCUS ON LONG-TERM CLIMATOLOGIES, TRENDS, AND FUTURE CHANGES ASSESSMENT.
          #for observational datasets analysis, a copy of this workflow can be made and adapted to fit the time periods of interest.
              #future changes quantification may then be replaced by anomalies and normals assessment!
###############################



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>SECTION1>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#set working directory and enable clean connection between paths and filenames using r"(path)"
wkDir <- r"(C:\Users\alida.thiombiano\OneDrive - University of Calgary\WMO-ETSCI_Assessment_R-Workflow\PART2_ETSCIprocessing\codes)"
slash <- "\\" #clean connection

#source then run the user-defined function which load libraries and source user-defined functions called within this script
source(paste0(wkDir, slash, "libLoad_funcSource_func.R"))
	#note that functions developed for trend analysis by Andrew Tefs (UC-HAL: https://ucalgary.ca/labs/hydrological-analysis/team) are within anoter directory
trendFunc_dir <- r"(C:\Users\alida.thiombiano\OneDrive - University of Calgary\WMO-ETSCI_Assessment_R-Workflow\PART2_ETSCIprocessing\codes\UCHAL_MKtrend_Rfunctions)"
libLoad_funcSource_func()



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>SECTION2>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#Set paths to directory

    #Path to ET-SCI netCDF
ETSCI_dir <- r"(C:\Users\alida.thiombiano\OneDrive - University of Calgary\WMO-ETSCI_Assessment_R-Workflow\PART1_ETSCIcalculation\ETSCIfromGriddedTS\GriddedETSCI)"

    #Path to analysis outputs
outpath1 <- r"(C:\Users\alida.thiombiano\OneDrive - University of Calgary\WMO-ETSCI_Assessment_R-Workflow\PART2_ETSCIprocessing\spatialAnalysis)"
outpath2 <- r"(C:\Users\alida.thiombiano\OneDrive - University of Calgary\WMO-ETSCI_Assessment_R-Workflow\PART2_ETSCIprocessing\temporalAnalysis)"

#create an empty list and store the time periods to be assessed
    #in the present case study, the ET-SCI were driven by climate models from CMIP6 with a temporal coverage from 1950 to 2100
tperiods <- vector("list", 6)
tperiods[[1]] <- c(1950:2100) #full temporal coverage ; 
tperiods[[2]] <- c(1950:2014) #historical period
tperiods[[3]] <- c(2015:2100) #future period
tperiods[[4]] <- c(1985:2014) #baseline period
tperiods[[5]] <- c(2041:2070) #near future
tperiods[[6]] <- c(2071:2100) #far future

tper_strings <- c("1950-2100", "1950-2014", "2015-2100", "1985-2014", "2041-2070", "2071-2100")
    #get their position for further processing
nfull <- 1; nhist <- 2; nfut <- 3; nref <- 4; nfut1 <- 5; nfut2 <- 6

#defining time periods for further statistical analyses
    #for trend analysis, time evolution plot, anomalies calculation (temporal scale)
fullPer <- tperiods[[nfull]]; fullHistPer <- tperiods[[nhist]]; fullFutPer <- tperiods[[nfut]]
        #derive names strings
fullPer_label <- tper_strings[nfull]; fullHist_label <- tper_strings[nhist]; fullFut_label <- tper_strings[nfut]
    #to quantify projected future changes (spatial and temporal scales)
refPer <- tperiods[[nref]]; nearFut <- tperiods[[nfut1]]; farFut <- tperiods[[nfut2]] 
        #derive the names strings
refLabel <- tper_strings[nref]; fut1Label <- tper_strings[nfut1]; fut2Label <- tper_strings[nfut2]
    #define variable for temporal assessment of future changes
nyears <- 30

###################################CAUTION!!!
#exploratory analysis of data show an extra year (i.e., 2101 or 152yrs, or 12x152=1824months)
    #this is straightforward addressed for indices' original annual time series
    #for indices with monthly time series like SPEI, below variable will be used: set condition for removing this extra year if necessary
extraYear <- "YES"
extraYr <- 1 #last year will be removed during numerical analysis
###################################


#Path to the study region shapefile
regShp_dir <- r"(C:\Users\alida.thiombiano\OneDrive - University of Calgary\WMO-ETSCI_Assessment_R-Workflow\PART2_ETSCIprocessing\targetRegion_shp)"
#study region shapefile name
regShp_fname <- "SMM_polygons.shp" #here, it's the combined St Mary and Milk Rivers basins
#read the region shapefile that will be used to crop the netCDF to the extent of the study region 
reg_polygon <- st_read(dsn=paste0(regShp_dir, slash, regShp_fname), quiet = TRUE) 
    #uncomment the following lines to overview target region if necessary
##plot(reg_polygon, max.plot = 1) #view the whole target region
##plot(reg_polygon[1, ], max.plot = 1, col = "transparent") #1 (2) is Milk river west (east) tributaries sub-basin; 3 is St Mary river sub-basin

# get the shapefile projection and extent
reg_crs <- st_crs(reg_polygon); reg_ext <- extent(reg_polygon)


##################################CAUTION
#There is a small portion without data for the Milk river east tributaries sub-basin
#therefore, only polygons 1 and 3 will be analyzed here
polygon_pos <- c(1, 3) #this will be used for the loop over the target regions
reg_names <- c("Milk River Basin", "St-Mary River Basin")
polygon_ext <- c(extent(reg_polygon[1, ]), extent(reg_polygon[3, ]))

nreg <- length(reg_names)

#for good visualization of maps, consider below setting for these specific target regions
latlon_breaks <- vector("list", 2)
    #Milk (East)
latlon_breaks[[1]][["lon"]] <- seq(-113.5, -109.75, 0.25) ; latlon_breaks[[1]][["lat"]] <- seq(48.5, 49.5, 0.25)
    #St-Mary
latlon_breaks[[2]][["lon"]] <- seq(-114, -112.5, 0.25) ; latlon_breaks[[2]][["lat"]] <- seq(48.50, 49.75, 0.25)

#setting labels for geographic coordinates (values called as.character() to avoid W and N beside values and give more space to the plot)
regLon_label <- "Longitude" ; regLat_label <- "Latitude Nord"

# #Note: in condition where data are available for the entire region, 
#         #to loop over the regions, this would be the format (see below commented statement)
# # name the whole region and its sub-basins (if applicable)
# reg_names <- c("SMM", "Milk River (West)", "Milk River (East)", "St. Mary")
# for (r in 1:length(reg_names)) {
#   if (r == 1) {
#     shp <- reg_polygon
#   } else {
#     shp <- reg_polygon[r-1, ]
#   }
# }
####################################

#For trend temporal analysis, temporal evolution plots, temporal long-term descriptive statistics, and temporal future changes assessment
    #set default setting for trend analysis
mkAlpha <- 0.05 #alphga value for Mann-Kendall tests (0.05 == 95% significance)
acLag   <- 1 #number of timesteps to lag data
mkDigit <- 3 #used to round slope and intercept values
MKtrendVar <- c("pvalue", "slope", "intercept") #variables names used for trend analysis' output saving
perSet <- tper_strings  #the User may select a set of time periods within the available list; e.g., c(fullPer_label, fullHist_label, fullFut_label)
tempLTVar <- c("tempLTmean", "tempLTsd") #variables names used for temporal long-term descriptive statistics calculation' output saving
    
    #create a folder to save trend analysis outputs, temporal evolution plots, temporal long-term descriptive statistics, and projected future changes results
trendFolder <- "trendResults" ; tsPlotFolder <- "tsPlots"; tempLTstats <- "tpLTstats"; futDeltFolder <- "futDeltas"
        #create folders and sub-folders
dir.create(file.path(outpath2, slash, trendFolder))
dir.create(file.path(outpath2, slash, tempLTstats))
for (r in 1:length(reg_names)) {
  dir.create(file.path(outpath2, slash, trendFolder, slash, reg_names[r]))
  dir.create(file.path(outpath2, slash, tempLTstats, slash, reg_names[r]))
  for (p in 1:length(tper_strings)) {
    dir.create(file.path(outpath2, slash, trendFolder, slash, reg_names[r], slash, tper_strings[p]))
    dir.create(file.path(outpath2, slash, tempLTstats, slash, reg_names[r], slash, tper_strings[p]))
  }
}

dir.create(file.path(outpath2, slash, tsPlotFolder))
for (r in 1:length(reg_names)) {
  dir.create(file.path(outpath2, slash, tsPlotFolder, slash, reg_names[r]))
}

dir.create(file.path(outpath2, slash, futDeltFolder))
for (r in 1:length(reg_names)) {
  dir.create(file.path(outpath2, slash, futDeltFolder, slash, reg_names[r]))
}

    #set the linewidth for both the time series and the trend lines plots
trendLw <- 0.75; TSplotLw <- 0.50
    #settings for the shading region showing data range from all data sources assessed
polyCol <- "skyblue2" #color to fill the polygon
colAlpha <- 0.4 # transparency to plot range of minimum to maximum (NO GREATER THAN 0.66)
    #color for trend lines (full historical and future time periods), and for index annual means
trendcolH <- "green4" ; trendcolF <- "orange" ; indTScol <- "black"
    #text size for temporal evolution plots
txtSize0 <- 12; txtSize1 <- 8; txtSize2 <- 8
#CAUTION: (1) regarding temporal evolution plots (step2 of section5): HARD CODING is used to aggregate subplots for each index based on the number of data sources involved!!!
        # (2) text size, colors, and other settings for graphs are based on preliminary illustrations of results.

    #For temporal scale analyses: projected future changes assessment
temporalStats <- c("mean", "sd", "min", "max")    
    #set strings to generate filenames for projected changes output Tables
annDeltasMetrics <- c("annChangeMean", "annChangeSd", "annChangeMin", "annChangeMax")
    #set color and other variables for projected change distributions boxplots
deltasBpCol <- c("cadetblue4", "chocolate1") #for near and far future respectively
bp_xaxis <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12") #because 12 data sources are assessed herein
bp_xlabel <- "Data sources" ; bp_txt1 <- 12 ; bp_txt2 <- 10 ; bp_txt3 <- 14 ; bp_hjust <- 0 
bp_fname <- "annChange_" ; ind_bpScale <- 2

#setting to derive multimodel ensemble time series which may be used to run trend analysis from annual means and medians if applicable
    #SUGGESTION to add a multimodel ensemble analysis perspective ONLY for models sharing the same future scenarios!!!
EnsTS_cNames <- c("years", "mean", "median", "min", "max")
multiModel <- "YES" #if yes, call and run the function that run trend analysis for single time series; if "NO", nothing will be done. 
    #CAUTION: note that in this script, the ensemble time series are generated from the list of data sources assessed!
              #other aggregations would required modification within the present script (e.g., considering a subset ensemble)



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>SECTION3: Indices general information>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      #CAUTION: VARIABLES GENERATED STAND AS IT IS AS LONG AS SIMILAR NETCDF OUTPUTS ARE OBTAINED FROM INDICES CALCULATION (WHICH IS EXTERNAL TO THIS SCRIPT)
          #Hence, the User should check the validity of things in this section.

#Generate the labels of the gridded ET-SCI and of the original datasets sources
    #get the full list of the original datasets used to calculate the gridded ET-SCI
climDatSource_list <- list.files(ETSCI_dir) 
    #get the full list of the ET-SCI filenames from one data source
allETSCI_files <- list.files(path=paste0(ETSCI_dir, slash, climDatSource_list[1]), pattern="*.nc", all.files=FALSE, full.names=FALSE)
    #create empty vectors for both variables
ETSCI_labels <- character(length(allETSCI_files))
climDatSource_names <- character(length(climDatSource_list))
####################
    #recall that the ET-SCI file naming should follow the Climpact netCDF naming convention
fname_struct <- strsplit(allETSCI_files[1], '_') #check that there are 6 sub-strings
####################
    #loop over the ET-SCI files to retrieve each index label (NOTE HARD CODING BASED ON EXPECTED NAMING CONVENTION)
for(i in 1:length(allETSCI_files)) {
  ETSCI_labels[i] <- paste0((strsplit(allETSCI_files[i], '_')[[1]][1]), "_", (strsplit(allETSCI_files[i], '_')[[1]][2]))
}
    #loop over the original dataset to retrieve each source name (NOTE HARD CODING BASED ON EXPECTED NAMING CONVENTION)
for(d in 1:length(climDatSource_list)) {
  climDatSource_names[d] <- paste0((strsplit(climDatSource_list[d], '_')[[1]][1]), "_", (strsplit(climDatSource_list[d], '_')[[1]][2]))
}

###############################VERY IMPORTANT NOTE RELATED TO DATA SOURCES USED!!!###############################
#The data used to implement this workflow, include low, medium and high emissions scenarios (i.e., LES, MES, HES)
    #to enable a proper ordering of the data sources, below position will be called when it comes to plot indices as driven by each data source (subplots)
    #CAUTION: HARD CODED; these positions derived from preliminary review of data sources original order when read and stored into a variable
pos_LES <- c(1, 3, 5, 9) ; pos_MES <- 7 ; pos_HES <- c(2, 6, 10, 11, 4, 8, 12)
#reorderedPos_climDatSources <- c(pos_LES, pos_MES, pos_HES)
    #Setting condition to automate the process
dsReorder <- "YES"
if (dsReorder == "YES") {
  reorderedPos_climDatSources <- c(pos_LES, pos_MES, pos_HES) #go with the re-arranged order
} else {
  reorderedPos_climDatSources <- c(1:length(climDatSource_list)) #keep the original order
}
###############################################################################################################

    #get the ET-SCI with annual time series only (Recall that ET-SCI time series are annual ('ANN') and/or monthly ('MON') ones)
        #however, add the drought indices 'spei' and 'spi' for which monthly time series alone are generated.
            #for these indices, annual aggregation should be further done to get annual time series for any month of interest.
    #using 'grep' function here (https://statisticsglobe.com/grep-grepl-r-function-example)
find_ETSCIann <- grep(paste(c("ANN","spei_MON","spi_MON"), collapse = "|"), ETSCI_labels) #to keep this applicable below, keep the files in the order they came from indices calculation (alphabetic in this case)
annETSCI_labels <- ETSCI_labels[find_ETSCIann]

#organize in categories, the list of indices with annual time series 
    #beforehand, the User should explore the indices netcdf to know the data structure (e.g., 'hw_ANN...', 'spei_MON...') using for instance Panoply (https://www.giss.nasa.gov/tools/panoply/)
    
        #heatwave and colwave indices (a pre-processing of data is necessary)
hwCw_indLabel <- "hw_ANN" #this the index label in the netcdf filename
hwCw_thresDef <- c("ecf", "ehf", "tn90", "tx90") #these are the available thresholds used to define heatwaves
hwCw_ind <- c("cwn_ecf","cwd_ecf", "cwf_ecf", "cwa_ecf", "cwm_ecf", "hwn_ehf", "hwd_ehf", "hwf_ehf", "hwa_ehf",  "hwm_ehf", 
                  "hwn_tn90", "hwd_tn90", "hwf_tn90", "hwa_tn90",  "hwm_tn90", "hwn_tx90", "hwd_tx90", "hwf_tx90", "hwa_tx90",  "hwm_tx90") #these are the variables names of the heatwave indices

        #drought indices (a pre-processing of data is necessary)
drought_indLabel <- c("spei_MON", "spi_MON") #these are the indices label in the netcdf filename

        #general indices: no pre-processing of the indices time series is necessary!
general_indLabel <- annETSCI_labels[!annETSCI_labels %in% c(hwCw_indLabel, drought_indLabel)]
general_ind <- character(length(general_indLabel))
for (g in 1:length(general_indLabel)) {
  general_ind[g] <- strsplit(general_indLabel[g], '_')[[1]][1] #CAUTION: HARD CODED (preliminary overview of indices naming)
}

#associate the unit for each index (explore https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md#appendixa)
genInd_units <- c("[days]", "[degree-days]", "[days]", "[days]", "[days]", "[°C]", "[days]", "[degree-days]", 
                 "[days]", "[degree-days]", "[days]", "[mm]", "[days]", "[days]", "[days]", "[mm]",
                 "[%]", "[mm]", "[%]", "[mm]", "[mm]", "[mm]", "[mm/day]", "[days]",
                 "[days]", "[days]", "[days]", "[days]", "[°C]", "[%]", "[%]", "[days]", 
                 "[days]", "[days]", "[°C]", "[°C]", "[°C]", "[days]", "[%]", "[events]",
                 "[%]", "[events]", "[days]", "[days]", "[%]", "[°C]", "[°C]", "[°C]", "[days]", "[days]")

droughtInd_units <- "[unitless]" #spei and spi are unitless.

    #the units of the heawave/coldwave indices (i.e., aspects) depend on the definition chosen and the units respect the following order of the heatwave aspects:
          #1. number (hwn); 2. duration (hwd); 3. frequency (hwf); 4. amplitude (hwa); 5. magnitude (hwm)
hwCwInd_units <- array(data = NA, dim = c(length(hwCw_thresDef), 5)) #CAUTION: HARD CODED (knowledge that 5 aspects of heatwave events are included in the ET-SCI)
for(k in 1:length(hwCw_thresDef)) {
  if(hwCw_thresDef[k] == "ehf" | hwCw_thresDef[k] == "ecf") {
    hwCwInd_units[k, ] <- c("[events]", "[days]", "[days]", "[°C^2]", "[°C^2]")
  } else {
    hwCwInd_units[k, ] <- c("[events]", "[days]", "[days]", "[°C]", "[°C]")
  }
}

#associate the definition of each index (explore https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md#appendixa)
hwInd_def <- c("Number of heatwaves (HWN)", "Duration of heatwaves (HWD)", "Frequency of heatwaves (HWF)", 
               "Amplitude of heatwaves (HWA)", "Magnitude of heatwaves (HWM)")
cwInd_def <- c("Number of coldwaves (CWN)", "Duration of coldwaves (CWD)", "Frequency of coldwaves (CWF)", 
               "Amplitude of coldwaves (CWA)", "Magnitude of coldwaves (CWM)")

genInd_def <- c("Consecutive dry days (CDD)", "Cooling degree days (CDDcold18)", "Cold spell duration indicator (CSDI)", "Cold spell duration indicator* (CSDI5)",
                "Consecutive wet days (CWD)", "Diurnal temperature range (DTR)", "Frost days (FD)", "Growing degree days (GDDgrow10)",
                "Growing season length (GSL)", "Heating degree days (HDDheat18)", "Ice days (ID)", "Total wet-day precipitation (PRCPTOT)", 
                "Heavy rain days (R10mm)", "Very heavy rain days (R20mm)", "Extreme rain days (R30mm)", "Total precipitation from very wet days (R95p)",
                "Contribution from very wet days (R95pTOT)", "Total precipitation from extremely wet days (R99p)", "Contribution from extremely wet days (R99pTOT)",
                "Maximum 1-day precipitation (RX1day)", "Maximum 5-days precipitation (RX5day)", "Maximum 7-days precipitation (RX7day)", 
                "Daily precipitation intensity (SDII)", "Summer days (SU)", "Days when mean temperature ≥ 10°C (TMge10)", "Days when mean temperature ≥ 5°C (TMge5)", 
                "Days when mean temperature < 10°C (TMlt10)", "Days when mean temperature < 5°C (TMlt5)",
                "Mean of daily mean temperature (TMm)", "Amount of cold nights (TN10p)", "Amount of warm nights (TN90p)", 
                "Days when minimum temperature < 2°C (TNlt2)", "Days when minimum temperature < -2°C (TNltm2)", "Days when minimum temperature < -20°C (TNltm20)", 
                "Mean of daily minimum temperature (TNm)", "Min. of daily minimum temperature (TNn)", "Max. of daily minimum temperature (TNx)", 
                "Tropical nights (TR)", "Amount of cool days (TX10p)", "Consecutive hot days and nights (TX3TN3)", "Amount of hot days (TX90p)", 
                "Consecutive cold days and nights (TXb3TNb3)", "Days when maximum temperature ≥ 30°C (TXge30)", "Days when maximum temperature ≥ 35°C (TXge35)", 
                "Fraction of days with above-median temperature (TXgt50p)", 
                "Mean of daily maximum temperature (TXm)", "Min. of daily maximum temperature (TXn)", "Max. of daily maximum temperature (TXx)", 
                "Warm spell duration indicator (WSDI)", "Warm spell duration indicator* (WSDI5)")

#for drought related indices (SPEI and SPI), specify the months of interest
      #annual aggregation will be further done for these months with respect to the adequate time scale
#################CAUTION: for any change to the time-scales and months of interest here, some adjustments may be necessary in Section4#######################
          #here, the following months and time scales are chosen to represent regional water availability 
          #over the winter (DJF), spring (MAM), summer (JJA), fall (SON), and water-year (Oct-Sept) periods (e.g. Dibike et al., 2017: https://doi.org/10.1002/joc.4912)
spei_spi_mth <- c(2, 5, 8, 11, 9) #respectively for february, may, august, november, and september annual time series!!! 
      #set the jump for annual values retrieval
spei_spi_jump <- 12 #keep this as it is
      #based on the ET-SCI netCDF for SPEI/SPI (check the netcdf files with for instance Panoply; https://www.giss.nasa.gov/tools/panoply/), 
          #data are structured in a 4D format (from original netcdf: lon x lat x time x scale)
          #scale dimension: 1 is for the 3-months time scale; 2=6-months; 3=12-months
data_dim <- 4 #spei or spi data are on the 4th dimension of the ncfile
tscale <- c(1, 3) #data are at different time-scales (here the 3- and 12-months are of interest)
      #associate definition for selected drought indices
spiSet_def <- c("SPI 3-months for February (Feb_SPI3)", "SPI 3-months for May (May_SPI3)",
                "SPI 3-months for August (Aug_SPI3)", "SPI 3-months for November (Nov_SPI3)",
                "SPI 12-months for September (Sep_SPI12)")

speiSet_def <- c("SPEI 3-months for February (Feb_SPEI3)", "SPEI 3-months for May (May_SPEI3)",
                 "SPEI 3-months for August (Aug_SEPI3)", "SPEI 3-months for November (Nov_SPEI3)",
                 "SPEI 12-months for September (Sep_SPEI12)")
      #set the label to be used in the list storing data to be retrieved
spei_set_label <- c("Feb_SPEI3", "May_SPEI3", "Aug_SPEI3", "Nov_SPEI3", "Sep_SPEI12")
spi_set_label <- c("Feb_SPI3", "May_SPI3", "Aug_SPI3", "Nov_SPI3", "Sep_SPI12")
      #get their position
tscale3_pos <- c(1:4); tscale12_pos <- 5
      #the 3- and 12-months are the time-scales of interest here 
drought_ind <- c("spei3", "spei12", "spi3", "spi12") #label to be given to data that will be retrieved and save in a list
spei_id <- c(1:2); spi_id <- c(3:4)


###########CAUTION: SPECIFIC CASES##############################################
#Keep in mind that for the following percentile-based indices describing an exceedance rate,
    # their annual time series already represent changes relative to the baseline period;
    # thus, to quantify mean or median future changes or std in future changes for these indices, 
    # derive the statistics of interest directly from the index time series over the time period of interest.
exceedRate_indices <- c("tn10p", "tn90p", "tx10p", "tx90p", "txgt50p")

#list the indices for which to round statistics following calculation
count_indices <- c("cdd", "csdi", "csdi5", "cwd", "fd", "gsl", "id", "r10mm", "r20mm", "r30mm", 
                   "su", "tmge10", "tmge5", "tmlt10", "tmlt5",  "tnlt2", "tnltm2", "tnltm20", "tr", "tx3tn3",
                   "txb3tnb3", "txge30", "txge35", "wsdi", "wsdi5",
                   "cwn_ecf","cwd_ecf", "cwf_ecf", "hwn_ehf", "hwd_ehf", "hwf_ehf",
                   "hwn_tn90", "hwd_tn90", "hwf_tn90", "hwn_tx90", "hwd_tx90", "hwf_tx90")

#for the following indices (which are all in the general category), 
      # relative changes should be calculated, and for all the rest, absolute changes.
indPos1 <- c(12, 16, 18, 20, 21, 22, 23); relativeDelInd <- general_ind[indPos1]; relativeDelInd_unit <- "[%]"
################################################################################

#combine the definitions and units of indices assessed (the indices should match the ones to be derived from step1 of section4)
    #this will be useful in section 4 for the spatial maps
        #finding position of the indices (CAUTION: this indexation required knowledge of indices order from pre-test runs!!!)
id_genInd <- c(1:10, 31:43, 48:74)
id_cwInd <- c(11:15); id_hwInd_ehf <- c(16:20); id_hwInd_tn90p <- c(21:25); id_hwInd_tx90p <- c(26:30)
id_droughtInd <- c(44:47)
        #creating a vector to store them
allInd_def <- character(length(general_indLabel) + length(drought_ind) + length(hwCw_ind))
allInd_unit <- allInd_def
        #retrieve information of interest from all groups of indices
allInd_def[id_genInd] <- genInd_def
allInd_def[id_cwInd] <- cwInd_def
allInd_def[id_hwInd_ehf] <- hwInd_def; allInd_def[id_hwInd_tn90p] <- hwInd_def; allInd_def[id_hwInd_tx90p] <- hwInd_def
allInd_def[id_droughtInd] <- drought_ind

allInd_unit[id_genInd] <- genInd_units
allInd_unit[id_cwInd] <- hwCwInd_units[1, ]
allInd_unit[id_hwInd_ehf] <- hwCwInd_units[2, ]; allInd_unit[id_hwInd_tn90p] <- hwCwInd_units[3, ]; allInd_unit[id_hwInd_tx90p] <- hwCwInd_units[4, ]
allInd_unit[id_droughtInd] <- droughtInd_units

        #considering specific sub-indices for drought category, generate similar variables for temporal analysis (done in Section5)
            #again, CAUTION: this indexation required knowledge of indices order from pre-test runs!!!
id2_genInd <- c(1:10, 31:43, 54:80)
id2_cwInd <- c(11:15)
id2_hwInd_ehf <- c(16:20); id2_hwInd_tn90p <- c(21:25); id2_hwInd_tx90p <- c(26:30)
id2_spei <- c(44:48); id2_spi <- c(49:53)

indFullSet_def <- character(length(general_indLabel) + length(hwCw_ind) + length(spei_set_label) + length(spi_set_label))
indFullSet_unit <- indFullSet_def

indFullSet_def[id2_genInd] <- genInd_def
indFullSet_def[id2_cwInd] <- cwInd_def
indFullSet_def[id2_hwInd_ehf] <- hwInd_def; indFullSet_def[id2_hwInd_tn90p] <- hwInd_def; indFullSet_def[id2_hwInd_tx90p] <- hwInd_def
indFullSet_def[id2_spei] <- speiSet_def; indFullSet_def[id2_spi] <- spiSet_def

indFullSet_unit[id2_genInd] <- genInd_units
indFullSet_unit[id2_cwInd] <- hwCwInd_units[1, ]
indFullSet_unit[id2_hwInd_ehf] <- hwCwInd_units[2, ]; indFullSet_unit[id2_hwInd_tn90p] <- hwCwInd_units[3, ]; indFullSet_unit[id2_hwInd_tx90p] <- hwCwInd_units[4, ]
indFullSet_unit[id2_spei] <- droughtInd_units; indFullSet_unit[id2_spi] <- droughtInd_units


#generate the legend title for the maps to be created (recall that these are the statistics selected herein for the spatial analysis)
legTitle <- c("mean", "median", "min", "max")

#get the position of precipitation-based indices position from the indices assessed (see list of indices derived in step2 of section4)
pos_precipInd <- c(1, 5, 32:47) #this is useful to set spatial maps color ('viridis' for these, and 'magma' for the remaining indices (i.e., temperature-based))
    #define Viridis Color Palettes for assessed indices
IndCol <- c("magma", "viridis") #respectively for temperature- and precipitation-based indices
colOrder <- c(1, -1) #1 is color default direction; -1 is reverse order
    #get position for indices that need specific color direction setting
pos_cdd <- 1 #consecutive dry days index
pos_coldInd <- c(3, 4, 7, 31) #cold temperature indices (csdi, csdi5, fd, id)

#define number of decimal for statistics computed when using "round" function
rd1 <- 0 #no decimal (count based-indices; this is also the default setting of round function without specification)
rd2 <- 2 #2 decimals (indices with units in % or mm, or degree)
rd3 <- 2 #4 decimals (SPI&SPEI indices)

#for spatial maps, time series should be averaged for each latxlon grid
gridDim <- c(1, 2) #latxlon  in this case study 
#for spatially averaged time series over the target region, time is on the 3rd dimension in this case study
avgDim <- 3 
    #CAUTION: this order (latxlonxtime) is the one we get after transforming netCDF into a raster, are further numerical analyses use that order!!!

#some general setting for spatial maps to enable easy modification
    #CAUTION: text size, colors, and other settings for graphs are based on preliminary illustrations of results
figTitle_sz <- 10 #size of the figure main title
subTitle_sz <- 7 #size for the subplot titles (by default, the subtitle is centered; so a variable can be created to set the desired position)
xyTicks_sz <- 6 #size for the x and y axis numbers
xyLabel_sz <- 8 #size for the x and y axis labels
shp_lwd <- 1 #region polygon line width
shp_col <- "darkgray" #region polygon line color
figT_hjust <- 0.5 #[figure main title justification: 0 for left-justified; 0.5 for center; 1 for right-justification]
txt_style <- "bold" #font of the title and subtitles
axis_ang <- 90 #x-axis numbers' angle
legPos <- "right" #plotted raster legend position
  #setting for optimal visual of the saved ggplot graph
fscale <- 3 #scaling (add the "scale" argument to manage the physical size of the plot, text, labels: https://stat545.com/save-figs.html)
    #NOTE: the size of the output plot can still change from one run to another; 
            #thus, after playing with this in a preliminary illustration of results, advice to set up width = ..., height = ..., units = ".." for a better control on plot size
    #output graph size
myWidth = 4200 ; myHeight = 4200 ; myUnits = "px"
ppi <- 600 #resolution

#create a subfolder to save maps illustrating long-term future changes quantified
Results_FutDeltas <- "FutChanges"
dir.create(file.path(outpath1, slash, Results_FutDeltas))

  #create within these folders, subfolders to separate outputs for each target region
      #and within these folders, create subfolders to separate results given each descriptive statistics
for (f in 1:length(reg_names)) {
  dir.create(file.path(outpath1, slash, Results_FutDeltas, slash, reg_names[f]))
  for (ff in 1:length(legTitle)) {
    dir.create(file.path(outpath1, slash, Results_FutDeltas, slash, reg_names[f], slash, legTitle[ff]))
  }
}

  #set option to generate or not the long-term climatologies maps (descriptive statistics)
getClimMaps <- "YES"

# #///////////////////////////////////////////////////////////////////////////////////////////
#EXPLORE INDICES netCDF files (OPTIONAL!!!)
    #to do so, open the R script "exploreIndices_ncfiles_Rscript.R";
    #and run code lines step by step to explore for instance, data contents and formats.
# #///////////////////////////////////////////////////////////////////////////////////////////



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>SECTION4:SPATIAL ANALYSIS>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Goal: create maps of projected future annual changes statistics (here, mean, median, min and max are the selected metrics)
        #and of ET-SCI long-term annual climatologies (similar metrics); 
            #note though that this is OPTIONAL: the function "getSpatialClimatologieMaps_func.R" is called at the appropriate section and step, if set as so!

################################################################################[EXPECTED TIME (FOR CURRENT TARGET REGIONS): LESS THAN 10MIN!]
      #step1: get for each data source and target region, the ET-SCI gridded data, and store them in a list 
annETSCI_array <- vector("list", length(climDatSource_list))

#you may assess running time for below loop
tstart1 <- Sys.time()

#loop over data source
for(s in 1:length(climDatSource_list)) {
  
  #get the ET-SCI filenames associated to each data source and retained the one for the annual scale assessment
  ETSCI_ncf <- list.files(path=paste0(ETSCI_dir, slash, climDatSource_list[s]), pattern="*.nc", all.files=FALSE, full.names=FALSE)
  annETSCI_ncf <- ETSCI_ncf[find_ETSCIann]
  
  #loop over target region
  for (r in 1:length(polygon_pos)) {
    
    #get the shapefile info
    targetReg <- reg_polygon[polygon_pos[r], ]
    reg <- reg_names[r]
    
    #loop over indices
    for(j in 1:length(annETSCI_ncf)) {
      
      #get the index '.nc' fine, label, and full filename (with path to directory)
      ncf <- annETSCI_ncf[j]
      indLabel <- paste0((strsplit(ncf, '_')[[1]][1]), "_", (strsplit(ncf, '_')[[1]][2])) #CAUTION: HARD CODED (preliminary overview of indices naming)
      ncf_name <- paste0(ETSCI_dir, slash, climDatSource_list[s], slash, ncf)
      
      #check to which category the index belongs, and proceed accordingly
          #condition1: general indices
      if(is.element(indLabel, general_indLabel) == TRUE) {
        #get the index variable name as it appears in the raster object to be created
        indVar <- strsplit(ncf, '_')[[1]][1] #CAUTION: HARD CODED (preliminary overview of indices naming)
        #Call function created at this end
        annETSCI_array[[s]][[reg]][[indVar]] <- ncfile_cropping_func(targetReg, ncf_name, indVar) #fewer arguments given for 3D ncfiles than 4D
      } #end of condition1
      
          #condition2: heatwave indices
      else if (is.element(indLabel, hwCw_indLabel) == TRUE) {
        #CAUTION: because in this case, multiple variables are embedded in the netcdf,
        #loop over the 20 heatwaves indices to process data
        for (a in 1:length(hwCw_ind)) {
          indVar <- hwCw_ind[a]
          #Call function created at this end
          annETSCI_array[[s]][[reg]][[indVar]]<- ncfile_cropping_func(targetReg, ncf_name, indVar) #fewer arguments given for 3D ncfiles than 4D
        } #end of the heatwave list of indices
      } #end of condition2
      
          #condition3: drought indices
      else if (is.element(indLabel, drought_indLabel) == TRUE) {
        indVar <- strsplit(ncf, '_')[[1]][1] #CAUTION: HARD CODED (preliminary overview of indices naming)
        #set condition with respect to the drought index
        if (indVar == "spei") {
          indset <- drought_ind[spei_id]
        }
        else if (indVar == "spi") {
          indset <- drought_ind[spi_id]
        }
        #retrieve data for selected time-scales
        for (b in 1:length(tscale)) {
          #Call function created at this end
          annETSCI_array[[s]][[reg]][[indset[b]]] <- ncfile_cropping_func(targetReg, ncf_name, indVar, data_dim, tscale[b]) #more arguments given for 4D nfciles
        } #end of time-scales
      } #end of condition3
      
    } #end of indices netcdf files
  } #end of target regions
} #end of data source
tend1 <- Sys.time()


################################################################################[EXPECTED TIME (FOR CURRENT TARGET REGIONS): LESS THAN 2MIN!]
      #step2: for each index, retrieve annual time series (gridded and spatially average) for the time period of interest, 
              #then calculate the long-term selected statistics (mean, median, min and max) over each one of them
#get the names of the two different cropped datasets created in above step
cropDat_list <- names(annETSCI_array[[1]][[reg]][[indVar]]) #CAUTION: HARD CODED (based on programming structure adopted herein)

#get the names of all the indices in the previous created list
annInd_list <- names(annETSCI_array[[1]][[reg]]) #CAUTION: HARD CODED (based on programming structure adopted herein)

#create a list to save indices gridded annual time series with respect to data source, target region, and time period
annETSCI_gridTS <- vector("list", length(climDatSource_list))

#create a list to save long-term gridded statistics calculated
annETSCI_gridStat <- vector("list", length(climDatSource_list))

#create a list to save indices spatially averaged annual time series with respect to data source, target region, and time period
annETSCI_spAvgTS <- vector("list", length(climDatSource_list))

#create a list to save gridded long-term projected future annual changes stats
annETSCI_deltas_gridStat <- vector("list", length(climDatSource_list))

#create a list to save indices min., max., mean and median basin-wide single annual time series as obtained from gridded data (ONLY FOR THE FULL TIME PERIOD)
      #to further generate same colobar limit for original driven-long-term descriptive statistics' spatial maps
annETSCI_spOrigDataRange <- vector("list", length(climDatSource_list))

#create a list to save indices min., max., mean and median basin-wide single annual change values as obtained from near and far future Deltas gridded data 
      #to further generate same colobar limit for changes' spatial maps
annETSCI_spDeltasDataRange <- vector("list", length(climDatSource_list))


#you may assess running time for below loop
tstart2 <- Sys.time()

#loop over data sources, target regions, indices and time periods
for (ss in 1:length(annETSCI_array)) {
  
  #target regions
  for (rr in 1:length(reg_names)) {
    reg <- reg_names[rr] #region name
    
    #indices
    for (ii in 1:length(annInd_list)) {
      
      #get index label and data
      indVar2 <- annInd_list[ii]
      
      #get arrays (CAUTION: HARD CODING BASED ON DATA FORMATTING KNOWLEDGE!!!)
      datCrop <- annETSCI_array[[ss]][[reg]][[indVar2]][[cropDat_list[1]]]
      datCropMask <- annETSCI_array[[ss]][[reg]][[indVar2]][[cropDat_list[2]]]
        #get dimension to check for time length (important for drought indices monthly time series aggregation)
      arrayDim <- dim(datCrop) #here, lat x lon x time
      
      #set conditions for values to be round regarding the index category
      if (is.element(indVar2, count_indices) == TRUE) {
        rd <- rd1
      } else if (is.element(indVar2, drought_ind) == TRUE) {
        rd <- rd3
      } else {
        rd <- rd2
      }
      
      #time periods
      for (p in 1:length(tperiods)) {
        
        #derive the years position of each period of interest in the study full time coverage (HARD CODED: this period is first in the created list!)
        per <- which(tperiods[[1]] %in% tperiods[[p]]) #period interval for data retrieval | CAUTION: HARD CODED (based on knowledge of the created list' structure)
        perName <- tper_strings[p]
        
        #set conditions for group of indices
            #conditions1&2: general indices & heatwave indices
        if (is.element(indVar2, general_ind) == TRUE || is.element(indVar2, hwCw_ind) == TRUE) {
          #gridded annual time series
          annETSCI_gridTS[[ss]][[reg]][[indVar2]][[perName]] <- datCrop[, , per]
          #long-term gridded annual statistics (Call function created at this end)
          dat0 <- datCrop[, , per]
          annETSCI_gridStat[[ss]][[reg]][[indVar2]][[perName]] <- generateSpatialStats_func(dat0, gridDim, rd)
          #spatially averaged annual time series (Call function created at this end)
          dat00 <- datCropMask[, , per]
          annETSCI_spAvgTS[[ss]][[reg]][[indVar2]][[perName]] <- generateSpatialAvgTS_func(dat00, avgDim, rd)
          
          ############################################
          #get min, max, mean and median single annual time series for each index, to further generate original data range from the cropped (but not masked) gridded annual time series
                #CAUTION: HARD CODED (interest for the study full time period which is first in the created list)
          if (p == 1) {
            annETSCI_spOrigDataRange[[ss]][[reg]][[indVar2]][[perName]] <- generateSpatialDataRange_func(dat0, avgDim, rd)
          }
          ###########################################
          
          #clear variable to avoid error
          rm(dat0, dat00)
             
        } #end of conditions1&2
        
            #condition3: drought indices
        else if (is.element(indVar2, drought_ind) == TRUE) {
          
          #set condition with respect to the drought index (CAUTION: HARD CODING BASED ON SELECTED TIME SCALES FOR SPI AND SPEI!!!)
          if (indVar2 == "spei3" || indVar2 == "spi3") {
            #loop over the set of month for which to compile annual time series
            for (t in 1:length(tscale3_pos)) {
              if (indVar2 == "spei3") {
                indvar3 <- spei_set_label[t]
              } else if (indVar2 == "spi3") {
                indvar3 <- spi_set_label[t]
              }
              #set the sequence interval (recall to remove the extra year if applicable)
              if (extraYear == "YES") {
                timeSeq <- head(seq(from=spei_spi_mth[t], to=arrayDim[avgDim], by=spei_spi_jump), - extraYr)
              } else {
                timeSeq <- seq(from=spei_spi_mth[t], to=arrayDim[avgDim], by=spei_spi_jump)
              }
              #get the data for each specific month
              dat1 <- datCrop[, , timeSeq]
              dat2 <- datCropMask[, , timeSeq] #for each month of interest
              #gridded annual time series
              annETSCI_gridTS[[ss]][[reg]][[indVar2]][[indvar3]][[perName]] <- dat1[, , per]
              #long-term gridded annual statistics (Call function created at this end)
              dat0 <- dat1[, , per]
              annETSCI_gridStat[[ss]][[reg]][[indVar2]][[indvar3]][[perName]] <- generateSpatialStats_func(dat0, gridDim, rd)
              #spatially averaged annual time series (Call function created at this end)
              dat00 <- dat2[, , per]
              annETSCI_spAvgTS[[ss]][[reg]][[indVar2]][[indvar3]][[perName]] <- generateSpatialAvgTS_func(dat00, avgDim, rd)
              
              ############################################
              #get min, max, mean and median single annual time series for each index, to further generate original data range from the cropped (but not masked) gridded annual time series
                  #CAUTION: HARD CODED (interest for the study full time period which is first in the created list)
              if (p == 1) {
                annETSCI_spOrigDataRange[[ss]][[reg]][[indVar2]][[indvar3]][[perName]] <- generateSpatialDataRange_func(dat0, avgDim, rd)
              }
              ###########################################
              
              #clear variable to avoid error
              rm(timeSeq, dat1, dat2, dat0, dat00)
              
            } #end of 3-months time scale
          } #end of 3-months time scale
          
          else if (indVar2 == "spei12" || indVar2 == "spi12") {
            for (t in 1:length(tscale12_pos)) {
              #setting things for automatic processing when there are more than on index here
              if (length(tscale12_pos) == 1) {
                tt <- tscale12_pos
              } else {
                tt <- tscale12_pos[t]
              }
              #consider here as well the specificity of the drought indices
              if (indVar2 == "spei12") {
                indvar3 <- spei_set_label[tt]
              } else if (indVar2 == "spi12") {
                indvar3 <- spi_set_label[tt]
              }
              #set the sequence interval (recall to remove the extra year if applicable)
              if (extraYear == "YES") {
                timeSeq <- head(seq(from=spei_spi_mth[tt], to=arrayDim[avgDim], by=spei_spi_jump), - extraYr)
              } else {
                timeSeq <- seq(from=spei_spi_mth[tt], to=arrayDim[avgDim], by=spei_spi_jump)
              }
              #get the data for each specific month
              dat1 <- datCrop[, , timeSeq]
              dat2 <- datCropMask[, , timeSeq] #for each month of interest
              #gridded annual time series
              annETSCI_gridTS[[ss]][[reg]][[indVar2]][[indvar3]][[perName]] <- dat1[, , per]
              #long-term gridded annual statistics (Call function created at this end)
              dat0 <- dat1[, , per]
              annETSCI_gridStat[[ss]][[reg]][[indVar2]][[indvar3]][[perName]] <- generateSpatialStats_func(dat0, gridDim, rd)
              #spatially averaged annual time series (Call function created at this end)
              dat00 <- dat2[, , per]
              annETSCI_spAvgTS[[ss]][[reg]][[indVar2]][[indvar3]][[perName]] <- generateSpatialAvgTS_func(dat00, avgDim, rd)
              
              ############################################
              #get min, max, mean and median single annual time series for each index, to further generate original data range from the cropped (but not masked) gridded annual time series
                  #CAUTION: HARD CODED (interest for the study full time period which is first in the created list)
              if (p == 1) {
                annETSCI_spOrigDataRange[[ss]][[reg]][[indVar2]][[indvar3]][[perName]] <- generateSpatialDataRange_func(dat0, avgDim, rd)
              }
              ###########################################
              
              #clear variable to avoid error
              rm(timeSeq, dat1, dat2, dat0, dat00)
              
              } #end of 12-months time scale
          } #end of 12-months time scale
          
        } #end of condition3 

        #clear variable to avoid error
        rm(per)
        
      } #end of time periods

      #clear variable after each index to avoid error
      rm(datCrop, datCropMask)
      
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Quantifying projected future annual changes>>>>>>>>>>>>>>>>>>>>>    
      #having now data from all time periods, quantify projected future annual changes relative to baseline conditions
          #loop over each grid point (lon x lat), and calculate the changes stats (mean, median, min, max) using the gridded annual time series
      
              #condtion1: drought indices specific case (i.e., sub-group of indices from the main index); again, CAUTION: HARD CODING BASED ON SELECTED TIME SCALES FOR SPI AND SPEI!!!
      if(is.element(indVar2, drought_ind) == TRUE) {
        
        #get sub-goup of indices names with respect to SPEI or SPI
        if (indVar2 == "spei3" || indVar2 == "spi3") {
          drought_subset <- names(annETSCI_gridTS[[ss]][[reg]][[indVar2]])
        } #end of 3-months time scale indices (here 4 indices)
        else if (indVar2 == "spei12" || indVar2 == "spi12") {
          drought_subset <- names(annETSCI_gridTS[[ss]][[reg]][[indVar2]])
        } #end of 12-months time scale indices (here 1 index)
        
        #now, loop over this subset
        for (w in 1:length(drought_subset)) {
          indVar4 <- drought_subset[w]
          
          #get input data regarding the index category (note that here, data dimension is lat x lon x time)
          dat_ref <- annETSCI_gridTS[[ss]][[reg]][[indVar2]][[indVar4]][[refLabel]] #baseline
          dat_Fut1 <- annETSCI_gridTS[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]] #near future
          dat_Fut2 <- annETSCI_gridTS[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]] #far future
          
          #create empty matrices to store long-term statistics over the future projected changes (two periods here)
          outMat_mean_F1 <- array(data = NA, dim = c(dim(dat_ref)[1], dim(dat_ref)[2])) #CAUTION: HARD CODED (preliminary review of data structure)
          outMat_median_F1 <- outMat_mean_F1
          outMat_min_F1 <- outMat_mean_F1
          outMat_max_F1 <- outMat_mean_F1
          
          outMat_mean_F2 <- outMat_mean_F1
          outMat_median_F2 <- outMat_mean_F1
          outMat_min_F2 <- outMat_mean_F1
          outMat_max_F2 <- outMat_mean_F1
          
          #iterate over data dimensions to calculate long-term min, max, mean and median changes
          for (lat in 1:dim(dat_ref)[1]) {
            for (lon in 1:dim(dat_ref)[2]) {
              
              #create empty matrix to store annual TS per grid: col1=ref; col2=fut1; col3=fut2; col4=fut1Deltas; col5=fut2Deltas
              outMat <- array(data = NA, dim = c(dim(dat_ref)[3], 5)) #CAUTION: HARD CODED (preliminary review of data structure + selected periods for future change assessment)
              outMat[, 1] <-  dat_ref[lat, lon, ]
              outMat[, 2] <-  dat_Fut1[lat, lon, ]
              outMat[, 3] <-  dat_Fut2[lat, lon, ]
              
              ##############################################
              #calculate the baseline mean condition (single value)
              ref_meanVal <- mean(outMat[, 1], na.rm=TRUE)
              ##############################################
              
                #loop over years
              for (y in 1:dim(dat_ref)[3]) {
                outMat[y, 4] <- outMat[y, 2] - ref_meanVal
                outMat[y, 5] <- outMat[y, 3] - ref_meanVal
              } #end of years
              
              #long-term stats over F1
              outMat_mean_F1[lat, lon] <- round(mean(outMat[, 4], na.rm = TRUE), rd)
              outMat_median_F1[lat, lon] <- round(median(outMat[, 4], na.rm = TRUE), rd)
              outMat_min_F1[lat, lon] <- round(min(outMat[, 4], na.rm = TRUE), rd)
              outMat_max_F1[lat, lon] <- round(max(outMat[, 4], na.rm = TRUE), rd)
              
              #long-term stats over F2
              outMat_mean_F2[lat, lon] <- round(mean(outMat[, 5], na.rm = TRUE), rd)
              outMat_median_F2[lat, lon] <- round(median(outMat[, 5], na.rm = TRUE), rd)
              outMat_min_F2[lat, lon] <- round(min(outMat[, 5], na.rm = TRUE), rd)
              outMat_max_F2[lat, lon] <- round(max(outMat[, 5], na.rm = TRUE), rd)
              
            } #end of longitude grids
          } #end of latitude grids
          
          #store long-term statistics of projected annual changes over selected periods
          annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]][["mean"]] <- outMat_mean_F1
          annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]][["median"]] <- outMat_median_F1
          annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]][["min"]] <- outMat_min_F1
          annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]][["max"]] <- outMat_max_F1
          
          annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]][["mean"]] <- outMat_mean_F2
          annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]][["median"]] <- outMat_median_F2
          annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]][["min"]] <- outMat_min_F2
          annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]][["max"]] <- outMat_max_F2
          
          #######################################
          #derive the range of values from the gridded annual changes
          annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]][["meanRange"]] <- c(min(outMat_mean_F1, na.rm = TRUE), max(outMat_mean_F1, na.rm = TRUE))
          annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]][["medianRange"]] <- c(min(outMat_median_F1, na.rm = TRUE), max(outMat_median_F1, na.rm = TRUE))
          annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]][["minRange"]] <- c(min(outMat_min_F1, na.rm = TRUE), max(outMat_min_F1, na.rm = TRUE))
          annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]][["maxRange"]] <- c(min(outMat_max_F1, na.rm = TRUE), max(outMat_max_F1, na.rm = TRUE))
          
          annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]][["meanRange"]] <- c(min(outMat_mean_F2, na.rm = TRUE), max(outMat_mean_F2, na.rm = TRUE))
          annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]][["medianRange"]] <- c(min(outMat_median_F2, na.rm = TRUE), max(outMat_median_F2, na.rm = TRUE))
          annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]][["minRange"]] <- c(min(outMat_min_F2, na.rm = TRUE), max(outMat_min_F2, na.rm = TRUE))
          annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]][["maxRange"]] <- c(min(outMat_max_F2, na.rm = TRUE), max(outMat_max_F2, na.rm = TRUE))
          #######################################
          
          #clear variable to avoid error
          rm(dat_ref, dat_Fut1, dat_Fut2)
          
        } #end of drought sub-index
        
      } #end of conditon1
 
              #condition2: all other indices
      else {
        dat_ref <- annETSCI_gridTS[[ss]][[reg]][[indVar2]][[refLabel]] #baseline
        dat_Fut1 <- annETSCI_gridTS[[ss]][[reg]][[indVar2]][[fut1Label]] #near future
        dat_Fut2 <- annETSCI_gridTS[[ss]][[reg]][[indVar2]][[fut2Label]] #far future
        
        #create empty matrices to store long-term statistics over the future projected changes (two periods here)
        outMat_mean_F1 <- array(data = NA, dim = c(dim(dat_ref)[1], dim(dat_ref)[2])) #CAUTION: HARD CODED (preliminary review of data structure)
        outMat_median_F1 <- outMat_mean_F1
        outMat_min_F1 <- outMat_mean_F1
        outMat_max_F1 <- outMat_mean_F1
        
        outMat_mean_F2 <- outMat_mean_F1
        outMat_median_F2 <- outMat_mean_F1
        outMat_min_F2 <- outMat_mean_F1
        outMat_max_F2 <- outMat_mean_F1
        
        #iterate over data dimensions to calculate long-term min, max, mean and median changes
        for (lat in 1:dim(dat_ref)[1]) {
          for (lon in 1:dim(dat_ref)[2]) {
            
            #create empty matrix to store annual TS per grid: col1=ref; col2=fut1; col3=fut2; col4=fut1Deltas; col5=fut2Deltas
            outMat <- array(data = NA, dim = c(dim(dat_ref)[3], 5)) #CAUTION: HARD CODED (preliminary review of data structure + selected periods for future change assessment)
            outMat[, 1] <-  dat_ref[lat, lon, ]
            outMat[, 2] <-  dat_Fut1[lat, lon, ]
            outMat[, 3] <-  dat_Fut2[lat, lon, ]
            
            ##############################################
            #calculate the baseline mean condition (single value)
            ref_meanVal <- mean(outMat[, 1], na.rm=TRUE)
            ##############################################
            
            #set the conditions for absolute changes, relative changes, and exceedance indices
                #conditionA: relative changes (i.e., precipitation-based indices in mm)
            if (is.element(indVar2, relativeDelInd) == TRUE) {
              #loop over years
              for (y in 1:dim(dat_ref)[3]) {
                outMat[y, 4] <- ((outMat[y, 2] - ref_meanVal) / ref_meanVal)*100
                outMat[y, 5] <- ((outMat[y, 3] - ref_meanVal) / ref_meanVal)*100
              } #end of years
            } #end of conditionA
            
                #conditionB: changes from exceedance indices
            else if (is.element(indVar2, exceedRate_indices) == TRUE) {
              outMat[, 4] <- outMat[, 2]
              outMat[, 5] <- outMat[, 3]
            } #end of conditionB
            
                #conditionC: absolute changes
            else {
              #loop over years
              for (y in 1:dim(dat_ref)[3]) {
                outMat[y, 4] <- outMat[y, 2] - ref_meanVal
                outMat[y, 5] <- outMat[y, 3] - ref_meanVal
              } #end of years
            } #end of conditionC
            
            #long-term stats over F1
            outMat_mean_F1[lat, lon] <- round(mean(outMat[, 4], na.rm = TRUE), rd)
            outMat_median_F1[lat, lon] <- round(median(outMat[, 4], na.rm = TRUE), rd)
            outMat_min_F1[lat, lon] <- round(min(outMat[, 4], na.rm = TRUE), rd)
            outMat_max_F1[lat, lon] <- round(max(outMat[, 4], na.rm = TRUE), rd)
            
            #long-term stats over F2
            outMat_mean_F2[lat, lon] <- round(mean(outMat[, 5], na.rm = TRUE), rd)
            outMat_median_F2[lat, lon] <- round(median(outMat[, 5], na.rm = TRUE), rd)
            outMat_min_F2[lat, lon] <- round(min(outMat[, 5], na.rm = TRUE), rd)
            outMat_max_F2[lat, lon] <- round(max(outMat[, 5], na.rm = TRUE), rd)

          } #end of longitude grids
        } #end of latitude grids
        
        #store long-term statistics of projected annual changes over selected periods
        annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[fut1Label]][["mean"]] <- outMat_mean_F1
        annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[fut1Label]][["median"]] <- outMat_median_F1
        annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[fut1Label]][["min"]] <- outMat_min_F1
        annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[fut1Label]][["max"]] <- outMat_max_F1
        
        annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[fut2Label]][["mean"]] <- outMat_mean_F2
        annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[fut2Label]][["median"]] <- outMat_median_F2
        annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[fut2Label]][["min"]] <- outMat_min_F2
        annETSCI_deltas_gridStat[[ss]][[reg]][[indVar2]][[fut2Label]][["max"]] <- outMat_max_F2
        
        #######################################
        #derive the range of values from the gridded annual changes
        annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut1Label]][["meanRange"]] <- c(min(outMat_mean_F1, na.rm = TRUE), max(outMat_mean_F1, na.rm = TRUE))
        annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut1Label]][["medianRange"]] <- c(min(outMat_median_F1, na.rm = TRUE), max(outMat_median_F1, na.rm = TRUE))
        annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut1Label]][["minRange"]] <- c(min(outMat_min_F1, na.rm = TRUE), max(outMat_min_F1, na.rm = TRUE))
        annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut1Label]][["maxRange"]] <- c(min(outMat_max_F1, na.rm = TRUE), max(outMat_max_F1, na.rm = TRUE))
        
        annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut2Label]][["meanRange"]] <- c(min(outMat_mean_F2, na.rm = TRUE), max(outMat_mean_F2, na.rm = TRUE))
        annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut2Label]][["medianRange"]] <- c(min(outMat_median_F2, na.rm = TRUE), max(outMat_median_F2, na.rm = TRUE))
        annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut2Label]][["minRange"]] <- c(min(outMat_min_F2, na.rm = TRUE), max(outMat_min_F2, na.rm = TRUE))
        annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut2Label]][["maxRange"]] <- c(min(outMat_max_F2, na.rm = TRUE), max(outMat_max_F2, na.rm = TRUE))
        #######################################
        
        #clear variable to avoid error
        rm(dat_ref, dat_Fut1, dat_Fut2)
        
      } #end of conditon2
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>end of projected changes quantification>>>>>>>>>>>
    
    } #end of indices
  } #end of regions
} #end of data sources
tend2 <- Sys.time()


#################################################
#get indices' values range common to all data sources, time periods & target regions for spatial maps colorbar limits setting 
allIndRangeList <- getIndicesDataRange_forSpatialMaps_func()

indRange_spClim <- allIndRangeList[["indRange_spClim"]] #used in a function (that generate spatial maps of long-term  climatologies) to define colobar limit
indRange_spDel <- allIndRangeList[["indRange_spDel"]] #used below, to define colobar limit for spatial maps of projected long-term annual changes

rm(allIndRangeList) #remove original variable
#################################################


################################################################################[EXPECTED TIME (FOR CURRENT TARGET REGIONS): ABOUT 40MIN!]
    #step3: create maps of ET-SCI long-term projected future annual changes (here, mean, median, min and max are the selected metrics)
            #use ggplot2 to map the raster object (https://stackoverflow.com/questions/33227182/how-to-set-use-ggplot2-to-map-a-raster)

#using Viridis Color Palettes predefined for temperature- and precipitation-based indices, 
  #create empty variable and store the sequence of the colors and their palette orientation for all assessed indices
annInd_list_colPal <- character(length(annInd_list))
annInd_list_coldir <- matrix(data = NA, nrow = 1, ncol = length(annInd_list))
  #loop
for (n in 1:length(annInd_list)) {
  if (is.element(n, pos_precipInd) == TRUE) {
    annInd_list_colPal[n] <- IndCol[2]
    annInd_list_coldir[n] <- colOrder[2]
  } else if (is.element(annInd_list[n], annInd_list[pos_coldInd]) == TRUE) {
    annInd_list_colPal[n] <- IndCol[1]
    annInd_list_coldir[n] <- colOrder[2]
  } else {
    annInd_list_colPal[n] <- IndCol[1]
    annInd_list_coldir[n] <- colOrder[1]
  }
} #end of loop

#loop over indices, regions, future periods assessed, and set of change metrics, 
        #and retrieve gridded data for all data sources, then generate maps with subplots (each subplot illustrating one data source)

    #get from one file, info necessary to iteration (CAUTION: HARD CODING based on knowledge of selected time scales for SPEI and SPI!!!)
futPers <- names(annETSCI_deltas_gridStat[[1]][[reg_names[1]]][[annInd_list[1]]])
statMetrics <- names(annETSCI_deltas_gridStat[[1]][[reg_names[1]]][[annInd_list[1]]][[futPers[1]]])
spei_scale3 <- names(annETSCI_deltas_gridStat[[1]][[reg_names[1]]][["spei3"]])
spei_scale12 <- names(annETSCI_deltas_gridStat[[1]][[reg_names[1]]][["spei12"]])
spi_scale3 <- names(annETSCI_deltas_gridStat[[1]][[reg_names[1]]][["spi3"]])
spi_scale12 <- names(annETSCI_deltas_gridStat[[1]][[reg_names[1]]][["spi12"]])

#you may assess running time for below loop
tstart3 <- Sys.time()

    #loop
for (iii in 1:length(annInd_list)) {
  indexLabel <- annInd_list[iii] #index label
  colPal <- annInd_list_colPal[iii] #color palette
  colDir <- annInd_list_coldir[iii] #color scale orientation
  
  for (h in 1:length(reg_names)) {
    xlon_breaks <- latlon_breaks[[h]][["lon"]]
    ylat_breaks <- latlon_breaks[[h]][["lat"]]
    regShp <- reg_polygon[polygon_pos[h], ]
    
    for (pp in 1:length(futPers)) {
      
      for (m in 1:length(statMetrics)) {
        
        if (indexLabel == annInd_list[pos_cdd]) {
          if (statMetrics[m] == "min") {
            colDir <- colOrder[2] #CAUTION: HARD CODING!!!
          } else {
            colDir <- colOrder[1] #CAUTION: HARD CODING!!!
          }
        }
        
        #folder where to store maps
        outFolder <- paste0(outpath1, slash, Results_FutDeltas, slash, reg_names[h], slash, legTitle[m])
        
        ######################        
        #set condition1: drought indices (CAUTION: HARD CODING BASED ON SELECTED TIME SCALES FOR SPI AND SPEI!!!)
        if(is.element(indexLabel, drought_ind) == TRUE) {
          
          #set legend title
          legTitle_full <- paste0(legTitle[m], " ", greeks("Delta"), "\n", droughtInd_units)
          
          #define figure title regarding SPEI or SPI
          if (indexLabel == "spei3") {
            indDef <- speiSet_def[tscale3_pos]
            drought_subset2 <- spei_scale3
          } else if (indexLabel == "spei12") {
            indDef <- speiSet_def[tscale12_pos]
            drought_subset2 <- spei_scale12
          } else if (indexLabel == "spi3") {
            indDef <- spiSet_def[tscale3_pos]
            drought_subset2 <- spi_scale3
          } else if(indexLabel == "spi12") {
            indDef <- spiSet_def[tscale12_pos]
            drought_subset2 <- spi_scale12
          }
          
          #loop over this subset
          for (ww in 1:length(drought_subset2)) {
            figName <- paste0(outFolder, slash, drought_subset2[ww], "_", futPers[pp], ".png")
            figTitle <- indDef[ww]
            colorbarLim <- indRange_spDel[[indexLabel]][[drought_subset2[ww]]][[statMetrics[m]]] #common limit for colorbar
            
            #and now over data sources (CAUTION: HERE, THE ORIGINAL ORDER "length(annETSCI_deltas_gridStat)" OF DATA SOURCES ARE REVISED TO GATHER low, medium and high emissions)
            for (dd in 1:length(reorderedPos_climDatSources)) {
              newPos <- reorderedPos_climDatSources[dd]
              #matrix
              inDat <- annETSCI_deltas_gridStat[[newPos]][[reg_names[h]]][[indexLabel]][[drought_subset2[ww]]][[futPers[pp]]][[statMetrics[m]]]
              #transform matrix into a raster
              inDat_rast <- raster(inDat, crs=reg_crs) #give region shapefile projection and extent
              extent(inDat_rast) <- c(polygon_ext[[h]]@xmin, polygon_ext[[h]]@xmax, polygon_ext[[h]]@ymin, polygon_ext[[h]]@ymax)
              #set the title for each subplot
              names(inDat_rast) <- climDatSource_names[newPos] #CAUTION: this does not display labels as expected; thus, labels to be arranged.
              #combined output raster for all data sources to generate a rasterStack
              if (dd==1) {
                output_ras <- inDat_rast
              } else {
                output_ras <- stack(output_ras, inDat_rast) #append the following raster objects to the 1st one
              }
              
              #clear variables to avoid errors
              rm(inDat, inDat_rast)
              
            } #end of data sources
            
            #generating maps for each sub-index (Call function created at this end)
            generateSpatialMaps_func(output_ras, colorbarLim)
            
            #clear variables to avoid errors
            rm(output_ras, colorbarLim)
            
          } #end of drought sub-index
        } #end of condition1
        
        ######################
        #set condition2: other indices (i.e., (is.element(indexLabel, general_ind) == TRUE || is.element(indexLabel, hwCw_ind) == TRUE))
        else {
          #set legend title, figure saving name, figure title
          if (is.element(indexLabel, relativeDelInd) == TRUE) {
            legTitle_full <- paste0(legTitle[m], " ", greeks("Delta"), "\n", relativeDelInd_unit)
          } else {
            legTitle_full <- paste0(legTitle[m], " ", greeks("Delta"), "\n", allInd_unit[iii])
          }
          
          figName <- paste0(outFolder, slash, indexLabel, "_", futPers[pp], ".png")
          figTitle <- allInd_def[iii]
          colorbarLim <- indRange_spDel[[indexLabel]][[statMetrics[m]]] #common limit for colorbar
          
          #and now over data sources
          for (dd in 1:length(reorderedPos_climDatSources)) {
            newPos <- reorderedPos_climDatSources[dd]
            #matrix
            inDat <- annETSCI_deltas_gridStat[[newPos]][[reg_names[h]]][[indexLabel]][[futPers[pp]]][[statMetrics[m]]]
            #transform matrix into a raster
            inDat_rast <- raster(inDat, crs=reg_crs) #give region shapefile projection and extent
            extent(inDat_rast) <- c(polygon_ext[[h]]@xmin, polygon_ext[[h]]@xmax, polygon_ext[[h]]@ymin, polygon_ext[[h]]@ymax)
            #set the title for each subplot
            names(inDat_rast) <- climDatSource_names[newPos] #CAUTION: this does not display labels as expected; thus, labels to be arranged.
            #combined output raster for all data sources to generate a rasterStack
            if (dd==1) {
              output_ras <- inDat_rast
            } else {
              output_ras <- stack(output_ras, inDat_rast) #append the following raster objects to the 1st one
            }
            
            #clear variables to avoid errors
            rm(inDat, inDat_rast)
            
          } #end of data sources
          
          #generating maps for each sub-index (Call function created at this end)
          generateSpatialMaps_func(output_ras, colorbarLim)
          
          #clear variables to avoid errors
          rm(output_ras, colorbarLim)
          
        } #end of condition2
        ######################

      } #end of descriptive stats
    } #end of assessed periods
  } #end of regions
} #end of indices        
tend3 <- Sys.time()


################################################################################[EXPECTED TIME (FOR CURRENT TARGET REGIONS): ABOUT 4H!]
    #step4 (OPTIONAL): create maps of ET-SCI long-term annual climatologies (here, mean, median, min and max are the selected metrics)
            #call the function that generate such results, if set-up so: "getSpatialClimatologieMaps_func.R" 
#you may assess running time for below loop
tstart <- Sys.time()
if (getClimMaps == "YES") {
  getSpatialClimatologieMaps_func()
}
tend <- Sys.time()

###############################################################################



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>SECTION5:TEMPORAL ANALYSIS>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

################################################################################[EXPECTED TIME (FOR CURRENT TARGET REGIONS): ABOUT 15MIN!]
#step1: assess trend over the full time coverage, the historical and future periods, and any other time period of interest;
          #additionally, from input datasets, re-organize spatially average time series regarding the full set of indices assessed here
          #then, retrieve mean, median, min, max annual time series from all data sources for further analyses.

    #create lists to save results of the trend analysis for each index with respect to data source, target region, and time period
MKtrend_pval <- list(); MKtrend_slope <- list(); MKtrend_interc <- list()
    #create lists to save annual time series for each data source, then the time series of the descriptive statistics generated from all data sources
allSim.allInd_spAvgTS <- list(); allSim.allInd_spAvgTS_stat <- list()

    #################for temporal long-term (LT) descriptive statistics#########
tempLT_mean <- list() ; tempLT_sd <- list()
    ############################################################################

    #you may assess running time for below loop
tstart4 <- Sys.time()

    #loop over the indices, regions, selected periods and data sources, 
        #to get input data, run trend analysis, and save results from all data sources, but for one region, period, and one index at the time
        #but also to aggregate all indices annual time series necessary for further analyses
for (i in 1:length(annInd_list)) {
  ind0 <- annInd_list[i]
  
  #set conditions for values to be round regarding the index category
  if (is.element(ind0, count_indices) == TRUE) {
    rd <- rd1
  } else if (is.element(ind0, drought_ind) == TRUE) {
    rd <- rd3
  } else {
    rd <- rd2
  }
  
  for (r in 1:length(reg_names)) {
    
    for (v in 1:length(perSet)) {
      
      #set conditions regarding indices group, and get indices label & input data looping through all data sources
      if(is.element(ind0, drought_ind) == TRUE) {
        indList <- names(annETSCI_spAvgTS[[d]][[reg_names[r]]][[ind0]])
        
        for (w in 1:length(indList)) {
          ind1 <- indList[w] #index label of interest
          
          #empty matrix to save data retrieved for further analyses
          allSimTS <- array(data = NA, dim = c(length(tperiods[[v]]), length(climDatSource_list)))

          #empty matrices to save pvalue, slope and intercept outputs
          MK_pValue  <- array(data = NA, dim = c(1, length(climDatSource_list))) 
          MK_slope <- MK_pValue; MK_intercept <- MK_pValue
          
          ###################################
          #empty matrices to save long-term (lt) temporal mean and standard deviation (Sd) statistics for each index and data sources
          tempSpAvgTS_ltMean <- array(data = NA, dim = c(1, length(climDatSource_list)))
          tempSpAvgTS_ltSd <- tempSpAvgTS_ltMean
          ###################################
          
          #loop over all data sources (CAUTION: HERE, THE ORIGINAL ORDER "length(annETSCI_deltas_gridStat)" OF DATA SOURCES ARE REVISED TO GATHER low, medium and high emissions)
          for (d in 1:length(reorderedPos_climDatSources)) {
            newPos <- reorderedPos_climDatSources[d]
            
            inputDat <- annETSCI_spAvgTS[[newPos]][[reg_names[r]]][[ind0]][[ind1]][[perSet[v]]]
            #call and run the trend analysis function
            mkTest <- general_test_mk_prewhiten(inputDat, acLag, mkAlpha)
            MK_pValue[1, d] <- mkTest[[1]]
            MK_slope[1, d] <- round(mkTest[[2]], mkDigit)
            MK_intercept[1, d] <- round(mkTest[[3]], mkDigit)
            #now get data for further analyses
            allSimTS[, d] <- inputDat
            
            #########get temporal long-term descriptive statistics##############
            tempSpAvgTS_ltMean[1, d] <- round(mean(inputDat, na.rm = TRUE), rd)
            tempSpAvgTS_ltSd[1, d] <- round(sd(inputDat, na.rm = TRUE), rd)
            ####################################################################
            
            #clear variables after iteration to avoid error
            rm(inputDat, mkTest)
          } #end of data sources
          
          #save results for each index within the adequate lists
          MKtrend_pval[[reg_names[r]]][[perSet[v]]][[ind1]] <- MK_pValue
          MKtrend_slope[[reg_names[r]]][[perSet[v]]][[ind1]] <- MK_slope
          MKtrend_interc[[reg_names[r]]][[perSet[v]]][[ind1]] <- MK_intercept
          
          allSim.allInd_spAvgTS[[reg_names[r]]][[perSet[v]]][[ind1]] <- allSimTS #INDIVIDUAL TIME SERIES FOR EACH DATA SOURCE
          
          ###########save the temporal long-term descriptive statistics#########
          tempLT_mean[[reg_names[r]]][[perSet[v]]][[ind1]] <- tempSpAvgTS_ltMean
          tempLT_sd[[reg_names[r]]][[perSet[v]]][[ind1]] <- tempSpAvgTS_ltSd
          ######################################################################
          
          
          ############################multimodel ensemble data aggregation######
          #now, get the mean, median, min and max annual time series from the pull of all data sources
              #empty matrices to save aggregated time series for further analyses
          allSimTS_stats <- array(data = NA, dim = c(length(tperiods[[v]]), 5)) ##CAUTION: HARD CODING (4 statistical metrics are considered here: mean, median, min, max)
          allSimTS_stats[, 1] <- tperiods[[v]] #CAUTION: HARD CODING (column 1 is for the years vector!)
              #loop over the annual time series
          for (y in 1:length(tperiods[[v]])) {
            #get yearly data vector
            data0 <- allSim.allInd_spAvgTS[[reg_names[r]]][[perSet[v]]][[ind1]][y, ]
            #calculate yearly stats with respect to the following conditions
            if (all(is.na(data0)) == TRUE) {
              allSimTS_stats[y, 2] <- NA
              allSimTS_stats[y, 3] <- NA
              allSimTS_stats[y, 4] <- NA
              allSimTS_stats[y, 5] <- NA
            } else {
              allSimTS_stats[y, 2] <- round(mean(data0, na.rm = TRUE), rd)
              allSimTS_stats[y, 3] <- round(median(data0, na.rm = TRUE), rd)
              #IMPORTANT: these are used to plot a shaded background depicting data range from all data sources evaluated
              allSimTS_stats[y, 4] <- round(min(data0, na.rm = TRUE), rd)
              allSimTS_stats[y, 5] <- round(max(data0, na.rm = TRUE), rd)
            }
            #clear variable after each iteration
            rm(data0)
          } #end of years
          #store descriptive stats annual time series in the dedicated list
          colnames(allSimTS_stats) <- EnsTS_cNames; allSimTS_stats <- as.data.frame(allSimTS_stats)
          allSim.allInd_spAvgTS_stat[[reg_names[r]]][[perSet[v]]][[ind1]] <- allSimTS_stats #THESE ARE ANNUAL TIME SERIES FROM MULTIMODELS PERSPECTIVE!!!
          ######################################################################
          
        } #end of drought indices sub-indices
      } #end of condition1 (drought indices)
      
      #condition2: remaining indices (general and heatwave ones)
      else {
        ind1 <- ind0 #index label of interest
        
        #empty matrix to save data retrieved for further analyses
        allSimTS <- array(data = NA, dim = c(length(tperiods[[v]]), length(climDatSource_list)))
        
        #empty matrices to save pvalue, slope and intercept outputs
        MK_pValue  <- array(data = NA, dim = c(1, length(climDatSource_list))) 
        MK_slope <- MK_pValue; MK_intercept <- MK_pValue
        
        ###################################
        #empty matrices to save long-term (lt) temporal mean and standard deviation (Sd) statistics for each index and data sources
        tempSpAvgTS_ltMean <- array(data = NA, dim = c(1, length(climDatSource_list)))
        tempSpAvgTS_ltSd <- tempSpAvgTS_ltMean
        ###################################
        
        #loop over all data sources
        for (d in 1:length(reorderedPos_climDatSources)) {
          newPos <- reorderedPos_climDatSources[d]
          
          inputDat <- annETSCI_spAvgTS[[newPos]][[reg_names[r]]][[ind1]][[perSet[v]]]
          #call and run the trend analysis function
          mkTest <- general_test_mk_prewhiten(inputDat, acLag, mkAlpha) ##CAUTION: HARD CODING below based on knowledge of the test's outputs
          MK_pValue[1, d] <- mkTest[[1]]
          MK_slope[1, d] <- round(mkTest[[2]], mkDigit)
          MK_intercept[1, d] <- round(mkTest[[3]], mkDigit)
          #now get data for further analyses
          allSimTS[, d] <- inputDat
          
          #########get temporal long-term descriptive statistics################
          tempSpAvgTS_ltMean[1, d] <- round(mean(inputDat, na.rm = TRUE), rd)
          tempSpAvgTS_ltSd[1, d] <- round(sd(inputDat, na.rm = TRUE), rd)
          ######################################################################
          
          #clear variables after iteration to avoid error
          rm(inputDat, mkTest)
        } #end of data sources
        
        #save results for each index within the adequate lists
        MKtrend_pval[[reg_names[r]]][[perSet[v]]][[ind1]] <- MK_pValue
        MKtrend_slope[[reg_names[r]]][[perSet[v]]][[ind1]] <- MK_slope
        MKtrend_interc[[reg_names[r]]][[perSet[v]]][[ind1]] <- MK_intercept
        
        allSim.allInd_spAvgTS[[reg_names[r]]][[perSet[v]]][[ind1]] <- allSimTS #INDIVIDUAL TIME SERIES FOR EACH DATA SOURCE
        
        ###########save the temporal long-term descriptive statistics###########
        tempLT_mean[[reg_names[r]]][[perSet[v]]][[ind1]] <- tempSpAvgTS_ltMean
        tempLT_sd[[reg_names[r]]][[perSet[v]]][[ind1]] <- tempSpAvgTS_ltSd
        ########################################################################
        
        
        ############################multimodel ensemble data aggregation########
        #now, get the mean, median, min and max annual time series from the pull of all data sources (for multimodel ensemble data aggregation)
            #empty matrices to save aggregated time series for further analyses
        allSimTS_stats <- array(data = NA, dim = c(length(tperiods[[v]]), 5)) #4 statistical metrics are considered here (mean, median, min, max)
        allSimTS_stats[, 1] <- tperiods[[v]]
            #loop over the annual time series
        for (y in 1:length(tperiods[[v]])) {
          #get yearly data vector
          data0 <- allSim.allInd_spAvgTS[[reg_names[r]]][[perSet[v]]][[ind1]][y, ]
          #calculate yearly stats with respect to the following conditions
          if (all(is.na(data0)) == TRUE) {
            allSimTS_stats[y, 2] <- NA
            allSimTS_stats[y, 3] <- NA
            allSimTS_stats[y, 4] <- NA
            allSimTS_stats[y, 5] <- NA
          } else {
            allSimTS_stats[y, 2] <- round(mean(data0, na.rm = TRUE), rd)
            allSimTS_stats[y, 3] <- round(median(data0, na.rm = TRUE), rd)
            #IMPORTANT: these are used to plot a shaded background depicting data range from all data sources evaluated
            allSimTS_stats[y, 4] <- round(min(data0, na.rm = TRUE), rd)  
            allSimTS_stats[y, 5] <- round(max(data0, na.rm = TRUE), rd)
          }
          #clear variable after each iteration
          rm(data0)
        } #end of years
        #store descriptive stats annual time series in the dedicated list
        colnames(allSimTS_stats) <- EnsTS_cNames; allSimTS_stats <- as.data.frame(allSimTS_stats)
        allSim.allInd_spAvgTS_stat[[reg_names[r]]][[perSet[v]]][[ind1]] <- allSimTS_stats #THESE ARE ANNUAL TIME SERIES FROM MULTIMODELS PERSPECTIVE!!!
        ########################################################################
        
      } #end of condition2
      
    } #end of selected periods
  } #end of regions
} #end of indices
tend4 <- Sys.time()

    #get the list of all indices assessed
indFullSet_label <- names(MKtrend_slope[[reg_names[r]]][[perSet[v]]])

    #concatenate indices labels and units for trend analysis & temporal long-term descriptive statistics SUMMARY TABLES!!!
indUnit_trend <- paste0(indFullSet_label, " ", indFullSet_unit)

    #now, loop over target regions and selected periods, to aggregate and save results for all indices and data sources
MKtrend_allResults <- list() ; tempLTstats_allResults <- list()

for (r in 1:length(reg_names)) {
  for (v in 1:length(perSet)) {
    
    #for trend analysis (CAUTION: HARD CODING based on trend output results!)
    outputFile1 <- paste0(outpath2, slash, trendFolder, slash, reg_names[r], slash, perSet[v], slash, MKtrendVar[1],"_", reg_names[r], "_", perSet[v], ".csv")
    outputFile2 <- paste0(outpath2, slash, trendFolder, slash, reg_names[r], slash, perSet[v], slash, MKtrendVar[2], "_", reg_names[r], "_", perSet[v], ".csv")
    outputFile3 <- paste0(outpath2, slash, trendFolder, slash, reg_names[r], slash, perSet[v], slash, MKtrendVar[3], "_", reg_names[r], "_", perSet[v], ".csv")
    
    mat_pvalue <- array(data = NA, dim = c(length(indFullSet_label), length(climDatSource_names)))
    mat_slope <- mat_pvalue; mat_intercept <- mat_pvalue
    
    #for temporal long-term descriptive statistics (CAUTION: HARD CODING based on defined descriptive statistics of interest!)
    outputFile4 <- paste0(outpath2, slash, tempLTstats, slash, reg_names[r], slash, perSet[v], slash, tempLTVar[1],"_", reg_names[r], "_", perSet[v], ".csv")
    outputFile5 <- paste0(outpath2, slash, tempLTstats, slash, reg_names[r], slash, perSet[v], slash, tempLTVar[2], "_", reg_names[r], "_", perSet[v], ".csv")
    
    mat_tpLTmean <- array(data = NA, dim = c(length(indFullSet_label), length(climDatSource_names))) ; mat_tpLTsd <- mat_tpLTmean
    
    
    for (z in 1:length(indFullSet_label)) {
      
      mat_pvalue[z, ] <- MKtrend_pval[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]]
      mat_slope[z, ] <- MKtrend_slope[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]]
      mat_intercept[z, ] <- MKtrend_interc[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]]
      
      mat_tpLTmean[z, ] <- tempLT_mean[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]]
      mat_tpLTsd[z, ] <- tempLT_sd[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]]
      
      
    } #end of full set of indices
    
    #save outputs as tables
    colnames(mat_pvalue) <- climDatSource_names[reorderedPos_climDatSources] 
    colnames(mat_slope) <- climDatSource_names[reorderedPos_climDatSources]
    colnames(mat_intercept) <- climDatSource_names[reorderedPos_climDatSources]
    rownames(mat_pvalue) <- indUnit_trend ; rownames(mat_slope) <- indUnit_trend; rownames(mat_intercept) <- indUnit_trend
    
    colnames(mat_tpLTmean) <- climDatSource_names[reorderedPos_climDatSources]
    colnames(mat_tpLTsd) <- climDatSource_names[reorderedPos_climDatSources]
    rownames(mat_tpLTmean) <- indUnit_trend; rownames(mat_tpLTsd) <- indUnit_trend
    
    write.csv(mat_pvalue, outputFile1)
    write.csv(mat_slope, outputFile2)
    write.csv(mat_intercept, outputFile3)
    
    write.csv(mat_tpLTmean, outputFile4)
    write.csv(mat_tpLTsd, outputFile5)
    
    MKtrend_allResults[[reg_names[r]]][[perSet[v]]][[MKtrendVar[1]]] <- as.data.frame(mat_pvalue)
    MKtrend_allResults[[reg_names[r]]][[perSet[v]]][[MKtrendVar[2]]] <- as.data.frame(mat_slope)
    MKtrend_allResults[[reg_names[r]]][[perSet[v]]][[MKtrendVar[3]]] <- as.data.frame(mat_intercept)
    
    tempLTstats_allResults[[reg_names[r]]][[perSet[v]]][[tempLTVar[1]]] <- as.data.frame(mat_tpLTmean)
    tempLTstats_allResults[[reg_names[r]]][[perSet[v]]][[tempLTVar[2]]] <- as.data.frame(mat_tpLTsd)
    
  } #end of selected periods
} #end of regions

#>>>>>>>>>>>>>>>>>>>Start of multimodel ensemble perspective>>>>>>>>>>>>>>>>>>>>
#run the trend analysis with multimodel ensemble annual means and medians if set as so!
if (multiModel == "YES") {
  #create a list to store results
  ensModels_MKtrend <- list()
  
  #loop over target region, time periods and indices to get input data and run trend analysis
  for (r in 1:length(reg_names)) {
    for (v in 1:length(perSet)) {
      
      #define output file
      outFile <- paste0(outpath2, slash, trendFolder, slash, reg_names[r], slash, perSet[v], slash, "multimodel_MKtrend_", reg_names[r], "_", perSet[v], ".csv")
      
      #organize results in an array for all indices: col1:3 for ensMeans; col4:6 for ensMedians
      ens_MKresults <- array(data = NA, dim = c(length(indFullSet_label), 6)) #col1:3=ensMeans; col2=ensMedians
      cNames <- c("pvalue_Emean", "slope_Emean", "intercept_Emean", "pvalue_Emedian", "slope_Emedian", "intercept_Emedian") #CAUTION: HARD CODING (based on analysis needs!)
      
      for (z in 1:length(indFullSet_label)) {
        
        #get input data
        ensAnnMeans <- allSim.allInd_spAvgTS_stat[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]]$mean
        ensAnnMedians <- allSim.allInd_spAvgTS_stat[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]]$median
        
        #call designed function that run trend analysis for single time series
        ensModels_MKtrend[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]][["annMeans"]] <- 
          trendAnalysis_func(ensAnnMeans, acLag, mkAlpha, mkDigit)
        ensModels_MKtrend[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]][["annMedians"]] <- 
          trendAnalysis_func(ensAnnMedians, acLag, mkAlpha, mkDigit)
        
        #compile results (CAUTION: HARD CODED based on analysis needs!)
        ens_MKresults[z, 1] <- ensModels_MKtrend[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]][["annMeans"]]$pvalue
        ens_MKresults[z, 2] <- ensModels_MKtrend[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]][["annMeans"]]$slope
        ens_MKresults[z, 3] <- ensModels_MKtrend[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]][["annMeans"]]$intercept
        ens_MKresults[z, 4] <- ensModels_MKtrend[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]][["annMedians"]]$pvalue
        ens_MKresults[z, 5] <- ensModels_MKtrend[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]][["annMedians"]]$slope
        ens_MKresults[z, 6] <- ensModels_MKtrend[[reg_names[r]]][[perSet[v]]][[indFullSet_label[z]]][["annMedians"]]$intercept
        
      } #end of indices
      
      #save outputs as tables
      colnames(ens_MKresults) <- cNames
      rownames(ens_MKresults) <- indUnit_trend
      write.csv(ens_MKresults, outFile)
      
    } #end of selected periods
  } #end of regions
} #end of condition
#>>>>>>>>>>>>>>>>>>>>>End of multimodel ensemble perspective>>>>>>>>>>>>>>>>>>>>


################################################################################[EXPECTED TIME (FOR CURRENT TARGET REGIONS): ABOUT +5MIN!]
#step2: produce temporal evolution plots of the indices annual time series for the full time period, then add the historical and future periods trend lines

    #years arrangement for upper and lower annual bounds polygon that show data range from all data sources assessed
polyX   <- c(seq(min(fullPer), max(fullPer)), rev(seq(min(fullPer), max(fullPer))))

#you may assess running time for below loop
tstart5 <- Sys.time()

#loop over indices and data sources annual time series and trend analysis results table for pvalue, and get data to produce plots
for (r in 1:length(reg_names)) {
  reg <- reg_names[r]
  
  for (z in 1:length(indFullSet_label)) {
    
    indexLabel <- indFullSet_label[z]
    
    #plot settings
    graphTitle <- paste0(indFullSet_def[z], " ", indFullSet_unit[z]) 
    graphFile <- paste0(outpath2, slash, tsPlotFolder, slash, reg, slash, reg, "_tempEvol_", indexLabel, ".png")
    
    #get the multimodel ensemble data range
    ens_minTS <- allSim.allInd_spAvgTS_stat[[reg]][[fullPer_label]][[indexLabel]]$min
    ens_maxTS <- allSim.allInd_spAvgTS_stat[[reg]][[fullPer_label]][[indexLabel]]$max
    polyEns <- cbind.data.frame(polyX, c(ens_minTS, rev(ens_maxTS)))
    colnames(polyEns) <- c("polyX", "polyY")
    
    #save subplot for each data source in a list for further aggregation
    Tfig <- list()
    
    for (d in 1:length(climDatSource_names)) {
      
      #NOTE: SINCE THE ORDER OF THE DATA SOURCES WERE REVISED for the temporal-scale analysis 
            #TO GATHER low, medium and high emissions, used below variable to further get the adequate data source name
      newPos <- reorderedPos_climDatSources[d] 
      
      #get MK trend results (CAUTION: HARD CODING based on trend output results!)
      MKpval_hist <- MKtrend_allResults[[reg]][[fullHist_label]][[MKtrendVar[1]]][z, d]
      MKpval_fut <- MKtrend_allResults[[reg]][[fullFut_label]][[MKtrendVar[1]]][z, d]
      
      MKslope_hist <- MKtrend_allResults[[reg]][[fullHist_label]][[MKtrendVar[2]]][z, d]
      MKslope_fut <- MKtrend_allResults[[reg]][[fullFut_label]][[MKtrendVar[2]]][z, d]
      
      MKinterc_hist <- MKtrend_allResults[[reg]][[fullHist_label]][[MKtrendVar[3]]][z, d]
      MKinterc_fut <- MKtrend_allResults[[reg]][[fullFut_label]][[MKtrendVar[3]]][z, d]
      
      #derive trend lines data for both periods
      trendHdat <- cbind.data.frame(fullHistPer, MKslope_hist * 1:length(fullHistPer) + MKinterc_hist)
      trendFdat <- cbind.data.frame(fullFutPer, MKslope_fut * 1:length(fullFutPer) + MKinterc_fut)
      colnames(trendHdat) <- c("years", "val"); colnames(trendFdat) <- c("years", "val")
      
      #set conditions for trend line
      if (is.na(MKpval_hist) == TRUE) {
        HtrendLine <- "dashed"
      } else if (MKpval_hist <= mkAlpha) {
        HtrendLine <- "solid"
      } else if (MKpval_hist > mkAlpha) {
        HtrendLine <- "dashed"
      }
      
      if (is.na(MKpval_fut) == TRUE) {
        FtrendLine <- "dashed"
      } else if (MKpval_fut <= mkAlpha) {
        FtrendLine <- "solid"
      } else if (MKpval_fut > mkAlpha) {
        FtrendLine <- "dashed"
      }
      
      #get annual time series for the full time coverage
      indexTS <- cbind.data.frame(fullPer, allSim.allInd_spAvgTS[[reg]][[fullPer_label]][[indexLabel]][, d])
      colnames(indexTS) <- c("years", "val")
      
      #plot
      subTitle <- climDatSource_names[newPos]
      
          #set condition for plotting for drought indices versus the others
              #CAUTION: preliminary analysis of data showed from some data sources, 3 levels of sign (-1, 0, 1); thus, a constraint is added below to limit the level to negative vs positive values
      if (is.element(indexLabel, spi_set_label) == TRUE || is.element(indexLabel, spei_set_label) == TRUE) {
        fig <- ggplot(indexTS, aes(x = years, y = val)) + 
          geom_bar(stat="identity", na.rm=TRUE, aes(colour = factor(sign(val), levels = c(-1, 1)), 
                                                    fill = factor(sign(val), levels = c(-1, 1)))) +
          scale_colour_manual(values = c("darkred", "blue")) + scale_fill_manual(values = alpha(c(polyCol, polyCol), 0.05)) + 
          theme_bw() + 
          geom_polygon(data = polyEns, mapping = aes(x = polyX, y = polyY), fill = polyCol, alpha = colAlpha) +
          geom_hline(yintercept = 0) + 
          geom_line(data = trendHdat, mapping = aes(x = years, y = val), colour = trendcolH, linetype = HtrendLine, linewidth = trendLw) +
          geom_line(data = trendFdat, mapping = aes(x = years, y = val), colour = trendcolF, linetype = FtrendLine, linewidth = trendLw) +
          labs(title = subTitle) + 
          theme(legend.position = "None", 
                plot.title = element_text(size = txtSize1, hjust = figT_hjust),
                axis.title.x=element_blank(), axis.title.y=element_blank(),
                axis.text.x = element_text(size = txtSize2), axis.text.y = element_text(size = txtSize2))
      } #end of drought indices condition
      else {
        fig <- ggplot() + theme_bw() + 
          geom_polygon(data = polyEns, mapping = aes(x = polyX, y = polyY), fill = polyCol, alpha = colAlpha) + 
          geom_line(data = indexTS, mapping = aes(x = years, y = val), colour = indTScol, linewidth = TSplotLw) +
          geom_line(data = trendHdat, mapping = aes(x = years, y = val), colour = trendcolH, linetype = HtrendLine, linewidth = trendLw) +
          geom_line(data = trendFdat, mapping = aes(x = years, y = val), colour = trendcolF, linetype = FtrendLine, linewidth = trendLw) +
          labs(title = subTitle) + 
          theme(plot.title = element_text(size = txtSize1, hjust = figT_hjust),
                axis.title.x=element_blank(), axis.title.y=element_blank(),
                axis.text.x = element_text(size = txtSize2), axis.text.y = element_text(size = txtSize2))
      } #end of other indices condition
      
          #store subplot
      Tfig[[d]] <- fig
      
      #clear variable after iteration to avoid errors
      rm(MKpval_hist, MKpval_fut, MKslope_hist, MKslope_fut, MKinterc_hist, MKinterc_fut, trendHdat, trendFdat, indexTS, fig)
      
    } #end of data sources
    
    #save graph for each index using function created at this end
    aggTSsubplots_func(Tfig)
    
    #clear variable after iteration to avoid errors
    rm(ens_minTS, ens_maxTS, polyEns)
    
  } #end of indices full list
} #end of regions
tend5 <- Sys.time()


################################################################################[EXPECTED TIME (FOR CURRENT TARGET REGIONS): LESS THAN 5MIN!]
#step3: quantify projected future annual changes for define future periods relative to baseline conditions,
        #then generate boxplots and tables to summarize results

#define steps for near and far future, time period common length vector, and shared time period label (CAUTION: HARD CODED based on number of future periods assessed; here, 2!)
idPosF1 <- seq(1, (length(climDatSource_names)*2), 2) ; idPosF2 <- seq(2, (length(climDatSource_names)*2), 2)
sharedFutLabel <- array(data <- NA, dim = c(1, 2*length(climDatSource_names)))
for (l in 1:length(climDatSource_names)) {
  newPos <- reorderedPos_climDatSources[l]
  sharedFutLabel[1, idPosF1[l]] <- paste0(climDatSource_names[newPos], "_F1")
  sharedFutLabel[1, idPosF2[l]] <- paste0(climDatSource_names[newPos], "_F2")
}

#derive indices units for annual change
indUnit_deltas <- array(data = NA, dim = c(length(indFullSet_label), 1))
indDef_deltas <- array(data = NA, dim = c(length(indFullSet_label), 1))
for (z in 1:length(indFullSet_label)) {
  indexLabel <- indFullSet_label[z] #index evaluated
  #set conditions for indices units related to change assessment
  if (is.element(indexLabel, spei_set_label) == TRUE || is.element(indexLabel, spi_set_label) == TRUE) {
    indUnit_deltas[z] <- paste0(indFullSet_label[z], " ", droughtInd_units)
    indDef_deltas[z] <- paste0(indFullSet_def[z], " ", droughtInd_units)
  } else if (is.element(indexLabel, relativeDelInd) == TRUE) {
    indUnit_deltas[z] <- paste0(indFullSet_label[z], " ", relativeDelInd_unit)
    indDef_deltas[z] <- paste0(indFullSet_def[z], " ", relativeDelInd_unit)
  } else {
    indUnit_deltas[z] <- paste0(indFullSet_label[z], " ", indFullSet_unit[z])
    indDef_deltas[z] <- paste0(indFullSet_def[z], " ", indFullSet_unit[z])
  }
} #end of indices

#create lists to save output results
annETSCI_deltas_temporalRes <- list() ; deltas_temporalStat <- list() ; deltasBoxplots <- list()

#you may assess running time for below loop
tstart6 <- Sys.time()

#loop over regions and indice, get input data data from all data sources and quantify annualchanges
for (r in 1:length(reg_names)) {
  reg <- reg_names[r]
  
  #initialize an empty matrix to store long-term temporal stats calculated from annual change distributions over near and far future periods
  deltas_temporalMean <- array(data = NA, dim = c(length(indFullSet_label), length(climDatSource_names)*2))
  deltas_temporalSd <- deltas_temporalMean ; deltas_temporalMin <- deltas_temporalMean ; deltas_temporalMax <- deltas_temporalMean
  
  #set output tables filenames (CAUTION: HARD CODING based on descriptive statistics of interest!)
  output1 <- paste0(outpath2, slash, futDeltFolder, slash, reg_names[r], slash, annDeltasMetrics[1],"_", reg, ".csv")
  output2 <- paste0(outpath2, slash, futDeltFolder, slash, reg_names[r], slash, annDeltasMetrics[2],"_", reg, ".csv")
  output3 <- paste0(outpath2, slash, futDeltFolder, slash, reg_names[r], slash, annDeltasMetrics[3],"_", reg, ".csv")
  output4 <- paste0(outpath2, slash, futDeltFolder, slash, reg_names[r], slash, annDeltasMetrics[4],"_", reg, ".csv")

  for (z in 1:length(indFullSet_label)) {
    indexLabel <- indFullSet_label[z] #index evaluated
    #individual saving of indices boxplots
    ind_bpOutFile <- paste0(outpath2, slash, futDeltFolder, slash, reg_names[r], slash, bp_fname, indexLabel, "_", reg, ".png")
    
    #set conditions for values to be round regarding the index category
    if (is.element(indexLabel, count_indices) == TRUE) {
      rd <- rd1
    } else if (is.element(indexLabel, spei_set_label) == TRUE || is.element(indexLabel, spi_set_label) == TRUE) {
      rd <- rd3
    } else {
      rd <- rd2
    }

    #initialize empty matrix to store calculated changes over near and far future periods, and store all results from all data sources
    deltasMat_all <- array(data = NA, dim = c(nyears, length(climDatSource_names)*2))
    
    #setting for boxplots
    subplotTitle <- paste0(greeks("Delta"), " ", indUnit_deltas[z]) #title for each graph (index label & unit)
    ##subplotTitle <- paste0(greeks("Delta"), " ", indDef_deltas[z]) #title for each graph (index definition with label & unit)

    #loop over data sources to calculate changes for each year (NOTE HERE that the the data sources were already reordered into the called list according to the User settings!!!)
    for (d in 1:length(climDatSource_names)) {

      refDat <- allSim.allInd_spAvgTS[[reg]][[refLabel]][[indexLabel]][, d] #baseline data
      fut1Dat <- allSim.allInd_spAvgTS[[reg]][[fut1Label]][[indexLabel]][, d] #near future data
      fut2Dat <- allSim.allInd_spAvgTS[[reg]][[fut2Label]][[indexLabel]][, d] #far future data
      
      #initialize empty matrix to store calculate changes for each data source (CAUTION: HARD CODING based on number of future periods evaluated; here 2!)
      deltasMat <- array(data = NA, dim = c(nyears, 2)) #col1 is for near future; col2 for far future
      
      ##############################################
      #calculate the baseline temporal mean condition (single value)
      ref_TpMeanVal <- mean(refDat, na.rm=TRUE)
      ##############################################

      #set condition for exceedance indices, relative changes indices, and absolute changes indices
      if (is.element(indexLabel, exceedRate_indices) == TRUE) {
        deltasMat[, 1] <- fut1Dat
        deltasMat[, 2] <- fut2Dat
      } #end of condition1

      else if (is.element(indexLabel, relativeDelInd) == TRUE) {
        for (y in 1:nyears) {
          deltasMat[y, 1] <- ((fut1Dat[y] - ref_TpMeanVal) / ref_TpMeanVal)*100
          deltasMat[y, 2] <- ((fut2Dat[y] - ref_TpMeanVal) / ref_TpMeanVal)*100
        } #end of period length
      } #end of condition2

      else {
        for (y in 1:nyears) {
          deltasMat[y, 1] <- fut1Dat[y] - ref_TpMeanVal
          deltasMat[y, 2] <- fut2Dat[y] - ref_TpMeanVal
        } #end of period length
      } #end of condition3
      
      #call function which calculate the temporal annual changes long-term descriptive statistics
      outF1 <- generateTemporalStats_func(deltasMat[, 1], rd)
      outF2 <- generateTemporalStats_func(deltasMat[, 2], rd)
          #then aggregate results for all indices and data sources
      deltas_temporalMean[z, idPosF1[d]] <- outF1$mean
      deltas_temporalSd[z, idPosF1[d]] <- outF1$sd
      deltas_temporalMin[z, idPosF1[d]] <- outF1$min
      deltas_temporalMax[z, idPosF1[d]] <- outF1$max
      
      deltas_temporalMean[z, idPosF2[d]] <- outF2$mean
      deltas_temporalSd[z, idPosF2[d]] <- outF2$sd
      deltas_temporalMin[z, idPosF2[d]] <- outF2$min
      deltas_temporalMax[z, idPosF2[d]] <- outF2$max
      
      #compiling all results for further generate boxplot per index for all data sources
      deltasMat_all[, idPosF1[d]] <- deltasMat[, 1]
      deltasMat_all[, idPosF2[d]] <- deltasMat[, 2]

      #clear variables after iteration to avoid errors
      rm(refDat, fut1Dat, fut2Dat, outF1, outF2)

    } #end of data sources
    
    #compile and save results as dataframe to generate boxplots of annual change distribution per index but for all data sources
    deltasBoxplotData <- cbind.data.frame(c(1:nyears), deltasMat_all)
    colnames(deltasBoxplotData) <- c("nyears", sharedFutLabel)
    annETSCI_deltas_temporalRes[[reg]][[indexLabel]] <- deltasBoxplotData
    
    #generate graph for each index and save in a list for further aggregation based on the User needs!
        #1st, convert plotting dataframe into long format for multiple variable plotting
    longDf <- melt(annETSCI_deltas_temporalRes[[reg]][[indexLabel]], id="nyears")
        #add a new column to distinguish near and far future
    longDf$period <- str_sub(longDf$variable, start = -2) #CAUTION: HARD CODING based on preliminary analysis!
    
        #generate and saved graph for each index and target region as individual files
    fig1 <- ggplot(data=longDf, aes(x=variable, y=value, color=period)) + theme_bw() + 
      geom_boxplot() + labs(title = subplotTitle, x = bp_xlabel) +
      scale_color_manual(values = deltasBpCol, name='Period', labels=c(fut1Label, fut2Label)) + 
      scale_x_discrete(breaks = sharedFutLabel[idPosF1], labels = bp_xaxis) +
      theme(legend.position = "bottom", legend.text = element_text(size=bp_txt1), legend.title = element_text(size = bp_txt1, face = "bold"),
            axis.title.x = element_text(size=bp_txt2), axis.title.y=element_blank(), 
            axis.text.x = element_text(size=bp_txt2, angle = axis_ang), axis.text.y = element_text(size=bp_txt1), 
            plot.title = element_text(face="bold", size=bp_txt3, hjust = bp_hjust))
    ggsave(ind_bpOutFile, fig1, scale = ind_bpScale, dpi = ppi)
    
        #generate and saved graph for each index and target region in a list for further aggregation based on User's needs
    fig <- ggplot(data=longDf, aes(x=variable, y=value, color=period)) + theme_bw() + 
      geom_boxplot() + labs(title = subplotTitle, x = bp_xlabel) +
      scale_color_manual(values = deltasBpCol, name='Period', labels=c(fut1Label, fut2Label)) + 
      scale_x_discrete(breaks = sharedFutLabel[idPosF1], labels = bp_xaxis) +
      theme(legend.position = "none", legend.text = element_text(size=bp_txt1), legend.title = element_text(size = bp_txt1, face = "bold"),
            axis.title.x = element_text(size=bp_txt2), axis.title.y=element_blank(), 
            axis.text.x = element_text(size=bp_txt2, angle = axis_ang), axis.text.y = element_text(size=bp_txt1), 
            plot.title = element_text(face="bold", size=bp_txt3, hjust = bp_hjust))
        #save plot in a list for each index and target region
    deltasBoxplots[[reg]][[indexLabel]] <- fig
    
    #clear variables after iteration to avoid errors
    rm(deltasBoxplotData, longDf, fig1, fig)

  } #end of indices full list
  
  #compile and save results as summary tables for all indices and data sources
  colnames(deltas_temporalMean) <- sharedFutLabel ; colnames(deltas_temporalSd) <- sharedFutLabel
  colnames(deltas_temporalMin) <- sharedFutLabel ; colnames(deltas_temporalMax) <- sharedFutLabel
  rownames(deltas_temporalMean) <- indUnit_deltas ; rownames(deltas_temporalSd) <- indUnit_deltas
  rownames(deltas_temporalMin) <- indUnit_deltas ; rownames(deltas_temporalMax) <- indUnit_deltas
  
  write.csv(deltas_temporalMean, output1)
  write.csv(deltas_temporalSd, output2)
  write.csv(deltas_temporalMin, output3)
  write.csv(deltas_temporalMax, output4)
  
  #save in a list, results for each target region from all indices and data sources (CAUTION: HARD CODING based on descriptive statistics of interest!)
  deltas_temporalStat[[reg]][[temporalStats[1]]] <- as.data.frame(deltas_temporalMean)
  deltas_temporalStat[[reg]][[temporalStats[2]]] <- as.data.frame(deltas_temporalSd)
  deltas_temporalStat[[reg]][[temporalStats[3]]] <- as.data.frame(deltas_temporalMin)
  deltas_temporalStat[[reg]][[temporalStats[4]]] <- as.data.frame(deltas_temporalMax)
  
} #end of regions
tend6 <- Sys.time()


#################start of indices boxplots aggregation##########################
#loop over target regions to aggregate indices given groups defined by User 
    #CAUTION: any required modifications should be done within the function called! 
bp_outFolder <- "aggDeltasBoxplots"

for (r in 1:length(reg_names)) {
  reg <- reg_names[r] #target region
  reg_figList <- deltasBoxplots[[reg]] #list of figure for target region
  dir.create(file.path(outpath2, slash, futDeltFolder, slash, reg, bp_outFolder)) #create a subfolder within each region to store aggregated graphs
  
  #call and run function which aggregate individual indices boxplots
  getTempDeltasAggBoxplots_func(reg_figList, reg)
  
  #clear variable
  rm(reg_figList)
  
} #end of regions

###################end of indices boxplots aggregation##########################


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>THE END!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #print time taken to run different sections
tend1 - tstart1 #crop of netcdf data to the extend of the target regions
tend2 - tstart2 #retrieval of gridded and spatially averaged annual time series, and calculation of long-term spatial climatologies and future changes 
tend3 - tstart3 #production of long-term future changes spatial maps
tend - tstart #production of long-term climatologies spatial maps
tend4 - tstart4 #run of trend analysis and aggregation of long-term climatologies: temporal scale analysis
tend5 - tstart5 #production of temporal evolution plots with trend line
tend6 - tstart6 #assessment of temporal future changes






