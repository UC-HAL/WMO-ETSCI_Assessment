# R version 4.3.1 (2023-06-16 ucrt) | October 2023

# FOR LOCAL MACHINE (e.g., physical computer; virtual desktop)

# This script enable to calculate gridded ET-SCI indices from a set of netCDF files using the Climpact R software.
    # reference: https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md

# It calls a function which is a conversion of the Climpact original script "climpact.ncdf.wrapper.r" 
    # See https://github.com/ARCCSS-extremes/climpact
    # Thus, many code lines from this original script are used herein to fit with the created function.

######################################################################################################################################################
# THE USER SHOULD OVERVIEW THE "ReadMe_Climpact&ETSCI_setting-calculation.txt" FILE FOR DETAILS ABOUT CLIMPACT SOFTWARE AND THE ET-SCI LIST OF INDICES

# THE USER SHOULD ADJUST THE PATHS WHERE NECESSARY AS WELL AS ANY INPUT INFORMATION WITH RESPECT TO THE RESEARCH STUDY UNDERTAKEN
######################################################################################################################################################


# THE USER INPUT ARGUMENTS INCLUDE THE FOLLOWING:
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Define and set working directory and add the folder and files necessary to the script execution
user_wk <- "C:/Users/alida.thiombiano/OneDrive - University of Calgary/WMO-ETSCI_Assessment_R-Workflow/PART1_ETSCIcalculation/ETSCIfromGriddedTS" #current User's physical computer

setwd(user_wk)

# Loading (and if necessary installing) the R packages required by Climpact
    # This required the contents of the 'climpact-master' folder created when downloading Climpact software
    # 1st, make sure to do this beforehand:
        # download Climpact ZIP folder by opening 'Code' green box at https://github.com/ARCCSS-extremes/climpact
        # extract the 'climpact-master' folder to a directory
    # then, add the following extra-steps to copy the software files to the User's working directory
ClimpactFolder <- "C:/Users/alida.thiombiano/OneDrive - University of Calgary/WMO-ETSCI_Assessment_R-Workflow/PART1_ETSCIcalculation/ClimpactSoftwareFiles/climpact-master"
climpactFiles <- list.files(ClimpactFolder)
file.copy(file.path(ClimpactFolder, climpactFiles), user_wk, overwrite = FALSE, recursive = TRUE)

source('server/climpact.master.installer.r')

###########################NOTE: NO NEED TO DO THIS ALWAYS (NECESSARY FOR 1ST TIME INSTALLATION OF THE PACKAGE)#############################
    # Installing the 'climdex.pcic.ncdf' package for the 1st time on local machine
        # This package is a modified version of 'climdex.pcic' to enable gridded indices calculation from netCDF file.
        # There is a special procedure to follow to install it;
        # and this enables the calculation of all the ET-SCI list of indices (https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md#appendixa).
install.packages("./server/pcic_packages/climdex.pcic.ncdf.climpact.tar.gz",repos=NULL,type="source")
############################################################################################################################################
    

# Loading the function to be called for the ET-SCI calculation
source("climpact.ncdf.wrapper_userFunc.R")

# Loading R packages
library(PCICt) #required to load the 'climdex.pcic.ncdf' package
library(climdex.pcic.ncdf) #modified version of 'climdex.pcic' package to enable gridded indices calculation from netCDF file
library(parallelly) #package used to detect available cores on the computer used

# Specify Climpact software directory (full pathname).
    # Leave as NULL because extra-steps were taken above to copy the software files to the User's working directory.
software_folder <- NULL

# Set the path to the ncfiles directory
ncfiles_dir <- paste0(user_wk, "/User_ncfiles/")

# Set the path to calculation outputs
user_outpath <- paste0(user_wk, "/GriddedETSCI")

# Set the number of cores for parallel processing (use FALSE for a single core, or provide the number)
ncores <- availableCores() - 6 #if running other tasks on physical computer, the User may not want to use all available cores

# From the netcdf metadata, get the variable names for the three climate datasets required to calculate indices with Climpact
    # Check these with App like Panoply (https://www.giss.nasa.gov/tools/panoply/)
var1 <- "pr"
var2 <- "tmax"
var3 <- "tmin"

# Define reference period for percentile-based indices calculation
baseline <- c(1985, 2014)

# Provide a list of the ET-SCI of interest (use NULL to calculate all indices, or provide a specific list (e.g., c("hw","tnn")))
indList <- NULL
##indList <- c("txx", "tnn", "prcptot") #this may be also useful to run a test

# The ET-SCI percentile based-threshold indices are by default calculated from the input dataset.
    # if an external thresholds file (calculated from another dataset) is to be used, specify and provide it (e.g., 'thresholds.test.1991-1997.nc') with the directory path.
    # the climpact original 'climpact.ncdf.thresholds.wrapper.r' script may also be modified and used for external thresholds calculation from netCDF data
threshFile <- NULL #no external thresholds file

# User's information
user_institution <- "University of Calgary"
user_contact <- "Alida Thiombiano (alida.thiombiano@ucalgary.ca)"

# Number of data values to process at once. If receiving "Error: rows.per.slice >= 1 is not TRUE", try increasing this to 20. You might have a large grid.
    # 10 is the default setting inside the original script;
nval=20
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


##ET-SCI's CALCULATION SECTION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# get the list of the input netCDF files without their directory path
ncf_List <- list.files(path=ncfiles_dir, pattern="*.nc", all.files=FALSE, full.names=FALSE)

# iterate through the list
for (f in 1:length(ncf_List)) {
  
  # filename with path to directory
  ncfile <- paste0(ncfiles_dir, ncf_List[f]) 
  
  # create individual sub-directories to store results for each ncfile
  t <- 3 #removing the last 3 characters of each ncfilename ('.nc')
  subFolder <- substr(ncf_List[f], start = 1, stop = (nchar(ncf_List[f]) - t))
  dir.create(file.path(user_outpath, subFolder))
  outpath <- paste0(user_outpath, "/", subFolder,"/")
  
  # Climpact's naming convention
  outFilename <- paste0("var_", "tres_", ncf_List[f]) 
  
  # call the created function (i.e., Climpact original script modified)
  climpact.ncdf.wrapper_userFunc(ncfile, var1, var2, var3, outpath, outFilename, 
                              user_institution, user_contact, baseline, ncores,
                              indList, threshFile, software_folder, nval)
  
  # clear variables after iteration
  rm(ncfile, subFolder, outpath, outFilename)
  
}


# Clean-up the climpact files copied to the working directory
unlink(file.path(user_wk, climpactFiles), recursive = TRUE)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#THE END!!!





