# R version 4.3.1 (2023-06-16 ucrt) | October 2023

# FOR LOCAL MACHINE (e.g., physical computer; virtual desktop)

# This script enable to calculate ET-SCI indices from single station-type time series using the Climpact R software.
    # reference: https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md

    # For a station-type single dataset, review procedure at https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md#calculatesingle
    # For a station-type batch of dataset, review procedure at https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md#calculatebatch

######################################################################################################################################################
# THE USER SHOULD OVERVIEW THE "ReadMe_Climpact&ETSCI_setting-calculation.txt" FILE FOR DETAILS ABOUT CLIMPACT SOFTWARE AND THE ET-SCI LIST OF INDICES

# THE USER SHOULD ADJUST THE PATHS WHERE NECESSARY AS WELL AS ANY INPUT INFORMATION WITH RESPECT TO THE RESEARCH STUDY UNDERTAKEN
######################################################################################################################################################


# THE USER INPUT ARGUMENTS INCLUDE THE FOLLOWING:
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Define and set working directory and add the folder and files necessary to the script execution
user_wk <- "C:/Users/alida.thiombiano/OneDrive - University of Calgary/WMO-ETSCI_Assessment_R-Workflow/PART1_ETSCIcalculation/ETSCIfromSingleTS" #current User's physical computer

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
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    

#ET-SCI's CALCULATION SECTION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Run the following two commands 
    # this opens the Climpact online portal (in this case a R session web-interface);
    # fill-in the required information, launch the calculation, and click "OK" in the pop-up box;
    # when calculation is completed on the interface, take note of the indications (shown on blue box) to retrieve outputs, then close the portal.

#loading 'shiny' R package for interactive web application
library(shiny)

#running the 'app.R'
runApp()
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



################################################################################
# after the interface is closed, come back in this script;

# then, move the outputs from the indicated default folder to the User's designated folder; 
    # note that retrieving the outputs is VERY IMPORTANT as the default folder is part of the files temporarily copied to the working directory

#then, clean-up the climpact files copied to the working directory
unlink(file.path(user_wk, climpactFiles), recursive = TRUE)
################################################################################


#THE END!!!





