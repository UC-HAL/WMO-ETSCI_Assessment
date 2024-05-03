# **WMO-ETSCI_Assessment_R-Workflow**

This repository proposes a workflow written in the R programming language, to help with a straightforward calculation and statistical analysis of hydroclimate indicators developed by the World Meteorological Organization (WMO)'s Expert Team on Sector-Specific Climate Indices (ET-SCI). These indicators are hereafter referred to as the ET-SCI. 

The workflow has two components, the second being executed upon completion of the first.


## **First component of the workflow**

The workflow starts with the calculation of the ET-SCI. 

The files to review and run this first component are stored in the folder **"PART1_ETSCIcalculation"**, and below is a general description of their contents:

- The sub-folder **"ClimpactSoftwareFiles"** contains folders, files, and codes necessary to install the Climpact software. These files are imported from the [ARCCSS-extremes / climpact](https://github.com/ARCCSS-extremes/climpact/) github repository, and incorporated into our workflow to enable a straightforward calculation of the ET-SCI, from netCDF or text files, using the Climpact software developed by the WMO's ET-SCI. 

- The sub-folder **"ETSCIfromGriddedTS"** contains folders, files, and codes necessary to run and store the calculation of the indices ***from netCDF files***. 
    - Two netCDF sample files are provided within the **"User_ncfiles"** folder, for the User to run a test.
    - The calculated indices are stored within the existing **"GriddedETSCI"** folder, created beforehand.
    - The User may run the calculation of the ET-SCI on local machine using the *"ETSCIfromGriddedTS_UserMainScript_LocalMachine.R"*, or on Linux server environment using the *"ETSCIfromGriddedTS_UserMainScript_LinuxServer.R"*. Note that the function *"climpact.ncdf.wrapper_userFunc"* is called within both scripts. An example of Linux Bash Shell script (*"SLURM_subjob_ETSCIfromGriddedTS.sh"*) is also provided for the Linux option.

- The sub-folder **"ETSCIfromSingleTS"** contains folders, files, and codes necessary to run and store the calculation of the indices ***from text files***. 
    - Two regionally averaged time series are provided as a sample within the **"inputSingleTS_txt"** folder, for the User to run a test.
    - The *"metadata_batchSingleTS.txt"* text file contains the required metadata associated to these sample files.
    - The User may run the calculation on local machine using the *"ETSCIfromSingleTS_UserMainScript_LocalMachine.R"*, and following the instruction to launch the calculation on the Climpact online portal, and retrieve the results where indicated. The User can also run the calculation on a Linux server environment using the Linux Bash Shell script *"SLURM_subjob_ETSCIfromSingleTS.sh"* provided as an example.

- The User should review the ***"ReadMe_Climpact&ETSCI_setting-calculation.txt"*** document for details about the Climpact software, its installation, the formatting it required for the climate input data, and the outputs formats of the indices it calculates. The images referred to in this file (*"Capture_climpact_RwebInterface.png"*; *"ETSCIfromSingleTS_Option2_ErrorMessage1.jpg"*; *"ETSCIfromSingleTS_Option2_ErrorMessage2.jpg"*) are stored in the **"PART1_ETSCIcalculation"** folder.

- The User should also review the ***"InstallingClimpact&Rpackages_LinuxTerminal.txt"*** document for the installation of Climpact and all other required R packages on their Linux server environment, before running any scripts for the calculation of the ET-SCI there.

- The *"climpact2-user-guide_Alexander and Herold_May2016.pdf"* within the **"PART1_ETSCIcalculation"** folder, is a resourceful document, similar to the [online resource](https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md), offering guidance for Climpact software application, and providing background and details for the ET-SCI.  


## **Second component of the workflow**

Following the calculation of the ET-SCI **ONLY** ***from netCDF files***, the User may use the contents of this second component of the workflow, which are stored in the **"PART2_ETSCIprocessing"** folder, to run **spatial and temporal statistical analyses of the ET-SCI at the annual scale**.

Below is a general description of the contents of the **"PART2_ETSCIprocessing"** folder:

- The existing sub-folders **"spatialAnalysis"** and **"temporalAnalysis"** are created beforehand to store respectively, the outputs of spatial and temporal analyses.

- The sub-folder **"codes"** contains all the scripts developed for the statistical processing of the gridded annual time series of the ET-SCI. The *"GridETSCI_SpatialTemporalAnalysis_annScale_mainScript.R"* constitutes the main script within which the functions from both the **"codes"** and **"UCHAL_MKtrend_Rfunctions"** folders are called. At the begining of each script and function, the developer has provided what each scripting does. Note also that the *"exploreIndices_ncfiles_Rscript.R"* script may be open and run for exploratory overview of the ET-SCI data.

    - The User should review the main script file *"GridETSCI_SpatialTemporalAnalysis_annScale_mainScript.R"* to have a global understanding of what it does, to adapt pathways to appropriate directories, and to change User-customized choices when and where necessary.
    - The User should also note that there are some hard coding within the following files, based on the developer knowledge of the study area and datasets used to implement the present workflow: 
        - script: "GridETSCI_SpatialTemporalAnalysis_annScale_mainScript.R"
        - script: "exploreIndices_ncfiles_Rscript.R"
        - function: "aggTSsubplots_func.R"
        - function: "getIndicesDataRange_forSpatialMaps_func.R"
        - function: "getTempDeltasAggBoxplots_func.R"
        - function: "trendAnalysis_func.R"

- The sub-folder **"targetRegion_shp"** is used to store the shapefile of the target region called within the main script *"GridETSCI_SpatialTemporalAnalysis_annScale_mainScript.R"*.
    - The shapefile of the St. Mary and Milk Rivers watershed is provided within this folder for the User to run a test in alignment with the sample of datasets provided for the ET-SCI's calculation.




