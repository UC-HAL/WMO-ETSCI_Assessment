###October 2023



#In this document, the User may learn about:

#(1) the Climpact R software and the indices it calculates;
#(2) how to install and run Climpact on both local machine and Linux server; 
#(3) the format of the climate input data required by Climpact;
#(4) the output formats of indices calculated by Climpact.



######################################################################(1)############################################################################
#As described on the following websites (https://climpact-sci.org/; https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md; https://climpact-sci.org/assets/climpact2-user-guide.pdf), 
#"Climpact is written in the R programming language and makes use of numerous third party R packages" to calculate indices of daily climate extremes (https://climpact-sci.org/indices/) to support adaptation and risk management. 
#the Github repository of Climpact may be explored at https://github.com/ARCCSS-extremes/climpact

#"Climpact uses the R packages 'climdex.pcic' and 'climdex.pcic.ncdf' as its core for calculating indices", and these packages and Climpact are developed by the Pacific Climate Impacts Consortium (https://www.pacificclimate.org/).

#Climpact can read data for a single site (e.g. a weather station; spatially average time series) in the form of a text file, or gridded data (e.g. from a climate model) in the form of netCDF files.

#The indices calculated by Climpact (https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md#appendixa; https://climpact-sci.org/indices/) 
#are known as the Expert Team on Sector-specific Climate Indices (ET-SCI) 
#and include the Expert Team on Climate Change Detection and Indices (ETCCDI) which are describe at http://etccdi.pacificclimate.org/list_27_indices.shtml.
	#among the ET-SCI, there are the so-called percentile threshold indices
		#the User may learn more about their calculation procedure at https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md#appendixc
	#regarding the heatwave and coldwave indices, their definition, aspects and calculation are detailed at https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md#appendixd

#Climpact's indices are robust standardized statistics derived from daily tmin., daily tmax., and daily precip., and enabling to
#detect, understand and monitor changes in mean, moderate and extreme climatic conditions with recurrence times of a year or shorter (e.g., Zhang et al. 2011; http://dx.doi.org/https://doi.org/10.1002/wcc.147), 
#evaluate climate models (e.g., Sillmann et al. 2013a,b; http://dx.doi.org/10.1002/jgrd.50203; http://dx.doi.org/10.1002/jgrd.50188; Braun et al. 2021; https://doi.org/10.1525/elementa.2021.00011) 
#generate meaningful climate information for decision-makers (e.g., Rajulapati et al. 2022; https://doi.org/10.1016/j.uclim.2022.101097).



######################################################################(2)############################################################################
#The Climpact original scripts are all listed on the Github repository of Climpact at https://github.com/ARCCSS-extremes/climpact, and are considered at this point.


############################Calculating indices from gridded data in netCDF format
#First, review procedure and requirements at the following link:
	#https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md#calculatenetcdf

#On a Linux server environment, 
	#follow the steps describe in "InstallingClimpact&Rpackages_LinuxTerminal.txt" to install Climpact.
	#then review and use (modify if and where necessary) the "climpact.ncdf.wrapper_userFunc.R", "ETSCIfromGriddedTS_UserMainScript_LinuxServer.R", and "SLURM_subjob_ETSCIfromGriddedTS.sh" scripts 
		#to run Climpact and calculate indices from netCDF input files.
		#The User may notice that the "climpact.ncdf.wrapper_userFunc.R" script is a conversion into a function of the Climpact original script "climpact.ncdf.wrapper.r"

#On a local machine, 
	#review and use (modify if and where necessary) the "climpact.ncdf.wrapper_userFunc.R" and "ETSCIfromGriddedTS_UserMainScript_LocalMachine.R" scripts, 
		#to install Climpact required packages and run the indices calculation from netCDF input files.

	
############################Calculating indices from station-type time series (data in '.txt' format)
#First, review procedure and requirements at the following links:
	#single dataset: https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md#calculatesingle
	#batch of dataset: https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md#calculatebatch 


#On a local machine, 
	#option1: use the online portal at https://ccrc-extremes.shinyapps.io/climpact/
		#when opened, fill-in the required information, launch the calculation, and click "OK" in the pop-up box.
		#when calculation is done, take note of the indications (shown on blue box) to retrieve outputs, and close the portal (in this case the open web-page).
	
	#option2: review then use (modify if and where necessary) the 'ETSCIfromSingleTS_UserMainScript_LocalMachine.R' script.
		#The User should keep in mind that this also open the Climpact online portal (in this case a R session web-interface),
		#and ask for the User to proceed as indicated with option 1; but, the User may also with option2, indicate the number of cores to use.
			#check the "Capture_climpact_RwebInterface.png" file for an overview.
		#CAUTION: with option2, errors warning message occur (check "ETSCIfromSingleTS_Option2_ErrorMessage1.jpg" and "ETSCIfromSingleTS_Option2_ErrorMessage2.jpg"), 
			#and some indices in the ET-SCI list are ultimately not calculated.
			#THUS, THE USER SHOULD GO WITH OPTION1.

#On a Linux server environment,
	#review and use (modify if and where necessary) the "SLURM_subjob_ETSCIfromSingleTS.sh" script to run it on a Linux server environment;
		#notice that the 'climpact.batch.stations.r' script is passed in the form of command line arguments and does not require modification.
	#CAUTION: as with the above referred option2, similar errors warning message occur.

#CONCLUDING ADVICE: THE USER SHOULD GO WITH OPTION1 TO CALCULATE INDICES FROM STATION-TYPE TIME SERIES!!!



######################################################################(3)############################################################################
#To calculate the ET-SCI list of indices, the User need daily time series of total precipitation (PR), minimum temperature (TN) and maximum temperature (TX)
	#these datasets may be single time series representation of a station, grid, region, ...
	#or gridded datasets in netCDF format from climate models outputs

#The User should prepare these input climate datasets to meet the formats required by Climpact:
	#see details at https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md#appendixb
	#see Climpact's sample data at https://github.com/ARCCSS-extremes/climpact/tree/master/www/sample_data or in the same directory from your local climpact folder

#Few things to keep in mind:
	#The Climate Data Operator (CDO) software (https://code.mpimet.mpg.de/projects/cdo/embedded/index.html) may be used for pre-processing netCDF files
		#e.g., merging or spliting files;  adjusting metadata like variables units; 
		#following Climpact's input and output files naming convention (i.e. "var_timeresolution_model_scenario_run_starttime-endtime.nc")
			#the first 2 parts ("var_tres_") will be automatically replaced by the index calculated (e.g., var=cdd) 
			#and its temporal resolution (tres=ANN and/or tres=MON)
			#e.g., climate simulation: "var_tres_BCC-CSM2-MR_hist-ssp126_NA_1950-2100.nc"; observation: "var_tres_Hybrid-BCABSKUS_Obs_NA_NA_1950_2019.nc"

	#for the station-type files, NO NEED TO FOLLOW THE CLIMPACT INPUT AND OUTPUT FILES NAMING CONVENTION



######################################################################(4)############################################################################
#The User should get familiar with the format of the indices time series as calculated by Climpact to understand the design of any further analyses.

#from station-type text files,
	#the indices annual and/or monthly time series are '.csv' files and located in the "indices" sub-directory of the output folder created by the system.
		#more details are available at https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md#outputstation

#from netCDF files,
	#gridded annual and/or monthly time series of the indices are produced in netCDF formats;
	#Tools like Panoply App (https://www.giss.nasa.gov/tools/panoply/) can be used to visualize the embedded data and metadata. 
	#example: https://github.com/ARCCSS-extremes/climpact/blob/master/www/user_guide/Climpact_user_guide.md#outputgridded



######################################################################################################################################################################
#With respect to the User's interest, the shell and R scripts provided for the Climpact setting and the ET-SCI calculation may be used for any target region globally.
######################################################################################################################################################################

