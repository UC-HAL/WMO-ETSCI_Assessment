# R version 4.3.1 (2023-06-16 ucrt) | October - December 2023

#Explore the netCDF format with respect to the indices categories

    #checking the full time period
nt <- length(tperiods[[1]])

    #selecting indices from "allETSCI_files" variable
cdd_pos <- 1 #cdd index
hw_pos <- 13 #heatwave group of indices
spei_pos <- 35 #drought group of time scales

    #exploring index "cdd" netcdf
ncf_cdd <- paste0(ETSCI_dir, slash, climDatSource_list[1], slash, allETSCI_files[cdd_pos])
ncin_cdd <- nc_open(ncf_cdd)
ncdat_cdd <- ncvar_get(ncin_cdd, "cdd") #see original data's dimension (here lon x lat x time)
TS_cdd <- ncdat_cdd[150, 120, ] #overview time series for one grid; the last value is a NA because the length exceeds "nt"
tvec_cdd <- as.numeric(format(nc.get.time.series(ncin_cdd),'%Y')) #checking years, there is indeed an extra row
        #to calculate the descriptive statistics over different long-term periods, this information will be considered to retrieve the proper rows.

    #exploring index "hw" netcdf
ncf_hw <- paste0(ETSCI_dir, slash, climDatSource_list[1], slash, allETSCI_files[hw_pos])
ncin_hw <- nc_open(ncf_hw)
ncdat_hw <- ncvar_get(ncin_hw, "hwn_ehf") #calling one index among the 20 heatwaves indices embedded within "hw_ANN"
TS_hw <- ncdat_hw[150, 120, ]
tvec_hw <- as.numeric(format(nc.get.time.series(ncin_hw),'%Y'))
        #conclusion: beside calling the heatwave indices individually, the format is the same as of the general indices

    #exploring index "spei" netcdf
ncf_spei <- paste0(ETSCI_dir, slash, climDatSource_list[1], slash, allETSCI_files[spei_pos])
ncin_spei <- nc_open(ncf_spei)
ncdat_spei <- ncvar_get(ncin_spei, "spei") #this is a 4D data (here lon x lat x time x scale; in the 4th dimension there 3 scales for respectively the 3-months, 6-months, 12-months moving windows)
TS_spei <- ncdat_spei[150, 120, , 1] #get the 3-months time series
tvec_spei <- as.numeric(format(nc.get.time.series(ncin_spei),'%Y'))
        #checking the last two years, there is indeed an extra year to do not count for during the descriptive statistics calculation
tail(TS_spei, 24)
tail(tvec_spei, 24)
        #conclusion: for each time series retrieved, the last 12 values should be omitted.

    #transform one netcdf into a rasterbrick object to explore the format
cdd_raster <- brick(x=ncf_cdd, varname="cdd", ncdf=TRUE)
cdd_raster;
        ###CAUTION: notice that there is a now a change in the data dimension:
                    #(lat x lon x time) from the raster versus (lon x lat x time) from the original netcdf.

    #clear all variables to save memory
rm(nt, cdd_pos, hw_pos, spei_pos, ncf_cdd, ncin_cdd, ncdat_cdd, TS_cdd, tvec_cdd,
   ncf_hw, ncin_hw, ncdat_hw, TS_hw, tvec_hw,
   ncf_spei, ncin_spei, ncdat_spei, TS_spei, tvec_spei, cdd_raster)


