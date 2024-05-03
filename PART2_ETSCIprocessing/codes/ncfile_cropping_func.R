#This function takes a netcdf file, 
#transform it into a rasterbrick object, 
#crop it to the extent of the shapefile of a target region, 
#mask the cropped raster, 
#and transform it back into an array for numerical calculation.

ncfile_cropping_func <- function (shp, nc_pathFilename, indexVarname, var_level, tscale_level) {
  
  #read/transform the netcdf file in a multi-layer raster object to shorten multi-layer file processing time (https://search.r-project.org/CRAN/refmans/raster/html/brick.html);
      #if more than one variable is embedded, specify the one to read, otherwise, a random one is automatically chosen
      #similarly, set condition with respect to 3D vs 4D ncfiles to specify the level to use for the 4th dimension (default is 1)
        #this specifically applies to the drought indices
  if ((is.element(indexVarname, c("spei", "spi")) == TRUE)) {
    rastData <- brick(x=nc_pathFilename, varname=indexVarname, lvar=var_level, level=tscale_level, ncdf=TRUE)
  } 
  else {
    rastData <- brick(x=nc_pathFilename, varname=indexVarname, ncdf=TRUE)
  }
  
  #assure both the raster object and the region shapefile are in the same coordinate systems
  raster_crs <- st_crs(rastData)
  shp2 <- st_transform(shp, raster_crs)
  
  #crop the masked raster to the extent of the shapefile to reduce data size and to save computer memory
  cropRaster <- crop(rastData, extent(shp2))
  
  #save the cropped raster for spatial maps
  arrayData_crop <- as.array(cropRaster)
  
  #mask the cropped raster object (transform area outside the polygon to NA val) to simplify operation for retrieving a specific area (https://rpubs.com/ricardo_ochoa/416711)
  maskRaster <- mask(cropRaster, shp2)
  
  #transform the cropped and masked raster back in array for numerical calculation
  arrayData_cropMask <- as.array(maskRaster) 
  
  #return the output of interest
  return(list("crop"=arrayData_crop, "cropMask"=arrayData_cropMask))
  
} #end of function

