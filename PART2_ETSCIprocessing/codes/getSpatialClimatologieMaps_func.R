#create spatial maps illustrating long-term climatologies (mean, median, min, max) and save them in the subfolder created at this end

getSpatialClimatologieMaps_func <- function() {
  
  #creating the a subfolder to save maps illustrating long-term climatologies calculated 
  Results_normals <- "Climatologies"
  dir.create(file.path(outpath1, slash, Results_normals))
  
  #creating within these folders, subfolders to separate outputs for each target region
  #and within these folders, create subfolders to separate results given each descriptive statistics
  for (f in 1:length(reg_names)) {
    dir.create(file.path(outpath1, slash, Results_normals, slash, reg_names[f]))
    for (ff in 1:length(legTitle)) {
      dir.create(file.path(outpath1, slash, Results_normals, slash, reg_names[f], slash, legTitle[ff]))
    }
  }
  
  #now, create maps of ET-SCI long-term annual climatologies (here, mean, median, min and max are the selected metrics)
  
  #loop over indices, regions, future periods assessed, and set of change metrics,
  #and pull gridded data for all data sources, then generate maps with subplots (each subplot illustrating one data source)
  
  
  for (iii in 1:length(annInd_list)) {
    indexLabel <- annInd_list[iii] #index label
    colPal <- annInd_list_colPal[iii] #color palette
    colDir <- annInd_list_coldir[iii] #color scale orientation
    
    for (h in 1:length(reg_names)) {
      xlon_breaks <- latlon_breaks[[h]][["lon"]]
      ylat_breaks <- latlon_breaks[[h]][["lat"]]
      regShp <- reg_polygon[polygon_pos[h], ]
      
      for (pp in 1:length(tper_strings)) {
        
        for (m in 1:length(statMetrics)) {
          
          if (indexLabel == annInd_list[pos_cdd]) {
            if (statMetrics[m] == "min") {
              colDir <- colOrder[2] #CAUTION: HARD CODING!!!
            } else {
              colDir <- colOrder[1] #CAUTION: HARD CODING!!!
            }
          }
          
          #folder where to store maps
          outFolder <- paste0(outpath1, slash, Results_normals, slash, reg_names[h], slash, legTitle[m])
          
          ######################
          #set condition1: drought indices (CAUTION: HARD CODING BASED ON SELECTED TIME SCALES FOR SPI AND SPEI!!!)
          if(is.element(indexLabel, drought_ind) == TRUE) {
            
            #set legend title
            legTitle_full <- paste0(legTitle[m], "\n", droughtInd_units)
            
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
              figName <- paste0(outFolder, slash, drought_subset2[ww], "_", tper_strings[pp], ".png")
              figTitle <- indDef[ww]
              colorbarLim <- indRange_spClim[[indexLabel]][[drought_subset2[ww]]][[statMetrics[m]]] #common limit for colorbar
              
              #and now over data sources (CAUTION: HERE, THE ORIGINAL ORDER "length(annETSCI_deltas_gridStat)" OF DATA SOURCES ARE REVISED TO GATHER low, medium and high emissions)
              for (dd in 1:length(reorderedPos_climDatSources)) {
                newPos <- reorderedPos_climDatSources[dd]
                #matrix
                inDat <- annETSCI_gridStat[[newPos]][[reg_names[h]]][[indexLabel]][[drought_subset2[ww]]][[tper_strings[pp]]][[statMetrics[m]]]
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
              
              #generating maps for each sub-index
                  #define spatial grid from the raster, and transform the output object in a dataframe
              ras_df_ggplt <- as.data.frame(as(output_ras, "SpatialPixelsDataFrame")) 
                  #check here about label mismatch (e.g., ras_df_ggplt$BCC.CSM2.MR_hist.ssp126)
              
                  #use reshape2 package to melt the dataframe regarding the coordinates (x,y)
              ras_df_ggplt_melt <- melt(ras_df_ggplt, id.vars = c('x','y'))
              
                      #CAUTION: check raster layer labels (subplots variables names) string that it represents meaningful labels
                      #fix the variables' names until it represents the desired label
              ras_df_ggplt_melt$variable <- gsub("\\.", "-", ras_df_ggplt_melt$variable) #'dot' replacement can be tricky; so you may check this link: https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot
              
                  #generate the ggplot (using "facet_wrap" to make a 1d long ribbon of panels: https://ggplot2.tidyverse.org/reference/facet_wrap.html)
              fig <- ggplot(data = ras_df_ggplt_melt) +
                geom_tile(aes(x = x, y = y, fill = value)) +
                geom_sf(data = regShp, fill = NA, color = shp_col, linewidth = shp_lwd) + 
                scale_fill_viridis_b(option = colPal, direction = colDir, limits = c(colorbarLim[1], colorbarLim[2])) +
                facet_wrap(.~factor(variable, levels = climDatSource_names[reorderedPos_climDatSources])) +
                theme_bw() +
                labs(title = figTitle, x = regLon_label, y= regLat_label, fill =  legTitle_full) +
                theme(legend.title = element_text(size = xyLabel_sz), legend.position = legPos,
                      axis.title.x = element_text(size = xyLabel_sz),
                      axis.title.y = element_text(size = xyLabel_sz),
                      axis.text.x = element_text(size = xyTicks_sz, angle = axis_ang),
                      axis.text.y = element_text(size = xyTicks_sz),
                      plot.title = element_text(face = txt_style, size = figTitle_sz, hjust = figT_hjust),
                      strip.text = element_text(size = subTitle_sz, face = txt_style)) +
                scale_x_continuous(limits = c(min(xlon_breaks), max(xlon_breaks)), breaks = xlon_breaks, labels = as.character(xlon_breaks)) +
                scale_y_continuous(limits = c(min(ylat_breaks), max(ylat_breaks)), breaks = ylat_breaks, labels = as.character(ylat_breaks))
              
                  #save the ggplot
              #ggsave(figName, fig, scale = fscale, dpi=ppi) #replace "scale = fscale" by "width = ..., height = ..., units = ".." for a better control on plot size
              ggsave(figName, fig, width = myWidth, height = myHeight, units = myUnits, dpi=ppi)
              
              
              #clear variables to avoid errors
              rm(output_ras)
              
            } #end of drought sub-index
          } #end of condition1
          
          ######################
          #set condition2: other indices (i.e., (is.element(indexLabel, general_ind) == TRUE || is.element(indexLabel, hwCw_ind) == TRUE))
          else {
            #set legend title, figure saving name, figure title
            legTitle_full <- paste0(legTitle[m], "\n", allInd_unit[iii])
            figName <- paste0(outFolder, slash, indexLabel, "_", tper_strings[pp], ".png")
            figTitle <- allInd_def[iii]
            colorbarLim <- indRange_spClim[[indexLabel]][[statMetrics[m]]] #common limit for colorbar
            
            #and now over data sources
            for (dd in 1:length(reorderedPos_climDatSources)) {
              newPos <- reorderedPos_climDatSources[dd]
              #matrix
              inDat <- annETSCI_gridStat[[newPos]][[reg_names[h]]][[indexLabel]][[tper_strings[pp]]][[statMetrics[m]]]
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
            
            #generating maps for each sub-index
                #define spatial grid from the raster, and transform the output object in a dataframe
            ras_df_ggplt <- as.data.frame(as(output_ras, "SpatialPixelsDataFrame")) 
                #check here about label mismatch (e.g., ras_df_ggplt$BCC.CSM2.MR_hist.ssp126)
            
                #use reshape2 package to melt the dataframe regarding the coordinates (x,y)
            ras_df_ggplt_melt <- melt(ras_df_ggplt, id.vars = c('x','y'))
            
                    #CAUTION: check raster layer labels (subplots variables names) string that it represents meaningful labels
                    #fix the variables' names until it represents the desired label
            ras_df_ggplt_melt$variable <- gsub("\\.", "-", ras_df_ggplt_melt$variable) #'dot' replacement can be tricky; so you may check this link: https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot
            
                #generate the ggplot (using "facet_wrap" to make a 1d long ribbon of panels: https://ggplot2.tidyverse.org/reference/facet_wrap.html)
            fig <- ggplot(data = ras_df_ggplt_melt) +
              geom_tile(aes(x = x, y = y, fill = value)) +
              geom_sf(data = regShp, fill = NA, color = shp_col, linewidth = shp_lwd) + 
              scale_fill_viridis_b(option = colPal, direction = colDir, limits = c(colorbarLim[1], colorbarLim[2])) +
              facet_wrap(.~factor(variable, levels = climDatSource_names[reorderedPos_climDatSources])) +
              theme_bw() +
              labs(title = figTitle, x = regLon_label, y= regLat_label, fill =  legTitle_full) +
              theme(legend.title = element_text(size = xyLabel_sz), legend.position = legPos,
                    axis.title.x = element_text(size = xyLabel_sz),
                    axis.title.y = element_text(size = xyLabel_sz),
                    axis.text.x = element_text(size = xyTicks_sz, angle = axis_ang),
                    axis.text.y = element_text(size = xyTicks_sz),
                    plot.title = element_text(face = txt_style, size = figTitle_sz, hjust = figT_hjust),
                    strip.text = element_text(size = subTitle_sz, face = txt_style)) +
              scale_x_continuous(limits = c(min(xlon_breaks), max(xlon_breaks)), breaks = xlon_breaks, labels = as.character(xlon_breaks)) +
              scale_y_continuous(limits = c(min(ylat_breaks), max(ylat_breaks)), breaks = ylat_breaks, labels = as.character(ylat_breaks))
            
                #save the ggplot
            #ggsave(figName, fig, scale = fscale, dpi=ppi) #replace "scale = fscale" by "width = ..., height = ..., units = ".." for a better control on plot size
            ggsave(figName, fig, width = myWidth, height = myHeight, units = myUnits, dpi=ppi)
            
            #clear variables to avoid errors
            rm(output_ras)
            
          } #end of condition2
          ######################
          
        } #end of descriptive stats
      } #end of assessed periods
    } #end of regions
  } #end of indices
  
} #end of function

