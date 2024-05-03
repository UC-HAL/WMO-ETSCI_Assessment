#this function takes a rasterstack, 
#process it, 
#and generate spatial maps for several data sources at once (as subplot).

generateSpatialMaps_func <- function(rastObject, colorbarLim) {
  
  #define spatial grid from the raster, and transform the output object in a dataframe
  ras_df_ggplt <- as.data.frame(as(rastObject, "SpatialPixelsDataFrame")) 
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
  ##ggsave(figName, fig, scale = fscale, dpi=ppi) #replace "scale = fscale" by "width = ..., height = ..., units = ".." for a better control on plot size
  ggsave(figName, fig, width = myWidth, height = myHeight, units = myUnits, dpi=ppi)

  
} #end of function

