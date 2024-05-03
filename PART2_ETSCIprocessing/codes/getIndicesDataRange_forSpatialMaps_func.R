
#the following function get indices' values range common to all data sources, time periods & target regions 
    #for spatial maps colorbar limits setting

getIndicesDataRange_forSpatialMaps_func <- function() {
  
  #variables called: "annETSCI_spOrigDataRange" for case 1; "annETSCI_spDeltasDataRange" for case2
           #case1: for long-term climatologies' spatial maps
  indRange_spClim <- list()
           #case2: for long-term projected changes' spatial maps
  indRange_spDel <- list()
  
  #define columns' position for data sources with respect to target region
  dsPos_reg <- vector("list", nreg) # for the number of target regions!
      #CAUTION: HARD CODED (herein, two target regions were considered!)
  dsPos_reg[[1]] <- seq(1, length(climDatSource_list), 1) #all data sources for target region n1
  dsPos_reg[[2]] <- seq(length(climDatSource_list)+1, length(climDatSource_list)*2, 1) #all data sources for target region n2
  
  #loop over set of indices
  for (ii in 1:length(annInd_list)) {
    indVar2 <- annInd_list[[ii]]
    
    ################################################################################  
    #for drought indices sub-group
    if(is.element(indVar2, drought_ind) == TRUE) {
      
      #get sub-goup of indices names with respect to SPEI or SPI
      if (indVar2 == "spei3" || indVar2 == "spi3") {
        drought_subset <- names(annETSCI_spOrigDataRange[[ss]][[reg]][[indVar2]])
      } #end of 3-months time scale indices (here 4 indices)
      else if (indVar2 == "spei12" || indVar2 == "spi12") {
        drought_subset <- names(annETSCI_spOrigDataRange[[ss]][[reg]][[indVar2]])
      } #end of 12-months time scale indices (here 1 index)
      
      for (w in 1:length(drought_subset)) {
        indVar4 <- drought_subset[w]
        
        #create empty matrices to retrieve data for each index: case1 (single time series for the full time period)
        outTable_origMean <- array(data = NA, dim = c(length(tperiods[[nfull]]), length(climDatSource_list)*nreg)) #columns are all data sources for all target regions
        outTable_origMedian <- outTable_origMean ; outTable_origMin <- outTable_origMean ; outTable_origMax <- outTable_origMean
        
        #create empty matrices to retrieve data for each index: case2 (single min&max values )
        outTable_DeltasMean <- array(data = NA, dim = c(4, length(climDatSource_list)*nreg)) #here, row1&2 (row3&4)for near&far future min (max) values |CAUTION: HARD CODING (based on analysis design!)
        outTable_DeltasMedian <- outTable_DeltasMean ; outTable_DeltasMin <- outTable_DeltasMean ; outTable_DeltasMax <- outTable_DeltasMean
        
        #loop over target regions and data sources, and compile data needed
        for (rr in 1:length(reg_names)){
          reg <- reg_names[rr] ; dsPos <- dsPos_reg[[rr]]
          for (ss in 1:length(climDatSource_list)) {
            
            outTable_origMean[, dsPos[ss]] <- annETSCI_spOrigDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fullPer_label]]$meanTS
            outTable_origMedian[, dsPos[ss]] <- annETSCI_spOrigDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fullPer_label]]$medianTS
            outTable_origMin[, dsPos[ss]] <- annETSCI_spOrigDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fullPer_label]]$minTS
            outTable_origMax[, dsPos[ss]] <- annETSCI_spOrigDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fullPer_label]]$maxTS
            
            outTable_DeltasMean[1, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]]$meanRange[1]
            outTable_DeltasMedian[1, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]]$medianRange[1]
            outTable_DeltasMin[1, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]]$minRange[1]
            outTable_DeltasMax[1, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]]$maxRange[1]
            
            outTable_DeltasMean[2, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]]$meanRange[1]
            outTable_DeltasMedian[2, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]]$medianRange[1]
            outTable_DeltasMin[2, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]]$minRange[1]
            outTable_DeltasMax[2, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]]$maxRange[1]
            
            outTable_DeltasMean[3, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]]$meanRange[2]
            outTable_DeltasMedian[3, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]]$medianRange[2]
            outTable_DeltasMin[3, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]]$minRange[2]
            outTable_DeltasMax[3, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut1Label]]$maxRange[2]
            
            outTable_DeltasMean[4, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]]$meanRange[2]
            outTable_DeltasMedian[4, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]]$medianRange[2]
            outTable_DeltasMin[4, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]]$minRange[2]
            outTable_DeltasMax[4, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[indVar4]][[fut2Label]]$maxRange[2]
            
          } #end of data sources
        } #end of target regions
        
        #get c(min, max) vector for each index
        #case1
        indRange_spClim[[annInd_list[[ii]]]][[indVar4]][["mean"]] <- c(floor(min(outTable_origMean, na.rm = TRUE)), 
                                                               ceiling(max(outTable_origMean, na.rm = TRUE)))
        
        indRange_spClim[[annInd_list[[ii]]]][[indVar4]][["median"]] <- c(floor(min(outTable_origMedian, na.rm = TRUE)), 
                                                                 ceiling(max(outTable_origMedian, na.rm = TRUE)))
        
        indRange_spClim[[annInd_list[[ii]]]][[indVar4]][["min"]] <- c(floor(min(outTable_origMin, na.rm = TRUE)), 
                                                              ceiling(max(outTable_origMin, na.rm = TRUE)))
        
        indRange_spClim[[annInd_list[[ii]]]][[indVar4]][["max"]] <- c(floor(min(outTable_origMax, na.rm = TRUE)), 
                                                              ceiling(max(outTable_origMax, na.rm = TRUE)))
        #case2
        indRange_spDel[[annInd_list[[ii]]]][[indVar4]][["mean"]] <- c(floor(min(outTable_DeltasMean[1:2, ], na.rm = TRUE)), 
                                                              ceiling(max(outTable_DeltasMean[2:4, ], na.rm = TRUE)))
        
        indRange_spDel[[annInd_list[[ii]]]][[indVar4]][["median"]] <- c(floor(min(outTable_DeltasMedian[1:2, ], na.rm = TRUE)), 
                                                                ceiling(max(outTable_DeltasMedian[2:4, ], na.rm = TRUE)))
        
        indRange_spDel[[annInd_list[[ii]]]][[indVar4]][["min"]] <- c(floor(min(outTable_DeltasMin[1:2, ], na.rm = TRUE)), 
                                                             ceiling(max(outTable_DeltasMin[2:4, ], na.rm = TRUE)))
        
        indRange_spDel[[annInd_list[[ii]]]][[indVar4]][["max"]] <- c(floor(min(outTable_DeltasMax[1:2, ], na.rm = TRUE)), 
                                                             ceiling(max(outTable_DeltasMax[2:4, ], na.rm = TRUE)))
        
      } #end of drought indices subset
    } #end of drought indices
    
    ################################################################################  
    #for remaining indices
    else {
      #create empty matrices to retrieve data for each index: case1 (single time series for the full time period)
      outTable_origMean <- array(data = NA, dim = c(length(tperiods[[nfull]]), length(climDatSource_list)*nreg)) #columns are all data sources for all target regions (here, there 2)
      outTable_origMedian <- outTable_origMean ; outTable_origMin <- outTable_origMean ; outTable_origMax <- outTable_origMean
      
      #create empty matrices to retrieve data for each index: case2 (single min&max values )
      outTable_DeltasMean <- array(data = NA, dim = c(4, length(climDatSource_list)*nreg)) #here, row1&2 (row3&4)for near&far future min (max) values |CAUTION: HARD CODING (based on analysis design!)
      outTable_DeltasMedian <- outTable_DeltasMean ; outTable_DeltasMin <- outTable_DeltasMean ; outTable_DeltasMax <- outTable_DeltasMean
      
      #loop over target regions and data sources, and compile data needed
      for (rr in 1:length(reg_names)){
        reg <- reg_names[rr] ; dsPos <- dsPos_reg[[rr]]
        for (ss in 1:length(climDatSource_list)) {
          
          outTable_origMean[, dsPos[ss]] <- annETSCI_spOrigDataRange[[ss]][[reg]][[indVar2]][[fullPer_label]]$meanTS
          outTable_origMedian[, dsPos[ss]] <- annETSCI_spOrigDataRange[[ss]][[reg]][[indVar2]][[fullPer_label]]$medianTS
          outTable_origMin[, dsPos[ss]] <- annETSCI_spOrigDataRange[[ss]][[reg]][[indVar2]][[fullPer_label]]$minTS
          outTable_origMax[, dsPos[ss]] <- annETSCI_spOrigDataRange[[ss]][[reg]][[indVar2]][[fullPer_label]]$maxTS
          
          outTable_DeltasMean[1, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut1Label]]$meanRange[1]
          outTable_DeltasMedian[1, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut1Label]]$medianRange[1]
          outTable_DeltasMin[1, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut1Label]]$minRange[1]
          outTable_DeltasMax[1, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut1Label]]$maxRange[1]
          
          outTable_DeltasMean[2, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut2Label]]$meanRange[1]
          outTable_DeltasMedian[2, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut2Label]]$medianRange[1]
          outTable_DeltasMin[2, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut2Label]]$minRange[1]
          outTable_DeltasMax[2, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut2Label]]$maxRange[1]
          
          outTable_DeltasMean[3, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut1Label]]$meanRange[2]
          outTable_DeltasMedian[3, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut1Label]]$medianRange[2]
          outTable_DeltasMin[3, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut1Label]]$minRange[2]
          outTable_DeltasMax[3, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut1Label]]$maxRange[2]
          
          outTable_DeltasMean[4, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut2Label]]$meanRange[2]
          outTable_DeltasMedian[4, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut2Label]]$medianRange[2]
          outTable_DeltasMin[4, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut2Label]]$minRange[2]
          outTable_DeltasMax[4, dsPos[ss]] <- annETSCI_spDeltasDataRange[[ss]][[reg]][[indVar2]][[fut2Label]]$maxRange[2]
          
        } #end of data sources
      } ##end of target regions
      
      #get c(min, max) vector for each index
      #case1
      indRange_spClim[[annInd_list[[ii]]]][["mean"]] <- c(floor(min(outTable_origMean, na.rm = TRUE)), 
                                                  ceiling(max(outTable_origMean, na.rm = TRUE)))
      
      indRange_spClim[[annInd_list[[ii]]]][["median"]] <- c(floor(min(outTable_origMedian, na.rm = TRUE)), 
                                                    ceiling(max(outTable_origMedian, na.rm = TRUE)))
      
      indRange_spClim[[annInd_list[[ii]]]][["min"]] <- c(floor(min(outTable_origMin, na.rm = TRUE)), 
                                                 ceiling(max(outTable_origMin, na.rm = TRUE)))
      
      indRange_spClim[[annInd_list[[ii]]]][["max"]] <- c(floor(min(outTable_origMax, na.rm = TRUE)), 
                                                 ceiling(max(outTable_origMax, na.rm = TRUE)))
      #case2
      indRange_spDel[[annInd_list[[ii]]]][["mean"]] <- c(floor(min(outTable_DeltasMean[1:2, ], na.rm = TRUE)), 
                                                 ceiling(max(outTable_DeltasMean[2:4, ], na.rm = TRUE)))
      
      indRange_spDel[[annInd_list[[ii]]]][["median"]] <- c(floor(min(outTable_DeltasMedian[1:2, ], na.rm = TRUE)), 
                                                   ceiling(max(outTable_DeltasMedian[2:4, ], na.rm = TRUE)))
      
      indRange_spDel[[annInd_list[[ii]]]][["min"]] <- c(floor(min(outTable_DeltasMin[1:2, ], na.rm = TRUE)), 
                                                ceiling(max(outTable_DeltasMin[2:4, ], na.rm = TRUE)))
      
      indRange_spDel[[annInd_list[[ii]]]][["max"]] <- c(floor(min(outTable_DeltasMax[1:2, ], na.rm = TRUE)), 
                                                ceiling(max(outTable_DeltasMax[2:4, ], na.rm = TRUE)))
      
    } #end of the other indices
    
    
  } #end of loop over full list of indices
  
  #return outputs of interest
  return(list("indRange_spClim"=indRange_spClim, "indRange_spDel"=indRange_spDel))
  
  
} #end of function



