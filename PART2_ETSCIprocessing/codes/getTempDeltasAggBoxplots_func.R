#this function aggregates individual indices boxplots given groups defined by User's need

#CAUTION: HARD CODING FOR GROUP SETTINGS AND SUBPLOTS AGGREGATION!!! (USER MAY MODIFY IT TO MEET OTHER NEEDS!)

getTempDeltasAggBoxplots_func <- function(fig_list, domain) {
  
  #define number of indices' groups of interest (here there are five groups evaluated!!!)
  study_indGroups <- vector("list", 5)
  
      #8 indices related to hot temperatures (hotTemp)
  study_indGroups[[1]][["indLabels"]] <- c("txx", "txn", "txm", "wsdi", "tx90p", "tx10p","su", "id")
  study_indGroups[[1]][["groupName"]] <- paste0(outpath2, slash, futDeltFolder, slash, domain, slash, 
                                                bp_outFolder, slash, bp_fname, "hotTemp_", domain, ".png")
  
      #8 indices related to cold temperatures (coldTemp)
  study_indGroups[[2]][["indLabels"]] <- c("tnn", "tnx", "tnm", "csdi", "tn10p", "tn90p","fd", "tr")
  study_indGroups[[2]][["groupName"]] <- paste0(outpath2, slash, futDeltFolder, slash, domain, slash, 
                                                bp_outFolder, slash, bp_fname, "coldTemp_", domain, ".png")
  
      #8 indices related to moderate temperatures and heatwaves (modTempHWehf)
  study_indGroups[[3]][["indLabels"]] <- c("tmm", "dtr", "gsl", "tmge5", "tmge10", "hwn_ehf", "hwd_ehf", "hwf_ehf")
  study_indGroups[[3]][["groupName"]] <- paste0(outpath2, slash, futDeltFolder, slash, domain, slash, 
                                                bp_outFolder, slash, bp_fname, "modTempHWehf_", domain, ".png")
  
      #8 indices related to general and extreme patterns in precipitation (genExtPrecip)
  study_indGroups[[4]][["indLabels"]] <- c("prcptot", "sdii", "cdd", "cwd", "r95ptot", "r99ptot","rx1day", "rx5day")
  study_indGroups[[4]][["groupName"]] <- paste0(outpath2, slash, futDeltFolder, slash, domain, slash, 
                                                bp_outFolder, slash, bp_fname, "genExtPrecip_", domain, ".png")
  
      #5 indices related towater balance (SPEI)
  study_indGroups[[5]][["indLabels"]] <- c("Feb_SPEI3", "May_SPEI3", "Aug_SPEI3",  "Nov_SPEI3",  "Sep_SPEI12")
  study_indGroups[[5]][["groupName"]] <- paste0(outpath2, slash, futDeltFolder, slash, domain, slash, 
                                                bp_outFolder, slash, bp_fname, "SPEI_", domain, ".png")
  
  #loop over the list, produce and save aggregated figures
  for (g in 1:length(study_indGroups)) {
    
    if (g == 5) {
      ggsave(study_indGroups[[g]][["groupName"]], ggarrange(fig_list[[study_indGroups[[g]][["indLabels"]][1]]],
                                                            fig_list[[study_indGroups[[g]][["indLabels"]][2]]],
                                                            fig_list[[study_indGroups[[g]][["indLabels"]][3]]],
                                                            fig_list[[study_indGroups[[g]][["indLabels"]][4]]],
                                                            fig_list[[study_indGroups[[g]][["indLabels"]][5]]],
                                                            nrow=3, ncol=2, common.legend = TRUE, legend="bottom"),
             width = 3800, height = 3400, units = myUnits, dpi=ppi)
    } #case1
    else {
      ggsave(study_indGroups[[g]][["groupName"]], ggarrange(fig_list[[study_indGroups[[g]][["indLabels"]][1]]],
                                                            fig_list[[study_indGroups[[g]][["indLabels"]][2]]],
                                                            fig_list[[study_indGroups[[g]][["indLabels"]][3]]],
                                                            fig_list[[study_indGroups[[g]][["indLabels"]][4]]],
                                                            fig_list[[study_indGroups[[g]][["indLabels"]][5]]],
                                                            fig_list[[study_indGroups[[g]][["indLabels"]][6]]],
                                                            fig_list[[study_indGroups[[g]][["indLabels"]][7]]],
                                                            fig_list[[study_indGroups[[g]][["indLabels"]][8]]],
                                                            nrow=4, ncol=2, common.legend = TRUE, legend="bottom"),
             width = 5800, height = 5100, units = myUnits, dpi=ppi)
    } #case2
    
    
  } #end of list of the groups of indices evaluated
  

} #end of function


