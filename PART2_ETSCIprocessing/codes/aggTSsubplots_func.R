#this function aggregates subplots from a list of figures,
#then save it 
    #replace "scale = fscale" by "width = ..., height = ..., units = ".." for a better control on plot size

#CAUTION: HARD CODING HERE BASED ON THE NUMBER OF DATA SOURCES INVOLVED

aggTSsubplots_func <- function(fig_list) {
  
  ggsave(graphFile, arrangeGrob(fig_list[[1]], fig_list[[2]], fig_list[[3]], fig_list[[4]], 
                                fig_list[[5]], fig_list[[6]], fig_list[[7]], fig_list[[8]], 
                                fig_list[[9]], fig_list[[10]], fig_list[[11]], fig_list[[12]], 
                                nrow = 3, ncol = 4, 
                                bottom = text_grob("Years", size = txtSize1, face = txt_style),  
                                top = text_grob(graphTitle, size = txtSize0, face = txt_style, hjust = figT_hjust)), 
         width = myWidth, height = myHeight, units = myUnits, dpi = ppi)
  
}