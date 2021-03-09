###############################################################################################
#####THIS R SCRIPT IS FOR PLOTTING INCUBATION RESULTS FROM ICP-MS##############################
#####Each incubation is saved as a CSV file, with a number, the number key is in the folder####
###############################################################################################

#####make a theme###########
theme_incubation <- theme(panel.background = element_rect(fill = "white", linetype = "solid", 
                                                          colour = "black"), 
                   legend.key = element_rect(fill = "white"), panel.grid.minor = element_blank(), 
                          panel.grid.major = element_blank(), 
                   axis.text=element_text(size=17),                 #axis numbers
                   plot.title=element_text(size=30, hjust=-0.35),   #Main title      --> to center title, hjust=0.5, to left align -0.2
                   legend.title = element_text(size = 16, colour = "black", angle = 0), #Legend title
                   legend.text = element_text(size = 14, colour = "black", angle = 0),  #Legend text
                   strip.text.x = element_text(size = 20, colour = "black", 
                                               angle = 0, vjust = 1),        #Facet x text size
                   strip.text.y = element_text(size = 20, colour = "black", angle = 270),
                   strip.background = element_blank(),
                   panel.border = element_blank())   

###=====Packages needed=====####
library(tidyverse) # for basic "tidy" workflows and data wrangling
library(lubridate) # for working with dates
library(cowplot) # for arranging multiple plots (and many more)
library(janitor) # lots of help when dealing with non-R-friendly datasets
library(ggplot2)
library(plyr)
library(dplyr)
library(grid)
library(ggforce)
library(ggpubr)
library(glue)

# Multiple plot function####

#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

####Choosing only SP points and no USP####

incubation.sp = incubation.join %>% dplyr::filter(!str_detect(treatment, "USP"))

#write.csv(incubation.sp, file = "SHEETS_FROM_EACH_SCRIPT/02_incubation_sp.csv")

spiked_conc = incubation.sp %>% aggregate(spiked_con ~ lake_name,. , mean )

#write.csv(spiked_conc, file = "SHEETS_FROM_EACH_SCRIPT/02_spiked_conc.csv")

# Make a figure for each lake in loops --------------------------------------------

plots = list()

var_list <- unique(incubation.sp$lake_name) #make a list of lake names

for (i in 1:length(var_list)) {
  incubation.temp <- incubation.sp[incubation.sp$lake_name==var_list[i],] #make a temporary file, with info about that lake
  spiked.temp <- spiked_conc[spiked_conc$lake_name==var_list[i],]
  
  methylation <- ggplot(data = incubation.temp, aes(x=as.numeric(time), 
                                                    y=X199_corrected_fixed_Conc, shape = Rep)) +
    geom_vline(xintercept = 24, colour='#a8dadc', size=2) +
    geom_vline(xintercept = 48, colour='#457b9d', size=2) +
    geom_line() +
    geom_point() +
    labs(x = "Time (hour)", y = 'Me198Hg Conc (ug/kg)') +
    theme_incubation
  
  demethylation = ggplot(data = incubation.temp, aes(x=as.numeric(time), 
                                                     y=X198_corrected_fixed_Conc, shape = Rep)) +
    geom_vline(xintercept = 24, colour='#a8dadc', size=2) +
    geom_vline(xintercept = 48, colour='#457b9d', size=2) +
    geom_hline(aes(yintercept = spiked.temp$spiked_con, colour="#E63946"), 
               data=spiked_conc, show.legend = F) + #
    expand_limits(x = 0, y = 0) +
    geom_line() +
    geom_point() +
    labs(x = "Time (hour)", y = 'Me198Hg Conc (ug/kg)') +
    theme_incubation
  
  
  
  figure <- ggarrange(methylation,demethylation, nrow = 2,
                      font.label = list(size = 25),
                      hjust = 0, vjust = 2,
                      label.x = 0.13) 
  figure <- annotate_figure(figure, top = text_grob(var_list[i], size = 14))
  
  ggsave(figure, file=glue(var_list[i]," Incubation",".svg"), dpi = 600, path = "Figures/")
  
  plots[[i]] = figure
  
}



# Average the reps and then make a new figure -----------------------------

incubation_mean = incubation.sp %>% dplyr::select(lake_name,time,treatment, X199_corrected_fixed_Conc, 
                                                  X198_corrected_fixed_Conc, ln.198.) %>%
  dplyr::group_by(lake_name, add=T) %>% dplyr::group_by(time, add=T) %>% 
  aggregate(., by=list(.$lake_name, .$time), FUN=mean) %>% 
  dplyr::select(.,-c(lake_name,time, treatment)) %>%
  rename(.,c("Group.2" = "time")) %>% 
  rename(.,c("Group.1" = "lake_name"))
  
#write.csv(incubation_mean, file = "SHEETS_FROM_EACH_SCRIPT/02_incubation_mean.csv")
  
spiked_MA300 = spike_all %>% dplyr::select(lake_name,corrected_MA300) %>% #take only the data from MA3000
  mutate(corrected_MA300 = if_else(is.na(corrected_MA300), 0, corrected_MA300)) #make sure that there in no NA, if there is, replace by 0, otherwise the loop will stop

spiked_ICP = spike_all %>% dplyr::select(lake_name,corrected_ICP) %>% #take only the data from MA3000
  mutate(corrected_ICP = if_else(is.na(corrected_ICP), 0, corrected_ICP)) #make sure that there in no NA, if there is, replace by 0, otherwise the loop will stop

# Make methylation and demethylation figures with means and MA3000 --------------------

plots = list() #make an empty list of plots

var_list <- unique(incubation_mean$lake_name) #make a list of lake names

spike_recovery <- matrix(ncol = 2)
colnames(spike_recovery) <-  c("lake_name","recovery") #make an empty matrix to add all the recovery results

for (i in 1:length(var_list)) {
  incubation.temp <- incubation_mean[incubation_mean$lake_name==var_list[i],] #make a temporary file, with info about that lake
  spiked.temp <- spiked_MA300[spiked_MA300$lake_name==var_list[i],] #make a temporary file with spiked data for specific lake
  spike <- pull(spiked.temp,corrected_MA300) #pull up the concentration
  initial.spike <- incubation.temp %>% subset(.,time==0) %>% pull(.,X198_corrected_fixed_Conc) #pull up the analysed cocentration of 198 at time 0
  recovery <- as.numeric(initial.spike/spike*100) %>% format(round(.,6)) #caculate the % recovery
  
  methylation <- ggplot(data = incubation.temp, aes(x=as.numeric(time), y=X199_corrected_fixed_Conc)) + 
    geom_vline(xintercept = 24, colour='#a8dadc', size=2) +
    geom_vline(xintercept = 48, colour='#457b9d', size=2) +
    geom_line() +
    geom_point() +
    labs(x = "Time (hour)", y = 'Me198Hg Conc (ug/kg)') +
    theme_incubation
  
  demethylation = ggplot(data = incubation.temp, aes(x=as.numeric(time), y=X198_corrected_fixed_Conc)) +
    geom_vline(xintercept = 24, colour='#a8dadc', size=2) +
    geom_vline(xintercept = 48, colour='#457b9d', size=2) +
    geom_hline(aes(yintercept = spiked.temp$corrected_MA300, colour="#E63946"), data=spiked_MA300, show.legend = F) + #
    expand_limits(x = 0, y = 0) +
    geom_line() +
    geom_point() +
    geom_label(aes(label = paste("recovery = ", recovery,"%", sep="")), x=40, y=c(spike), label.size = NA) +
    labs(x = "Time (hour)", y = 'Me198Hg Conc (ug/kg)') +
    theme_incubation
  
  spike_recovery <- rbind(spike_recovery, c(var_list[i], recovery))
  
  figure <- ggarrange(methylation,demethylation, nrow = 2,
                      font.label = list(size = 25),
                      hjust = 0, vjust = 2,
                      label.x = 0.13) 
  figure <- annotate_figure(figure, top = text_grob(var_list[i], size = 14))
  
  ggsave(figure, file=glue(var_list[i],"_Incubation_mean_MA300",".svg"), dpi = 600, path = "Figures/MA3000/")
  
  plots[[i]] = figure
  
}

write.csv(spike_recovery, file = "spike_recovery_MA3000.csv")

# Make methylation and demethylation figures with means and ICP--------------------

plots = list() #make an empty list of plots

var_list <- unique(incubation_mean$lake_name) #make a list of lake names

spike_recovery <- matrix(ncol = 2)
colnames(spike_recovery) <-  c("lake_name","recovery")

for (i in 1:length(var_list)) {
  incubation.temp <- incubation_mean[incubation_mean$lake_name==var_list[i],] #make a temporary file, with info about that lake
  spiked.temp <- spiked_ICP[spiked_ICP$lake_name==var_list[i],] #make a temporary file with spiked data for specific lake
  spike <- pull(spiked.temp,corrected_ICP) #pull up the concentration
  initial.spike <- incubation.temp %>% subset(.,time==0) %>% 
  pull(.,X198_corrected_fixed_Conc) #pull up the analysed cocentration of 198 at time 0
  recovery <- as.numeric(initial.spike/spike*100) %>% 
  format(round(.,6)) #caculate the % recovery
  
  methylation <- ggplot(data = incubation.temp, aes(x=as.numeric(time), 
                                                    y=X199_corrected_fixed_Conc)) + 
    geom_vline(xintercept = 24, colour='#a8dadc', size=2) +
    geom_vline(xintercept = 48, colour='#457b9d', size=2) +
    geom_line() +
    geom_point() +
    labs(x = "Time (hour)", y = 'Me198Hg Conc (ug/kg)') +
    theme_incubation
  
  demethylation = ggplot(data = incubation.temp, aes(x=as.numeric(time), 
                                                     y=X198_corrected_fixed_Conc)) +
    geom_vline(xintercept = 24, colour='#a8dadc', size=2) +
    geom_vline(xintercept = 48, colour='#457b9d', size=2) +
    geom_hline(aes(yintercept = spiked.temp$corrected_ICP, colour="#E63946"), 
               data=spiked_ICP, show.legend = F) + #
    expand_limits(x = 0, y = 0) +
    geom_line() +
    geom_point() +
    geom_label(aes(label = paste("recovery = ", recovery,"%", sep="")), 
               x=40, y=c(spike), label.size = NA) +
    labs(x = "Time (hour)", y = 'Me198Hg Conc (ug/kg)') +
    theme_incubation
  
  spike_recovery <- rbind(spike_recovery, c(var_list[i], recovery))
  
  figure <- ggarrange(methylation,demethylation, nrow = 2,
                      font.label = list(size = 25),
                      hjust = 0, vjust = 2,
                      label.x = 0.13) 
 figure <- annotate_figure(figure, top = text_grob(var_list[i], size = 14))
  
 ggsave(figure, file=glue(var_list[i],"_Incubation_mean_ICP",".pdf"), 
        dpi = 600, path = "Figures/ICP/")
  
  plots[[i]] = figure
  
}

write.csv(spike_recovery, file = "spike_recovery_ICP.csv")


####-------------------------------------------------- END ------------------------------------------------####
