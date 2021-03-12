###############################################################################################
#####THIS R SCRIPT IS FOR PLOTTING INCUBATION RESULTS FROM ICP-MS##############################
###################THIS SCRIPT IS TO MAKE FIGURES WITH REGRESSION LINES AND RATES##############
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

# Load Libraries ----------------------------------------------------------
library(tidyverse) # for basic "tidy" workflows and data wrangling
library(janitor) # lots of help when dealing with non-R-friendly datasets
library(lubridate) # for working with dates
library(ggplot2) #plotting
library(ggpubr) #probably don't need it



# ICP seems to be better for looking at the data recovery -----------------
# Will calculate Km and Kd for 24 hours -----------------------------------
# Make methylation and demethylation figures with means and ICP------------

incubation_mean_24 = incubation.sp %>% 
  dplyr::select(lake_name,time,treatment, X199_corrected_fixed_Conc, X198_corrected_fixed_Conc, ln.198.) %>%
  dplyr::group_by(lake_name, add=T) %>% dplyr::group_by(time, add=T) %>% 
  stats::aggregate.data.frame(., by = list(.$lake_name,.$time), FUN=mean) %>% 
  dplyr::select(.,-c(lake_name,treatment, time)) %>% 
  rename(.,c("Group.2" = "time")) %>% 
  rename(.,c("Group.1" = "lake_name")) %>% 
  dplyr::filter(time != 48) %>% #take out the 48 hour data 
  dplyr::mutate(ln.198. = replace(ln.198., ln.198. == "-Inf", log(0.0000001)))

#write.csv(incubation_mean_24, file = "SHEETS_FROM_EACH_SCRIPT/03_incubation_mean_24.csv")

spiked_ICP = spike_all %>% dplyr::select(lake_name,corrected_ICP) %>% #take only the data from MA3000
  mutate(corrected_ICP = if_else(is.na(corrected_ICP), 0, corrected_ICP)) #make sure that there in no NA, if there is, replace by 0, otherwise the loop will stop

plots = list() #make an empty list of plots

var_list <- unique(incubation_mean_24$lake_name) #make a list of lake names

rates <- matrix(ncol = 3)
colnames(rates) <-  c("lake_name","methylation","demethylation")

spike_recovery <- matrix(ncol = 2)
colnames(spike_recovery) <-  c("lake_name","recovery") #make an empty matrix to add all the recovery results

for (i in 1:length(var_list)) {
  incubation.temp <- incubation_mean_24 %>% filter(lake_name == var_list[i]) #make a temporary file, with info about that lake
  spiked.temp <- spiked_ICP[spiked_ICP$lake_name==var_list[i],] #make a temporary file with spiked data for specific lake
  spike <- pull(spiked.temp,corrected_ICP) #pull up the concentration
  initial.spike <- incubation.temp %>% subset(.,time==0) %>% pull(.,X198_corrected_fixed_Conc) #pull up the analysed cocentration of 198 at time 0
  recovery <- as.numeric(initial.spike/spike*100) %>% format(round(.,6)) #caculate the % recovery
  
  methylation <- ggplot(data = incubation.temp, aes(x=as.numeric(time), y=X199_corrected_fixed_Conc)) + 
    geom_vline(xintercept = 24, colour='#a8dadc', size=2) +
    geom_smooth(method = lm, formula = y~x, color = "#62BB35", fill = "lightgrey") +
    geom_point(aes(size=4), show.legend = FALSE) +
    labs(x = "Time (hour)", y = 'Me198Hg Conc (ug/kg)') +
    theme_incubation
  
  demethylation = ggplot(data = incubation.temp, aes(x=as.numeric(time), y=X198_corrected_fixed_Conc)) +
    geom_vline(xintercept = 24, colour='#a8dadc', size=2) +
    geom_hline(aes(yintercept = spiked.temp$corrected_ICP, colour="#E63946"), data=spiked_ICP, show.legend = F) + #
    expand_limits(x = 0, y = 0) +
    geom_smooth(method = lm, formula = y~x, fill = "lightgrey", color = "#4178BC") +
    geom_point(aes(size=4), show.legend = FALSE) +
    geom_label(aes(label = paste("recovery = ", recovery,"%", sep="")), x=40, y=c(spike), label.size = NA) +
    labs(x = "Time (hour)", y = 'Me198Hg Conc (ug/kg)') +
    theme_incubation
  

  methylation.reg <- lm(X199_corrected_fixed_Conc~as.numeric(time), data=incubation.temp)
                      
  methylation.reg <- summary(methylation.reg)$coefficients[2,1]

  demethylation.reg <- lm(ln.198.~as.numeric(time), data=incubation.temp)
  
  demethylation.reg <- summary(demethylation.reg)$coefficients[2,1]
  
  rates <- rbind(rates, c(var_list[i], methylation.reg, demethylation.reg))
  
  
  spike_recovery <- rbind(spike_recovery, c(var_list[i], recovery))
  
  figure <- ggarrange(methylation,demethylation, nrow = 2,
                      font.label = list(size = 25),
                      hjust = 0, vjust = 2,
                      label.x = 0.13) 
  figure <- annotate_figure(figure, top = text_grob(var_list[i], size = 14))
  
  ggsave(figure, file=glue(var_list[i],"_Incubation_mean_ICP",".svg"), dpi = 600, path = "Figures/ICP/")
  
  plots[[i]] = figure
  
}

write.csv(spike_recovery, file = "spike_recovery_ICP.csv")

write.csv(rates, file = "incubation_rates_24.csv")


# ICP seems to be better for looking at the data recovery -----------------
# Will calculate Km and Kd for 48 hours -----------------------------------
# Make methylation and demethylation figures with means and ICP------------

incubation_mean_48 = incubation.sp %>% 
  dplyr::select(lake_name,time,treatment, X199_corrected_fixed_Conc, X198_corrected_fixed_Conc, ln.198.) %>%
  dplyr::group_by(lake_name, add=T) %>%  dplyr::group_by(time, add=T) %>% 
  stats::aggregate(., by = list(.$lake_name,.$time), FUN=mean) %>% 
  dplyr::select(.,-c(lake_name,treatment, time)) %>% 
  rename(.,c("Group.2" = "time")) %>% 
  rename(.,c("Group.1" = "lake_name")) %>% 
  dplyr::mutate(ln.198. = replace(ln.198., ln.198. == "-Inf", log(0.0000001)))

#write.csv(incubation_mean_48, file = "SHEETS_FROM_EACH_SCRIPT/02_incubation_mean_48.csv")

spiked_ICP = spike_all %>% dplyr::select(lake_name,corrected_ICP) %>% #take only the data from MA3000
  dplyr::mutate(corrected_ICP = if_else(is.na(corrected_ICP), 0, corrected_ICP)) #make sure that there in no NA, if there is, replace by 0, otherwise the loop will stop

plots = list() #make an empty list of plots

var_list <- unique(incubation_mean_48$lake_name) #make a list of lake names

rates <- matrix(ncol = 3)
colnames(rates) <-  c("lake_name","methylation","demethylation")

spike_recovery <- matrix(ncol = 2)
colnames(spike_recovery) <-  c("lake_name","recovery") #make an empty matrix to add all the recovery results

for (i in 1:length(var_list)) {
  incubation.temp <- incubation_mean_48 %>% filter(lake_name == var_list[i]) #make a temporary file, with info about that lake
  spiked.temp <- spiked_ICP[spiked_ICP$lake_name==var_list[i],] #make a temporary file with spiked data for specific lake
  spike <- pull(spiked.temp,corrected_ICP) #pull up the concentration
  initial.spike <- incubation.temp %>% subset(.,time==0) %>% pull(.,X198_corrected_fixed_Conc) #pull up the analysed cocentration of 198 at time 0
  recovery <- as.numeric(initial.spike/spike*100) %>% format(round(.,6)) #caculate the % recovery
  
  methylation <- ggplot(data = incubation.temp, aes(x=as.numeric(time), y=X199_corrected_fixed_Conc)) + 
    geom_vline(xintercept = 24, colour='#a8dadc', size=2) +
    geom_smooth(method = lm, formula = y~x, color = "#62BB35", fill = "lightgrey") +
    geom_point(aes(size=4), show.legend = FALSE) +
    labs(x = "Time (hour)", y = 'Me198Hg Conc (ug/kg)') +
    theme_incubation
  
  demethylation = ggplot(data = incubation.temp, aes(x=as.numeric(time), y=X198_corrected_fixed_Conc)) +
    geom_vline(xintercept = 24, colour='#a8dadc', size=2) +
    geom_hline(aes(yintercept = spiked.temp$corrected_ICP, colour="#E63946"), data=spiked_ICP, show.legend = F) + #
    expand_limits(x = 0, y = 0) +
    geom_smooth(method = lm, formula = y~x, fill = "lightgrey", color = "#4178BC") +
    geom_point(aes(size=4), show.legend = FALSE) +
    geom_label(aes(label = paste("recovery = ", recovery,"%", sep="")), x=40, y=c(spike), label.size = NA) +
    labs(x = "Time (hour)", y = 'Me198Hg Conc (ug/kg)') +
    theme_incubation
  
  
  methylation.reg <- lm(X199_corrected_fixed_Conc~as.numeric(time), data=incubation.temp)
  
  methylation.reg <- summary(methylation.reg)$coefficients[2,1]
  
  demethylation.reg <- lm(ln.198.~as.numeric(time), data=incubation.temp)
  
  demethylation.reg <- summary(demethylation.reg)$coefficients[2,1]
  
  rates <- rbind(rates, c(var_list[i], methylation.reg, demethylation.reg))
  
  
  spike_recovery <- rbind(spike_recovery, c(var_list[i], recovery))
  
  figure <- ggarrange(methylation,demethylation, nrow = 2,
                      font.label = list(size = 25),
                      hjust = 0, vjust = 2,
                      label.x = 0.13) 
  figure <- annotate_figure(figure, top = text_grob(var_list[i], size = 14))
  
  ggsave(figure, file=glue(var_list[i],"_Incubation_mean_ICP",".svg"), dpi = 600, path = "Figures/ICP/")
  
  plots[[i]] = figure
  
}

write.csv(spike_recovery, file = "spike_recovery_ICP.csv")

write.csv(rates, file = "incubation_rates_48.csv")



# Calculating Rate Constants ----------------------------------------------
# The caclulations for this section of Km/Kd is found at Zhu et al (2018) -



incubation_rates_24 <- read_csv("incubation_rates_24.csv",local = locale(encoding = "latin1")) %>%
                        dplyr::select(-"X1") %>% 
                        na.omit(.) %>% merge(.,spiked_ICP, by = "lake_name", all=T) %>%
                        mutate(km_24 = ((methylation/corrected_ICP)*24)) %>% mutate(kd_24 = demethylation) %>% #this calculation comes from Zhu et al. 2018
                        mutate(methylation_24 = methylation) %>% mutate(demethylation_24 = demethylation) %>%
                        dplyr::select("lake_name","methylation_24","demethylation_24","km_24","kd_24", -c(methylation, demethylation)) %>% 
                        mutate(., kd_24 = ifelse(kd_24 >= 0, 0, kd_24)) %>%
                        mutate(kd_24 = (kd_24*-24)) %>%
                        mutate(., demethylation_24 = ifelse(demethylation_24 >= 0 , 0, demethylation_24))



incubation_rates_48 <- read_csv("incubation_rates_48.csv",local = locale(encoding = "latin1")) %>%
                        dplyr::select(-"X1") %>% na.omit(.)  %>% merge(.,spiked_ICP, by = "lake_name", all=T) %>%
                        mutate(km_48 = (methylation/corrected_ICP)*24) %>% mutate(kd_48 = demethylation) %>% #this calculation comes from Zhu et al. 2018
                        mutate(methylation_48 = methylation) %>% mutate(demethylation_48 = demethylation) %>%
                        dplyr::select("lake_name","methylation_48","demethylation_48","km_48","kd_48", -c(methylation, demethylation)) %>% 
                        mutate(., kd_48 = ifelse(kd_48 >= 0, 0, kd_48)) %>%
                        mutate(kd_48 = (kd_48*-24)) %>%
                        mutate(., demethylation_48 = ifelse(demethylation_48 >= 0 , 0, demethylation_48))

#methylation_table <- merge(incubation_rates_24, incubation_rates_48, by="lake_name", all=T)

#write.csv(methylation_table, file = "methylation_table.csv")

                        
########================================================== DONE ================================#################
                        
