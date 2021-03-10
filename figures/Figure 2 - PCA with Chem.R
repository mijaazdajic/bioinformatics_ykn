###############################################################################################
#####THIS R SCRIPT IS FOR MAKING A BIG COMPOSITE FIGURE OF RESULTS#########################
#####We will use data made in the last three scripts and a csv file with chem data#############
###############################################################################################


# Load Libraries ----------------------------------------------------------
library(tidyverse) # for basic "tidy" workflows and data wrangling
library(janitor) # lots of help when dealing with non-R-friendly datasets
library(lubridate) # for working with dates
library(ggplot2) #plotting
library(ggpubr) #probably don't need it
library(lmtest)
library(ggiraphExtra)
library(vegan)
library(ggfortify)
library(ggrepel)
library(ggpubr)


#####make a theme###########
theme_incubation <- theme(panel.background = element_rect(fill = "white", linetype = "solid", colour = "black"), 
                          legend.key = element_rect(fill = "white"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
                          axis.text=element_text(size=17),                                               #axis numbers
                          plot.title=element_text(size=30, hjust=-0.35),                                 #Main title      --> to center title, hjust=0.5, to left align -0.2
                          legend.title = element_text(size = 16, colour = "black", angle = 0),           #Legend title
                          legend.text = element_text(size = 14, colour = "black", angle = 0),            #Legend text
                          strip.text.x = element_text(size = 20, colour = "black", angle = 0, vjust = 1),        #Facet x text size
                          strip.text.y = element_text(size = 20, colour = "black", angle = 270),
                          strip.background = element_blank(),
                          panel.border = element_blank())   

# Load Water Chemistry Data -----------------------------------------------
# Make a mastersheet with chem and methylation ----------------------------


water_chem <- read_csv("CHEM_DATA/COMPLETE_MASTER_SHEET_ALL_WATER_CHEM.csv",local = locale(encoding = "latin1"))

distance <- read_csv("INFO_SHEETS/Distance_Lakes.csv", local = locale(encoding = "latin1"))

suva <- read_csv("CHEM_DATA/SUVA_Analysis.csv", local = locale(encoding = "latin1"))

years <- read_csv("INFO_SHEETS/Sampling_Year_Per_Lake.csv", locale = locale(encoding = "latin1"))

methylation_table <- read_csv("methylation_table.csv", local = locale(encoding = "latin1"))

methylation_table <- methylation_table %>% merge(.,years, by = "lake_name", all=T)

mastersheet <- left_join(water_chem, methylation_table, by = c("lake_name","Year")) %>% 
  left_join(., suva, by = "lake_name", all=T) %>% mutate(., proton_M = 10^(-pH)) %>% 
  filter(., !lake_name %in% c("POCKET")) %>% 
  transform(., Year = as.character(Year))


# Make a PCA --------------------------------------------------------------


df <- mastersheet %>% select(Sulphate, pH, DOC, Iron_Total, Total_Phophorus_UO, SUVA_percent)
pca_yk <- prcomp(df, scale. = TRUE)

pca_df1 <- data.frame(pca_yk$rotation[,1:2])



plot_2 <- autoplot(pca_yk, data = mastersheet, colour = 'Year', size = 8) +
  #geom_text_repel(aes(label = lake_name), colour = "black",size=3, point.padding = 0.25) +
  geom_segment(data=pca_df1, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="grey50", arrow=arrow(length=unit(0.01,"npc"))) +
#  geom_text(data=pca_df1, aes(x=PC1,y=PC2,label=rownames(pca_df1),
#                          hjust=0.5*(1-sign(PC1)),vjust=0.5*(1-sign(PC2))), 
#            color="grey50", size=6
#            ) +
  geom_text_repel(data = pca_df1, aes(label = rownames(pca_df1)), colour = "grey50",size=6, point.padding = 0.25) +
  coord_fixed() +
  theme_incubation +
  theme(legend.position = "right",
        legend.title = element_blank())


plot_2

ggsave(plot_2, file="Figure_2.svg", width = 10, height=10, units = "in", path = "Figures/")




