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
library(car) #for stats stuff like crPLots
library(lmtest) #for assumptions of reg
library(MASS) #for stepAIC
library(lmtest)
library(ggiraphExtra)


#####make a theme###########
theme_incubation <- theme(panel.background = element_rect(fill = "white", linetype = "solid", colour = "black"), 
                          legend.key = element_rect(fill = "white"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
                          axis.title = element_text(size = 30),
                          axis.text=element_text(size=40),                                               #axis numbers
                          plot.title=element_text(size=30, hjust=-0.35),                                 #Main title      --> to center title, hjust=0.5, to left align -0.2
                          legend.title = element_text(size = 16, colour = "black", angle = 0),           #Legend title
                          legend.text = element_text(size = 14, colour = "black", angle = 0),            #Legend text
                          strip.text.x = element_text(size = 50, colour = "black", angle = 0, vjust = 1),        #Facet x text size
                          strip.text.y = element_text(size = 50, colour = "black", angle = 270),
                          strip.background = element_blank(),
                          panel.border = element_blank(),
                          plot.margin = unit(c(1,4,1,1), "lines")
                          )

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
  filter(., !lake_name %in% c("VEE"))

# Make a figure with percent Hg and sulphate ------------------------------

plot.sul.percent <- ggplot(data = mastersheet, aes(x=Sulphate, y=Percent_Hg)) +
  geom_smooth(data=mastersheet, colour="#C6CDF7",method="glm", formula = y~log(x+1),
              method.args = list(family = gaussian(link = 'log')),
              fill = "#C6CDF7") +
  geom_point(size=8, colour = "#7294D4") +
  labs(x = "Sulfate in Water (mg/L)", y = 'Ratio of Methylmercury (%)') +
  expand_limits(x = 0, y = 0) +
  annotate('text', x = 60, y = 1,  label = "p-value<=0.05~~R^{2}==0.39" ,parse = TRUE,size=8) +
  theme_incubation 


plot.sul.percent


# Subset only the rates from the matersheet -------------------------------

subset_1 <- subset(mastersheet, select = c("lake_name", "Percent_Hg", "Nitrogen_Total", "DOC", "pH", "Sulphate", 
                                           "Iron_Total", "Total_Phophorus_UO", "Total_Arsenic_UO",
                                           "km_24", "kd_24", "SUVA_percent", "proton_M")) %>%
  na.omit(.) %>% column_to_rownames(., var="lake_name")



# Make a figure with Km and DOC -------------------------------------------


doc <- ggplot(data = mastersheet, aes(x=(DOC), y=(km_24))) + 
  labs(x = "DOC (mg/L)", y = 'Methylation Rate Constant (h^-1)') +
  geom_smooth(formula = 'y ~ log(x)', method = 'glm', colour="#C6CDF7",
              method.args = list(family = gaussian(link = 'log')), 
              fill = "#C6CDF7") +
  geom_point(size=8, color = "#7294D4", show.legend = FALSE) +
  xlim (0, 70) +
  annotate('text', x = 55, y = 1, label = "p-value<0.05~~R^{2}==0.76" ,parse = TRUE,size=8) +
  theme_incubation 

doc


# Seperate the sampels to the two different groups ------------------------

low.km <- mastersheet %>% filter(lake_name %in% c("RAT", "FRAME", "BCR07",
                                                  "DAVID", "NIVEN", "YKE1",
                                                  "VEE", "BC20", "YK42", "BC21",
                                                  "PROSPEROUS", "YK60",
                                                  "ICING", "BC18", "YK11")) 

low.km["Group"]="Group_1"


high.km <- mastersheet %>% filter(!lake_name %in% c("RAT", "FRAME", "BCR07",
                                                    "DAVID", "NIVEN", "YKE1",
                                                    "VEE", "BC20", "YK42", "BC21",
                                                    "PROSPEROUS", "YK60",
                                                    "ICING", "BC18", "YK11"))

high.km["Group"]="Group_2"


joined.mastersheet <- rbind(low.km, high.km) %>% drop_na(km_24)


# make the figure ---------------------------------------------------------


sul.group <- ggplot(data = joined.mastersheet, aes(x=Sulphate, y=km_24, 
                                                   label = lake_name, color = Group)) + 
  geom_point(size=8, show.legend = T) +
  #geom_text_repel() +
  labs(x = "Sulfate (mg/L)", y = 'Methylation Rate Constant (h^-1)') +
  theme_incubation +
  labs(color = "") +
  theme(legend.position = c(0.8,0.8),
        legend.text = element_text(size = 30, colour = "black", angle = 0),
        legend.title = element_text(size= 30, colour = "black", angle = 0)) +
  scale_color_manual(values = c("#7294D4", "#515c69"), labels = c("Group 1", "Group 2"))


sul.group


# make a figure that has kd and sulfate plus DOC as secondary axis --------

#add legend in post

doc.sul.kd <-  ggplot(mastersheet, show.legend = T)+
  geom_point(aes(x = DOC, y = kd_24, colour ="DOC" ), size = 8, show.legend = T)+                                   
  geom_point(aes(Sulphate, kd_24, colour = "Sulfate"), size = 8, show.legend = T)+                        #Adding points for second X-axis
  scale_color_manual(values = c("#7294D4", "#515c69")) + 
  scale_x_continuous("Sulfate (mg/L)", 
                     sec.axis = sec_axis(~ . *1, name = "DOC (mg/L)")) +
  scale_y_continuous("Demethylation Rate Constant (h^-1)") +
  theme_incubation +
  labs(color = "") +
  theme(legend.position = c(0.8,0.8),
        legend.text = element_text(size = 30, colour = "black", angle = 0),
        legend.title = element_text(size= 30, colour = "black", angle = 0))

doc.sul.kd




# Add them together -------------------------------------------------------

figure <- ggarrange(plot.sul.percent, sul.group, doc, doc.sul.kd, ncol = 2, nrow = 2,
                    labels = c("A", "B", "C", "D"),
                    font.label = list(size = 50, color = "#42434a"),
                    hjust = -0.5, vjust = 2,
                    label.x = 0.11, label.y = 0.90,
                    align="hv")


figure

ggsave("Figure_3(composite).png", plot=figure, width = 58.5, height = 45, 
       units = c("cm"), path = "Figures/")


########================================================== DONE ================================#################
