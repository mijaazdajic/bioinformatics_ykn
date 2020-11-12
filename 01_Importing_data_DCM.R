
# Importing data from DCM extraction --------------------------------------

###############################################################################################
#####THIS R SCRIPT IS FOR IMPORTING THE MOISTURE CONTENT AND ISOTOPE RESULTS FROM ICP-MS#######
#####Each incubation is saved as a CSV file, with a number, the number key is in the folder####
###############################################################################################

###=====Packages needed=====####
library(tidyverse) # for basic "tidy" workflows and data wrangling
library(lubridate) # for working with dates
library(cowplot) # for arranging multiple plots (and many more)
library(janitor) # lots of help when dealing with non-R-friendly datasets


###=====Import Data=====####
###before importing, clean up the excel sheet, only one row of headers###
##isotope data###

temp = list.files(path="RAW_DCM/.", pattern="*.csv")
All <- lapply(temp,function(i){
  read_csv(paste("RAW_DCM/",i,sep=""),local = locale(encoding = "latin1"))
})


dat = c()

for (i in 1:length(All)){
  #i=1
  dat <- bind_rows(dat, All[[i]])
}

incubation <- dat %>%  
  mutate(date = dmy(paste(Day, Month, Year, sep = "-"))) %>%
  dplyr::select(date, everything(), -c(Day, Month, Year)) %>%
  dplyr::filter(str_detect(Sample_Name, "SP")) %>%
  dplyr::filter(!str_detect(Sample_Name, "198")) %>%
  mutate(lake_name = str_extract(Sample_Name, pattern="[[:alnum:]]*(?=[:blank:])"),
         condition = str_extract(Sample_Name, pattern="(?<=[:blank:])[[:alnum:]]*")) %>%
  mutate(time = str_extract(condition, pattern="[[:digit:]]*$"),
         treatment = str_extract(condition, pattern="[[:alpha:]]*")) %>%
  dplyr::select("lake_name","treatment","time",everything(), -c(condition,Sample_Name)) %>%
  mutate_all(., toupper) %>% 
  mutate(spiked_con = (as.numeric(Volume)*as.numeric(Spike_Conc))/as.numeric(Weight_sample))

##moisture data###
temp = list.files(path="Moisture/.", pattern="*.csv")
Moisture <- lapply(temp,function(i){
  read_csv(paste("Moisture/",i,sep=""),local = locale(encoding = "latin1"))
})


dat = c()

for (i in 1:length(Moisture)){
  #i=1
  dat <- bind_rows(dat, Moisture[[i]])
} #when there is a error about characters then it is probably becuase you have one of the years writte as such: "20. 19"

Moisture = dat %>% 
  mutate(date = dmy(paste(Day, Month, Year, sep = "-"))) %>%
  dplyr::select(date, everything(), -c(Day, Month, Year)) %>%
  mutate(lake_name = str_extract(Sample_Name, pattern="[[:alnum:]]*(?=[:blank:])"),
         condition = str_extract(Sample_Name, pattern="(?<=[:blank:])[[:alnum:]]*")) %>%
  mutate(treatment = str_extract(condition, pattern="[[:alpha:]]*")) %>%
  dplyr::select("lake_name","treatment",everything(), -c(condition,Sample_Name)) %>%
  mutate_all(., toupper)

###=====Manipulating Data=====####

incubation.join = left_join(x=incubation, y=Moisture, by = c("lake_name", "treatment", "date", "Rep")) #join the moisture data to the incbation data

incubation.join <- incubation.join %>%                                      #correct the concentrations to the moistuer content
  mutate(moisture_factor = (100-as.numeric(Moisture_Content)) / 100) %>% 
  mutate(corrected_198 = (moisture_factor*as.numeric(`198_Hg_ug_kg`))) %>% 
  mutate(corrected_199 = (moisture_factor*as.numeric(`199_Hg_ug_Kg`))) %>% 
  mutate(corrected_202 = (moisture_factor*as.numeric(`202_Hg_ug_kg`)))

incubation.join <- incubation.join %>% 
  mutate(X199_corrected_fixed_Conc=ifelse(corrected_199<0, 0, corrected_199)) %>% 
  mutate(X198_corrected_fixed_Conc=ifelse(corrected_198<0, 0, corrected_198)) %>% 
  mutate(ln.198.=log(X198_corrected_fixed_Conc))

#write.csv(incubation.join, file = "SHEETS_FROM_EACH_SCRIPT/01_incubation_join.csv")


###=====Import Spike Data=====####
Moisture$Moisture_Content = as.numeric(Moisture$Moisture_Content) #make sure the dataset is numeric

moisture_mean <- Moisture %>% aggregate(Moisture_Content ~ lake_name,. , mean ) #average the treatments into one


spike_all = read_csv("Spike_Concentations.csv", col_names = T, local = locale(encoding = "latin1"), col_types = cols(.default = "d", lake_name = "c")) %>% #correct the concentrations to the moistuer content
  mutate_all(., toupper) %>%      left_join(x=., y=moisture_mean, by = c("lake_name")) %>% #make all row names capitalized, add the moisture data to the frame
  mutate(corrected_MA300 = (as.numeric(Volume)*as.numeric(Spike_MA))/as.numeric(Mass)) %>% #calculate the cocentation of spike using volume and mass 
  mutate(corrected_ICP = (as.numeric(Volume)*as.numeric(Spike_ICP))/as.numeric(Mass))

#write.csv(spike_all, file = "SHEETS_FROM_EACH_SCRIPT/01_spike_all.csv")

########================================================== DONE ================================#################


