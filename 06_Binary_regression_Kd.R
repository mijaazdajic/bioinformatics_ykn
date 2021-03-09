###############################################################################################
#######THIS R SCRIPT IS FOR MERGING METHYLATION AND CHEM DATA TO ANALYSE#######################
########this script is to identify variables sig correlated with Km in#########################
#####We will use data made in the last three scripts and a csv file with chem data#############
###########this script is used to deal with the large amount of zeros #########################


# Load Libraries ----------------------------------------------------------
library(tidyverse) # for basic "tidy" workflows and data wrangling
library(janitor) # lots of help when dealing with non-R-friendly datasets
library(lubridate) # for working with dates
library(ggplot2) #plotting
library(ggpubr) #probably don't need it
library(car) #for stats stuff like crPLots
library(lmtest) #for assumptions of reg
library(MASS) #for stepAIC
library(ggiraphExtra)


# Load Water Chemistry Data -----------------------------------------------
# Make a mastersheet with chem and methylation ----------------------------


water_chem <- read_csv("CHEM_DATA/COMPLETE_MASTER_SHEET_ALL_WATER_CHEM.csv",local = locale(encoding = "latin1"))

suva <- read_csv("CHEM_DATA/SUVA_Analysis.csv", local = locale(encoding = "latin1"))

years <- read_csv("INFO_SHEETS/Sampling_Year_Per_Lake.csv", locale = locale(encoding = "latin1"))

methylation_table <- read_csv("methylation_table.csv", local = locale(encoding = "latin1"))

methylation_table <- methylation_table %>% merge(.,years, by = "lake_name", all=T)

mastersheet <- left_join(water_chem, methylation_table, by = c("lake_name","Year")) %>% 
  left_join(., suva, by = "lake_name", all=T) %>% mutate(., proton_M = 10^(-pH))

#write.csv(mastersheet, file = "SHEETS_FROM_EACH_SCRIPT/04_mastersheet.csv")

# Subset only the rates from the matersheet -------------------------------

subset_1 <- subset(mastersheet, select = c("lake_name", "Percent_Hg", "Nitrogen_Total", "DOC", "pH", "Sulphate", 
                                           "Iron_Total", "Total_Phophorus_UO", "Total_Arsenic_UO",
                                           "km_24", "kd_24", "SUVA_percent", "proton_M")) %>%
  na.omit(.) 



# Seperate the data into two categories -----------------------------------

zero <- subset_1 %>% subset(., kd_24 == "0") %>% mutate(kd_observed=ifelse((kd_24 = 0),"NO","YES"))

non_zero <- subset_1 %>% subset(., kd_24 !=  "0") %>% mutate(kd_observed=ifelse((kd_24 = 0),11,"NO"))

categorical <- rbind(zero, non_zero)

# Do a regression on non zero dataset -------------------------------------

scatterplotMatrix(~logkd + log(Sulphate+1) + SUVA_percent + log(Total_Arsenic_UO) + 
                    log(proton_M) + log(DOC) + log(Iron_Total) + 
                    log(Total_Phophorus_UO+1), reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = non_zero) #decided to make sure they were normally distributed before anlaysis

non_zero <- non_zero %>% mutate(., logkd = sqrt(kd_24))


full.model <- lm(logkd ~ log(Sulphate+1) + SUVA_percent + log(Total_Arsenic_UO) + 
                   pH + log(DOC) + log(Iron_Total) + 
                   log(Total_Phophorus_UO+1), data = non_zero)
step.model <- stepAIC(full.model, direction = "both") #none
summary(step.model)
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=subset_1)
resettest(step.model, power = 2:3, type = "regressor", data = subset_1)
shapiro.test(step.model$residuals) #did not pass = p=0.03, not great




# Make each boxplot sperately ---------------------------------------------



# t-test to see if different ----------------------------------------------
#sulfate
categorical <- categorical %>% mutate(., logsul = log(Sulphate+1))


model.zero <- lm(logsul ~ kd_observed, data = categorical)
summary(model.zero)
leveneTest(logsul ~ kd_observed, data = categorical) #good
shapiro.test(model.zero$residuals) 

#SUVA

model.zero <- lm(SUVA_percent ~ kd_observed, data = categorical)
summary(model.zero)
leveneTest(logsul ~ kd_observed, data = categorical) #good
shapiro.test(model.zero$residuals) 

#Arsenic
model.zero <- lm(log(Total_Arsenic_UO) ~ kd_observed, data = categorical)
summary(model.zero)
leveneTest(logsul ~ kd_observed, data = categorical) #good
shapiro.test(model.zero$residuals) #meh

#ph
model.zero <- lm((l) ~ kd_observed, data = categorical)
summary(model.zero)
leveneTest(logsul ~ kd_observed, data = categorical) #good
shapiro.test(model.zero$residuals) 

#DOC
model.zero <- lm(DOC ~ kd_observed, data = categorical)
summary(model.zero)
leveneTest(logsul ~ kd_observed, data = categorical) #good
shapiro.test(model.zero$residuals) 

#Iron
model.zero <- lm(log(Iron_Total) ~ kd_observed, data = categorical)
summary(model.zero)
leveneTest(logsul ~ kd_observed, data = categorical) #good
shapiro.test(model.zero$residuals) 

#Phos
model.zero <- lm(log(Total_Phophorus_UO+1) ~ kd_observed, data = categorical)
summary(model.zero)
leveneTest(logsul ~ kd_observed, data = categorical) #good
shapiro.test(model.zero$residuals) 
