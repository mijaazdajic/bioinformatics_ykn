###############################################################################################
#####THIS R SCRIPT IS FOR MERGING METHYLATION AND CHEM DATA FOR STATS##########################
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
library(lmtest) #residual testing
library(ggiraphExtra)


#####make a theme###########
theme_incubation <- theme(panel.background = element_rect(fill = "white", 
                                                          linetype = "solid", colour = "black"), 
                          legend.key = element_rect(fill = "white"), 
                          panel.grid.minor = element_blank(), 
                          panel.grid.major = element_blank(), 
                          axis.text=element_text(size=17),                          #axis numbers
                          plot.title=element_text(size=30, hjust=-0.35),            #Main title      --> to center title, hjust=0.5, to left align -0.2
                          legend.title = element_text(size = 16, 
                                                      colour = "black", angle = 0),           #Legend title
                          legend.text = element_text(size = 14, 
                                                     colour = "black", angle = 0),            #Legend text
                          strip.text.x = element_text(size = 20, 
                                                      colour = "black", angle = 0, vjust = 1),  #Facet x text size
                          strip.text.y = element_text(size = 20, 
                                                      colour = "black", angle = 270),
                          strip.background = element_blank(),
                          panel.border = element_blank())   


# Load Water Chemistry Data -----------------------------------------------
# Make a mastersheet with chem and methylation ----------------------------


water_chem <- read_csv("CHEM_DATA/COMPLETE_MASTER_SHEET_ALL_WATER_CHEM.csv",
                       local = locale(encoding = "latin1")) #data sheet contains "lake_name" as first column and chem data as subsequent columns

suva <- read_csv("CHEM_DATA/SUVA_Analysis.csv", local = locale(encoding = "latin1")) #add SUVA data from seperate analysis

years <- read_csv("INFO_SHEETS/Sampling_Year_Per_Lake.csv", locale = locale(encoding = "latin1")) #info sheet with sampling year for each lake

methylation_table <- read_csv("methylation_table.csv", local = locale(encoding = "latin1")) #dataframe made in previous script

methylation_table <- methylation_table %>% merge(.,years, by = "lake_name", all=T) #add the years to methylation table

mastersheet <- left_join(water_chem, methylation_table, by = c("lake_name","Year")) %>% 
                left_join(., suva, by = "lake_name", all=T) %>% mutate(., proton_M = 10^(-pH)) %>% 
                filter(., !lake_name %in% c("VEE")) #merge all methylation and chem data together, take out VEE because the incubation experiment was not succesful

#write.csv(mastersheet, file = "SHEETS_FROM_EACH_SCRIPT/04_mastersheet.csv")

# Subset only the rates from the matersheet -------------------------------

subset_1 <- subset(mastersheet, select = c("lake_name", "Percent_Hg", "Nitrogen_Total", "DOC", "pH", "Sulphate", 
                                           "Iron_Total", "Total_Phophorus_UO", "Total_Arsenic_UO",
                                           "km_24", "kd_24", "SUVA_percent", "proton_M")) %>%
  na.omit(.) %>% column_to_rownames(., var="lake_name")



# Stepwise regression with Km_24 as dependent -----------------------------

# I think a stepwise regression with Km and other varaibles would  --------


scatterplotMatrix(~log(km_24) + log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + #start with scatterplot matrix
                    log(proton_M) + log(DOC) + log(Iron_Total) + 
                    Total_Phophorus_UO + log(Percent_Hg), reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = subset_1) #decided to make sure that raw data was normally distributed before anlaysis

subset_1 <- subset_1 %>% mutate(., logkm = log(km_24))

full.model <- lm(logkm ~ log(Sulphate+1) + log(Percent_Hg) + SUVA_percent + 
                   log(Total_Arsenic_UO) + log(proton_M) + log(DOC) + log(Iron_Total) + 
                   Total_Phophorus_UO, data = subset_1)
step.model <- stepAIC(full.model, direction = "both") #DOC 
summary(step.model)
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=subset_1) #passed
resettest(step.model, power = 2:3, type = "regressor", data = subset_1) #passed
shapiro.test(step.model$residuals) #passed

par(mfrow=c(1,2))
qqnorm(step.model$resid)
qqline(step.model$resid)
hist(step.model$resid,freq=F)
xseq=seq(-5,5,length=200)
lines(xseq,dnorm(xseq,sd=summary(step.model)$sigma))

########DOC is the best predictor of pooled km_24###############

# Show all the variables as dependent with Km_24 --------------------------

sul <- ggplot(data = mastersheet, aes(x=Sulphate, y=km_24, label = lake_name)) + 
  geom_point(aes(size=4), show.legend = FALSE) +
  #geom_text_repel() +
  labs(x = "Sulphate (mg/L)", y = 'Methylation Rate Constant (h^-1)') +
  theme_incubation

#ggsave("sul_Km.svg", path = "Figures/", plot = sul)

doc <- ggplot(data = mastersheet, aes(x=(DOC), y=(km_24))) + 
  labs(x = "DOC (mg/L)", y = 'Methylation Rate Constant (h^-1)') +
  geom_smooth(formula = 'y ~ log(x)', method = 'glm', 
              method.args = list(family = gaussian(link = 'log')), 
              fill = "#515c69", color = "#515c69") +
  geom_point(aes(size=4), show.legend = FALSE) +
  xlim (0, 70) +
  annotate('text', x = 60, y = 1, label = "p-value<0.05~~R^{2}==0.76" ,parse = TRUE,size=4) +
  theme_incubation

#ggsave("DOC_Km.svg", path = "Figures/", plot = doc)

ph <- ggplot(data = mastersheet, aes(x=proton_M, y=log(km_24))) + 
  geom_point(aes(size=4), show.legend = FALSE) +
  labs(x = "pH", y = 'Methylation Rate Constant (h^-1)') +
  theme_incubation

iron <-  ggplot(data = mastersheet, aes(x=log(Iron_Total), y=log(km_24))) + 
  geom_point(aes(size=4), show.legend = FALSE) +
  labs(x = "Iron (mg/L)", y = 'Methylation Rate Constant (h^-1)') +
  theme_incubation

suva.graph <-  ggplot(data = mastersheet, aes(x=SUVA_percent, y=log(km_24))) + 
  geom_point(aes(size=4), show.legend = FALSE) +
  labs(x = "SUVA", y = 'Methylation Rate Constant (h^-1)') +
  theme_incubation


arsenic <-  ggplot(data = mastersheet, aes(x=log(Total_Arsenic_UO), y=log(km_24))) + 
  geom_point(aes(size=4), show.legend = FALSE) +
  labs(x = "Arsenic (ng/L)", y = 'Methylation Rate Constant (h^-1)') +
  theme_incubation


phos <- ggplot(data = mastersheet, aes(x=log(Total_Phophorus_UO+1), y=log(km_24))) + 
  geom_point(aes(size=4), show.legend = FALSE) +
  labs(x = "Phosphorous", y = 'Methylation Rate Constant (h^-1)') +
  theme_incubation


allgraph <- ggarrange(sul, doc, ph, suva.graph, iron, phos, arsenic, nrow = 3, ncol = 3)

ggsave(allgraph, file = "methylation_rate_with_varaibles.svg", dpi = 600, path = "Figures/")


# Show all the variables as dependent with Kd_24 --------------------------
# stats for Kd is in other script -----------------------------------------

sul <- ggplot(data = mastersheet, aes(x=Sulphate, y=kd_24)) + 
  geom_point(aes(size=1), show.legend = FALSE) +
  labs(x = "Sulphate (mg/L)", y = 'Methylation Rate Constant (h^-1)') +
  theme_incubation

doc <- ggplot(data = mastersheet, aes(x=DOC, y=(kd_24))) + 
  geom_point(aes(size=1), show.legend = FALSE) +
  labs(x = "DOC (mg/L)", y = 'Methylation Rate Constant (h^-1)') +
  theme_incubation

ph <- ggplot(data = mastersheet, aes(x=pH, y=kd_24)) + 
  geom_point(aes(size=1), show.legend = FALSE) +
  labs(x = "pH", y = 'Methylation Rate Constant (h^-1)') +
  theme_incubation

iron <-  ggplot(data = mastersheet, aes(x=Iron_Total, y=kd_24)) + 
  geom_point(aes(size=1), show.legend = FALSE) +
  labs(x = "Iron (mg/L)", y = 'Methylation Rate Constant (h^-1)') +
  theme_incubation

suva.graph <-  ggplot(data = mastersheet, aes(x=SUVA_percent, y=kd_24)) + 
  geom_point(aes(size=1), show.legend = FALSE) +
  labs(x = "SUVA", y = 'Methylation Rate Constant (h^-1)') +
  theme_incubation


arsenic <-  ggplot(data = mastersheet, aes(x=Total_Arsenic_UO, y=kd_24)) + 
  geom_point(aes(size=1), show.legend = FALSE) +
  labs(x = "Arsenic (ng/L)", y = 'Methylation Rate Constant (h^-1)') +
  theme_incubation


phos <- ggplot(data = mastersheet, aes(x=Total_Phophorus_UO, y=kd_24)) + 
  geom_point(aes(size=1), show.legend = FALSE) +
  labs(x = "Phosphorous", y = 'Methylation Rate Constant (h^-1)') +
  theme_incubation


allgraph <- ggarrange(sul, doc, ph, suva.graph, iron, phos, arsenic, nrow = 3, ncol = 3)

ggsave(allgraph, file = "demethylation_rate_with_varaibles.svg", dpi = 600, path = "Figures/")



################################### ONLY CHEM DATA ################################################


# check which enviro best describe changes in MeHg percent ----------------
scatterplotMatrix(~log(Percent_Hg) + log(Sulphate+1) + SUVA_percent + log(Total_Arsenic_UO) + 
                    log(proton_M) + log(DOC) + log(Iron_Total) + 
                    log(Total_Phophorus_UO+1), reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = mastersheet) #decided to make sure they were normally distributed before anlaysis


mastersheet <- mastersheet %>% mutate(., logper = log(Percent_Hg))
#remove pocket since it doesn't have SUVA
mastersheet.np <- subset(mastersheet, lake_name!="POCKET")


full.model <- lm(logper ~ log(Sulphate+1) + SUVA_percent + log(Total_Arsenic_UO) + 
                   log(proton_M) + log(DOC) + log(Iron_Total) + 
                   log(Total_Phophorus_UO+1), data = mastersheet.np)
step.model <- stepAIC(full.model, direction = "both") #none
summary(step.model)
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=mastersheet)
resettest(step.model, power = 2:3, type = "regressor", data = subset_1)
shapiro.test(step.model$residuals) 



#make a more simple model with only sul
model.sul <- lm(logper ~ log(Sulphate+1), data = mastersheet.np)
summary(model.sul)
bptest(model.sul,varformula = ~fitted.values(model.sul), studentize=T, data=mastersheet.np)#passed
resettest(model.sul, power = 2:3, type = "regressor", data = mastersheet.np)#passed
shapiro.test(model.sul$residuals) #passed

anova(step.model, model.sul) #p-value = 0.04986

#make a more simple model with only sul and ph
model.sulph <- lm(logper ~ log(Sulphate+1) + log(proton_M), data = mastersheet.np)
summary(model.sulph)
bptest(model.sulph,varformula = ~fitted.values(model.sulph), studentize=T, data=mastersheet.np)#passed
resettest(model.sulph, power = 2:3, type = "regressor", data = mastersheet.np)#passed
shapiro.test(model.sulph$residuals) #passed

anova(model.sulph, model.sul) # p-value = 0.0456

#make a more simple model with only sul and doc
model.suldoc <- lm(logper ~ log(Sulphate+1) + log(DOC), data = mastersheet.np)
summary(model.suldoc)
bptest(model.suldoc,varformula = ~fitted.values(model.suldoc), studentize=T, data=mastersheet.np)#passed
resettest(model.suldoc, power = 2:3, type = "regressor", data = mastersheet.np)#passed
shapiro.test(model.suldoc$residuals) #passed

anova(model.suldoc, model.sul) #p-value = 0.2911

#make a more simple model with only sul and doc
model.phdoc <- lm(logper ~ log(proton_M) + log(DOC), data = mastersheet.np)
summary(model.phdoc)
bptest(model.phdoc,varformula = ~fitted.values(model.phdoc), studentize=T, data=mastersheet.np)#passed
resettest(model.phdoc, power = 2:3, type = "regressor", data = mastersheet.np)#doesnt pass
shapiro.test(model.phdoc$residuals) #passed

anova(model.phdoc, model.sul) #p-value = 3.1469


# Make a figure with percent Hg and sulphate ------------------------------

plot.sul <- ggplot(data = mastersheet, aes(x=log(Sulphate+1), y=log(Percent_Hg))) +
  geom_smooth(data=mastersheet,
              colour="#C6CDF7",method="lm", formula = y~x, se=T, show.legend = F) +
  geom_point(size=5, colour = "#7294D4") +
  labs(x = "ln(Sulfate in Water) (mg/L)", y = 'ln(Ratio of Methylmercury)') +
  expand_limits(x = 0, y = 0) +
  annotate('text', x = 3, y = 1,  label = "p-value<=0.05~~R^{2}==0.39" ,parse = TRUE,size=4) +
  theme_incubation
  

plot.sul


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################



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

# low km sample group treatments ------------------------------------------
# Stepwise regression with Km_24 as dependent -----------------------------

#start with scatterplot matrix
scatterplotMatrix(~sqrt(km_24) + log(Sulphate+1) + log(SUVA_percent) + sqrt(Total_Arsenic_UO) + 
                    sqrt(proton_M) + sqrt(DOC) + sqrt(Iron_Total) + 
                    sqrt(Total_Phophorus_UO) + log(Percent_Hg), reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = low.km) #decided to make sure they were normally distributed before anlaysis


subset.low <- low.km %>% mutate(., sqrtkm = sqrt(km_24))

full.model <- lm(sqrtkm ~ log(Sulphate+1) + log(Percent_Hg) + SUVA_percent + log(Total_Arsenic_UO) + 
                   sqrt(proton_M) + log(DOC) + 
                   Total_Phophorus_UO, data = subset.low) #took off total iron, do not have enough DFs to run all the independent variables
step.model <- stepAIC(full.model, direction = "both") #DOC, percent Hg
summary(step.model)
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=low.km)
resettest(step.model, power = 2:3, type = "regressor", data = low.km)
shapiro.test(step.model$residuals) #passed

par(mfrow=c(1,2))
qqnorm(step.model$resid)
qqline(step.model$resid)
hist(step.model$resid,freq=F)
xseq=seq(-5,5,length=200)
lines(xseq,dnorm(xseq,sd=summary(step.model)$sigma))


# make a reduced model with doc and percent Hg ----------------------------

subset.low <- subset.low %>% mutate(., logDOC = log(DOC)) %>% 
  mutate(., logHg = log(Percent_Hg))

model.doc.hg <- lm(sqrtkm ~ logDOC + logHg, data = subset.low) #there is no sig interaction
summary(model.doc.hg)
bptest(model.doc.hg,varformula = ~fitted.values(model.doc.hg), studentize=T, data=low.km)
resettest(model.doc.hg, power = 2:3, type = "regressor", data = low.km)
shapiro.test(model.doc.hg$residuals) #passed

anova(step.model, model.doc.hg) #there is no significant difference between the two models, will choose the more simple one

crPlots(model.doc.hg)


# Make a graph of the model -----------------------------------------------

doc <- ggplot(data = subset.low, aes(x=log(DOC), y=sqrtkm, label = lake_name)) +
  geom_point(aes(size=1), show.legend = FALSE) +
  labs(x = "log(DOC) (mg/L)", y = 'Methylation Rate Constant (h^-1)') +
  stat_smooth(method = 'lm', formula = y ~ x, se = T, aes(colour = "darkgrey")) +
  annotate('text', x = 2.4, y = 0.62, label = "p-value<0.05~~R^{2}==0.76" ,parse = TRUE,size=4) +
  theme_incubation
doc

hg <- ggplot(data = subset.low, aes(x=log(Percent_Hg), y=sqrtkm, label = lake_name)) +
  geom_point(aes(size=1), show.legend = FALSE) +
  labs(x = "log(MeHg Percent) (%)", y = 'Methylation Rate Constant (h^-1)') +
  stat_smooth(method = 'lm', formula = y ~ x, se = T, aes(colour = "darkgrey")) +
  theme_incubation
hg




model.graph <- ggPredict(model.doc.hg, interactive = T) 

model.graph

# high km sample group treatments ------------------------------------------
# Stepwise regression with Km_24 as dependent -----------------------------

#start with scatterplot matrix
scatterplotMatrix(~log(km_24) + log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                    log(proton_M) + sqrt(DOC) + sqrt(Iron_Total) + 
                    sqrt(Total_Phophorus_UO) + log(Percent_Hg), reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = high.km) #decided to make sure they were normally distributed before anlaysis


subset.high <- high.km %>% mutate(., logkm = log(km_24))

full.model <- lm(logkm ~ log(Sulphate+1) + log(Percent_Hg) + SUVA_percent + log(Total_Arsenic_UO) + 
                   sqrt(proton_M) + log(DOC),
                    data = subset.high) #took off total iron, and TP, do not have enough DFs to run all the independent variables
step.model <- stepAIC(full.model, direction = "both") #DOC, percent Hg
summary(step.model)
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=high.km)
resettest(step.model, power = 2:3, type = "regressor", data = high.km)
shapiro.test(step.model$residuals) #passed

par(mfrow=c(1,2))
qqnorm(step.model$resid)
qqline(step.model$resid)
hist(step.model$resid,freq=F)
xseq=seq(-5,5,length=200)
lines(xseq,dnorm(xseq,sd=summary(step.model)$sigma))


# make a reduced model with only sulphate  ----------------------------

subset.high <- subset.high %>% mutate(., logSulphate = log(Sulphate + 1))

model.sul <- lm(logkm ~ logSulphate, data = subset.high)
summary(model.sul)
bptest(model.sul,varformula = ~fitted.values(model.sul), studentize=T, data=high.km)
resettest(model.sul, power = 2:3, type = "regressor", data = high.km)
shapiro.test(model.sul$residuals) #passed

anova(step.model, model.sul) #there is no significant difference between the two models, will choose the more simple one

crPlots(model.sul)

sul <- ggplot(data = high.km, aes(x=log(Sulphate+1), y=log(km_24), label = lake_name)) +
geom_point(aes(size=1), show.legend = FALSE) +
labs(x = "log(Sulphate + 1 ) (mg/L)", y = 'Methylation Rate Constant (h^-1)') +
  stat_smooth(method = 'lm', formula = y ~ x, se = T, aes(colour = "darkgrey")) +
theme_incubation
sul

model <- ggPredict(model.sul, interactive = T) +
  labs(x = "Sulphate (mg/L)", y = 'Methylation Rate Constant (h^-1)')


########================================================== DONE ================================#################
