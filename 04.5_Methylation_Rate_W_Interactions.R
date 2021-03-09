###############################################################################################
#####THIS R SCRIPT IS FOR MERGING METHYLATION AND CHEM DATA TO ANALYSE#########################
#####We will use data made in the last three scripts and a csv file with chem data#############
################# Regression of Methylation rate with interactions ############################



# Load libraries ----------------------------------------------------------

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



# Load Water Chemistry Data -----------------------------------------------
# Make a mastersheet with chem and methylation ----------------------------


water_chem <- read_csv("CHEM_DATA/COMPLETE_MASTER_SHEET_ALL_WATER_CHEM.csv",local = locale(encoding = "latin1"))

suva <- read_csv("CHEM_DATA/SUVA_Analysis.csv", local = locale(encoding = "latin1"))

years <- read_csv("INFO_SHEETS/Sampling_Year_Per_Lake.csv", locale = locale(encoding = "latin1"))

methylation_table <- read_csv("methylation_table.csv", local = locale(encoding = "latin1"))

methylation_table <- methylation_table %>% merge(.,years, by = "lake_name", all=T)

mastersheet <- left_join(water_chem, methylation_table, by = c("lake_name","Year")) %>% 
  left_join(., suva, by = "lake_name", all=T) %>% mutate(., proton_M = 10^(-pH)) %>% 
  drop_na(km_24)

#write.csv(mastersheet, file = "SHEETS_FROM_EACH_SCRIPT/04_mastersheet.csv")

# Subset only the rates from the matersheet -------------------------------

subset_1 <- subset(mastersheet, select = c("lake_name", "Percent_Hg", "Nitrogen_Total", "DOC", "pH", "Sulphate", 
                                           "Iron_Total", "Total_Phophorus_UO", "Total_Arsenic_UO",
                                           "km_24", "kd_24", "SUVA_percent", "proton_M")) %>%
  na.omit(.) %>% column_to_rownames(., var="lake_name")



# Stepwise regression with Km_24 as dependent -----------------------------

# I think a stepwise regression with Km and other varaibles would  --------

# test with ph interaction ------------------------------------------------

#start with scatterplot matrix
scatterplotMatrix(~log(km_24) + log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                    log(proton_M) + log(DOC) + log(Iron_Total) + 
                    Total_Phophorus_UO + log(Percent_Hg), reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = subset_1) #decided to make sure they were normally distributed before anlaysis

subset_1 <- subset_1 %>% mutate(., logkm = log(km_24))

full.model <- lm(logkm ~ log(Sulphate+1) + log(Percent_Hg) + SUVA_percent + 
                   log(Total_Arsenic_UO) + pH + log(DOC) + log(Iron_Total) + 
                   Total_Phophorus_UO +
                   pH:log(Sulphate+1) +
                   pH:log(DOC) +
                   pH:log(Iron_Total), data = subset_1)
step.model <- stepAIC(full.model, direction = "both") #DOC is again the only variable
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





# test with sulfate interaction -------------------------------------------

#start with scatterplot matrix
scatterplotMatrix(~log(km_24) + log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                    log(proton_M) + log(DOC) + log(Iron_Total) + 
                    Total_Phophorus_UO + log(Percent_Hg), reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = subset_1) #decided to make sure they were normally distributed before anlaysis

subset_1 <- subset_1 %>% mutate(., logkm = log(km_24))

full.model <- lm(logkm ~ log(Sulphate+1) + log(Percent_Hg) + SUVA_percent + 
                   log(Total_Arsenic_UO) + pH + log(DOC) + log(Iron_Total) + 
                   Total_Phophorus_UO +
                   pH:log(Sulphate+1) +
                   log(Sulphate+1):log(DOC) +
                   log(Sulphate+1):log(Iron_Total), data = subset_1)
step.model <- stepAIC(full.model, direction = "both") #no significant
summary(step.model)
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=subset_1)
resettest(step.model, power = 2:3, type = "regressor", data = subset_1)
shapiro.test(step.model$residuals) #passed


# test with sulfate interaction -------------------------------------------

#start with scatterplot matrix
scatterplotMatrix(~log(km_24) + log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                    log(proton_M) + log(DOC) + log(Iron_Total) + 
                    Total_Phophorus_UO + log(Percent_Hg), reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = subset_1) #decided to make sure they were normally distributed before anlaysis

subset_1 <- subset_1 %>% mutate(., logkm = log(km_24))

full.model <- lm(logkm ~ log(Sulphate+1) + log(Percent_Hg) + SUVA_percent + 
                   log(Total_Arsenic_UO) + pH + log(DOC) + log(Iron_Total) + 
                   Total_Phophorus_UO +
                   pH:log(DOC) +
                   log(Sulphate+1):log(DOC) +
                   log(DOC):log(Iron_Total), data = subset_1)
step.model <- stepAIC(full.model, direction = "both") #no significant
summary(step.model)
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=subset_1)
resettest(step.model, power = 2:3, type = "regressor", data = subset_1)
shapiro.test(step.model$residuals) #passed



###################################################################################################
###################################################################################################
# Seperate the sampels to the two different groups ------------------------

low.km <- mastersheet %>% filter(lake_name %in% c("RAT", "FRAME", "BCR07",
                                                  "DAVID", "NIVEN", "YKE1",
                                                  "VEE", "BC20", "YK42", "BC21",
                                                  "PROSPEROUS", "YK60",
                                                  "ICING", "BC18"))


high.km <- mastersheet %>% filter(!lake_name %in% c("RAT", "FRAME", "BCR07",
                                                    "DAVID", "NIVEN", "YKE1",
                                                    "VEE", "BC20", "YK42", "BC21",
                                                    "PROSPEROUS", "YK60",
                                                    "ICING", "BC18"))



# low km sample group treatments ------------------------------------------
# Stepwise regression with Km_24 as dependent -----------------------------

#start with scatterplot matrix
scatterplotMatrix(~sqrt(km_24) + log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                    proton_M + sqrt(DOC) + sqrt(Iron_Total) + 
                    log(Total_Phophorus_UO+1) + log(Percent_Hg),  reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = low.km) #decided to make sure they were normally distributed before anlaysis


subset.low <- low.km %>% mutate(., sqrtkm = sqrt(km_24))

full.model <- lm(sqrtkm ~ log(Sulphate) + log(Percent_Hg) + 
                    proton_M + log(DOC), data = subset.low) #took off total iron, do not have enough DFs to run all the independent variables
step.model <- stepAIC(full.model, direction = "both") 
summary(step.model) #DOC, percent Hg
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=low.km)
resettest(step.model, power = 2:3, type = "regressor", data = low.km) #too many factors
shapiro.test(step.model$residuals) #passed


par(mfrow=c(1,2))
qqnorm(step.model$resid)
qqline(step.model$resid)
hist(step.model$resid,freq=F)
xseq=seq(-5,5,length=200)
lines(xseq,dnorm(xseq,sd=summary(step.model)$sigma))

#make a more simple model with only sul and doc
model.nodoc <- lm(sqrtkm ~ log(Percent_Hg) + log(DOC),  data = subset.low)
summary(model.nodoc)
bptest(model.nodoc,varformula = ~fitted.values(model.nodoc), studentize=T, data=subset.low)#passed
resettest(model.nodoc, power = 2:3, type = "regressor", data = subset.low)#doesnt pass
shapiro.test(model.nodoc$residuals) #passed

anova(model.nodoc, step.model) #the simpler model is not significantly different than the model above, so we will use it



# high km sample group treatments ------------------------------------------
# Stepwise regression with Km_24 as dependent -----------------------------

#start with scatterplot matrix
scatterplotMatrix(~sqrt(km_24) + log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                    log(proton_M) + sqrt(DOC) + sqrt(Iron_Total) + 
                    log(Total_Phophorus_UO+1) + log(Percent_Hg), reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = high.km) #decided to make sure they were normally distributed before anlaysis


subset.high <- high.km %>% mutate(., sqrtkm = sqrt(km_24))

full.model <- lm(sqrtkm ~ log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                   log(proton_M) + sqrt(DOC) + log(Percent_Hg), data = subset.high) #took off total iron and phosphorous, do not have enough DFs to run all the independent variables
step.model <- stepAIC(full.model, direction = "both") #sulfate and suva
summary(step.model)
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=high.km)
resettest(step.model, power = 2:3, type = "regressor", data = high.km) #too many factors
shapiro.test(step.model$residuals) #passed


par(mfrow=c(1,2))
qqnorm(step.model$resid)
qqline(step.model$resid)
hist(step.model$resid,freq=F)
xseq=seq(-5,5,length=200)
lines(xseq,dnorm(xseq,sd=summary(step.model)$sigma))

#make a more simple model with only sul and doc
model.nodoc <- lm(sqrtkm ~ log(Sulphate + 1),  data = subset.high)
summary(model.nodoc)
bptest(model.nodoc,varformula = ~fitted.values(model.nodoc), studentize=T, data=subset.high)#passed
resettest(model.nodoc, power = 2:3, type = "regressor", data = subset.high)#doesnt pass
shapiro.test(model.nodoc$residuals) #passed

anova(model.nodoc, step.model) #the simpler model is not significantly different than the model above, so we will use it


#######################################t-tests################################################

# T-tests with the two groups of Km ---------------------------------------

# Seperate the data into two categories -----------------------------------

group.1 <- mastersheet %>% subset(., lake_name %in% c("RAT", "FRAME", "BCR07",
                                                   "DAVID", "NIVEN", "YKE1",
                                                   "VEE", "BC20", "YK42", "BC21",
                                                   "PROSPEROUS", "YK60",
                                                    "ICING", "BC18")) %>% add_column(grouping = "ONE")

group.2 <- mastersheet %>% filter(!lake_name %in% c("RAT", "FRAME", "BCR07",
                                                 "DAVID", "NIVEN", "YKE1",
                                                 "VEE", "BC20", "YK42", "BC21",
                                                 "PROSPEROUS", "YK60",
                                                 "ICING", "BC18")) %>% add_column(grouping = "TWO")

categorical <- rbind(group.1, group.2)




#sulfate
categorical <- categorical %>% mutate(., logsul = log(Sulphate+1))


model <- lm(logsul ~ grouping, data = categorical)
summary(model) #significant
leveneTest(logsul ~ grouping, data = categorical) #good
shapiro.test(model$residuals) 

#SUVA
categorical <- categorical %>% mutate(., logSUVA = log(SUVA_percent))

model.zero <- lm(logSUVA ~ grouping, data = categorical)
summary(model.zero) #no
leveneTest(logSUVA ~ grouping, data = categorical) #good
shapiro.test(model.zero$residuals) 

#Arsenic
model.zero <- lm(log(Total_Arsenic_UO) ~ grouping, data = categorical)
summary(model.zero) #no
leveneTest(log(Total_Arsenic_UO) ~ grouping, data = categorical) #good
shapiro.test(model.zero$residuals) #meh

#ph
model.zero <- lm(sqrt(proton_M) ~ grouping, data = categorical)
summary(model.zero) #no
leveneTest(sqrt(proton_M) ~ grouping, data = categorical) #good
shapiro.test(model.zero$residuals) 

#DOC
model.zero <- lm(DOC ~ grouping, data = categorical)
summary(model.zero)#no
leveneTest(DOC ~ grouping, data = categorical) #good
shapiro.test(model.zero$residuals) 

#Iron
model.zero <- lm(log(Iron_Total) ~ grouping, data = categorical)
summary(model.zero)
leveneTest(log(Iron_Total) ~ grouping, data = categorical) #not good
shapiro.test(model.zero$residuals)

#cannot pass assumptions, will do non-perametric test
wilcox.test((Iron_Total) ~ grouping, data = categorical) #not significant
stripchart((Iron_Total) ~ grouping, data = categorical, meth="stack")


#Phos
model.zero <- lm(log(Total_Phophorus_UO+1) ~ grouping, data = categorical)
summary(model.zero) #not significant
leveneTest(log(Total_Phophorus_UO+1) ~ grouping, data = categorical) #good
shapiro.test(model.zero$residuals) 




























#########################################################################################################
#########################################################################################################
#########################################################################################################

# k-means clustering of data ----------------------------------------------


# subset only km24 data and sluphate data ---------------------------------

kdata <- subset_1


# Subset into two groups --------------------------------------------------
set.seed(25)
clusters <- kmeans(kdata[,c("Sulphate","km_24")], 2)

# Save the cluster number in the dataset as column 'Borough'
kdata$group <- as.factor(clusters$cluster)

# Inspect 'clusters'
str(clusters)

#seperate the dataset into two using kmeans clustering

low <- kdata %>% dplyr::filter(.,group==1)

high <- kdata %>% dplyr::filter(group == 2)





# stepwise regession for low ----------------------------------------------

# Stepwise regression with Km_24 as dependent -----------------------------

#start with scatterplot matrix
scatterplotMatrix(~log(km_24) + log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                    log(proton_M) + sqrt(DOC) + sqrt(Iron_Total) + 
                    sqrt(Total_Phophorus_UO) + log(Percent_Hg), reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = low) #decided to make sure they were normally distributed before anlaysis


subset.low <- low %>% mutate(., logkm = log(km_24))

full.model <- lm(logkm ~ log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                   log(proton_M) + log(DOC) + sqrt(Iron_Total) + 
                   sqrt(Total_Phophorus_UO), data = subset.low) #took off total iron, do not have enough DFs to run all the independent variables
step.model <- stepAIC(full.model, direction = "both") #DOC, %MeHg
summary(step.model)
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=low.km)
resettest(step.model, power = 2:3, type = "regressor", data = low.km) 
shapiro.test(step.model$residuals) #passed



#make a simpler model with just DOC
model.doc <- lm(logkm ~ log(DOC), data = subset.low) # DOC is significant
summary(model.doc)
bptest(model.doc,varformula = ~fitted.values(model.doc), studentize=T, data=low.km)
resettest(model.doc, power = 2:3, type = "regressor", data = low.km) 
shapiro.test(model.doc$residuals) #passed

anova(step.model, model.doc) #not sig different


doc <- ggplot(data = low.km, aes(x=log(DOC), y=log(km_24), label = lake_name)) + 
  geom_point(aes(size=4), show.legend = FALSE) +
  geom_text_repel() +
  labs(x = "ln(DOC) (mg/L)", y = 'ln(Methylation Rate Constant) (h^-1)') +
  stat_smooth(method = 'lm', formula = y ~ x, se = T, aes(colour = "darkgrey")) +
  theme_incubation

doc


# stepwise regession for high ----------------------------------------------

# Stepwise regression with Km_24 as dependent -----------------------------

#start with scatterplot matrix
scatterplotMatrix(~km_24 + log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                    log(proton_M) + sqrt(DOC) + sqrt(Iron_Total) + 
                    sqrt(Total_Phophorus_UO) + log(Percent_Hg), reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = high) #decided to make sure they were normally distributed before anlaysis


subset.high <- high %>% mutate(., logkm = log(km_24))

full.model <- lm(logkm ~ log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                   log(proton_M) + sqrt(DOC) + log(Iron_Total) + 
                   sqrt(Total_Phophorus_UO), data = subset.high) #took off total iron, do not have enough DFs to run all the independent variables
step.model <- stepAIC(full.model, direction = "both") #DOC, %MeHg
summary(step.model)
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=low.km)
resettest(step.model, power = 2:3, type = "regressor", data = low.km) 
shapiro.test(step.model$residuals) #passed





# Try simple regression with certain terms --------------------------------


# sulphate ----------------------------------------------------------------

model.sul <- lm(logkm ~ log(Sulphate + 1), data = subset.high)
summary(model.sul) #not significant


# DOC     ----------------------------------------------------------------

model.doc <- lm(logkm ~ log(DOC), data = subset.high)
summary(model.doc) #not significant



# Iron    ----------------------------------------------------------------

model.doc <- lm(logkm ~ log(Iron_Total), data = subset.high)
summary(model.doc) #not significant



# Subset into three groups --------------------------------------------------
set.seed(25)
clusters <- kmeans(kdata[,c("Sulphate","km_24")], 3)

# Save the cluster number in the dataset as column 'Borough'
kdata$group <- as.factor(clusters$cluster)

# Inspect 'clusters'
str(clusters)


#seperate the dataset into two using kmeans clustering

low <- kdata %>% dplyr::filter(.,group==2)

mid <- kdata %>% dplyr::filter(group == 1)

high <- kdata %>% dplyr::filter(group == 3)



# stepwise regession for low ----------------------------------------------

# Stepwise regression with Km_24 as dependent -----------------------------

#start with scatterplot matrix
scatterplotMatrix(~km_24 + log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                    log(proton_M) + sqrt(DOC) + sqrt(Iron_Total) + 
                    sqrt(Total_Phophorus_UO) + log(Percent_Hg), reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = low) #decided to make sure they were normally distributed before anlaysis


subset.low <- low %>% mutate(., logkm = log(km_24))

full.model <- lm(km_24 ~ log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                   log(proton_M) + log(DOC) + sqrt(Iron_Total) + 
                   sqrt(Total_Phophorus_UO), data = subset.low) #took off total iron, do not have enough DFs to run all the independent variables
step.model <- stepAIC(full.model, direction = "both") #no significant factor
summary(step.model)
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=low.km)
resettest(step.model, power = 2:3, type = "regressor", data = low.km) 
shapiro.test(step.model$residuals) #passed


# stepwise regession for mid ----------------------------------------------

# Stepwise regression with Km_24 as dependent -----------------------------

#start with scatterplot matrix
scatterplotMatrix(~log(km_24) + log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                    proton_M + sqrt(DOC) + sqrt(Iron_Total) + 
                    sqrt(Total_Phophorus_UO) + log(Percent_Hg), reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = mid) #decided to make sure they were normally distributed before anlaysis


subset.mid <- mid %>% mutate(., logkm = log(km_24))

full.model <- lm(logkm ~ log(Sulphate+1) +
                   log(proton_M) + log(DOC) + sqrt(Iron_Total), data = subset.mid) #took off total iron, do not have enough DFs to run all the independent variables
step.model <- stepAIC(full.model, direction = "both") #no significant factor
summary(step.model)
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=subset.mid)
resettest(step.model, power = 2:3, type = "regressor", data = mid.km)  #not enough points
shapiro.test(step.model$residuals) #passed


# stepwise regession for high ----------------------------------------------

# Stepwise regression with Km_24 as dependent -----------------------------

#start with scatterplot matrix
scatterplotMatrix(~log(km_24) + log(Sulphate+1) + log(SUVA_percent) + log(Total_Arsenic_UO) + 
                    proton_M + sqrt(DOC) + sqrt(Iron_Total) + 
                    sqrt(Total_Phophorus_UO) + log(Percent_Hg), reg.line = lm,
                  smooth = TRUE, span = 0.5,
                  diagonal = "denisty", data = high) #decided to make sure they were normally distributed before anlaysis


subset.high <- high %>% mutate(., logkm = log(km_24))

full.model <- lm(logkm ~ log(Sulphate+1) +
                   log(proton_M) + log(DOC), data = subset.high) #not enough points to do a stepwise
step.model <- stepAIC(full.model, direction = "both") #no significant factor
summary(step.model)
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=subset.high)
resettest(step.model, power = 2:3, type = "regressor", data = high.km)  #not enough points (n=4)
shapiro.test(step.model$residuals) #passed