###############################################################################################
#######THIS R SCRIPT IS MULTIVARAITE ANLAYSES FOR WATER CHEM AND RATES#########################
#####We will use data made in the last three scripts and the mastersheet#######################
###############################################################################################

# Load Libraries ----------------------------------------------------------
library(tidyverse) # for basic "tidy" workflows and data wrangling
library(janitor) # lots of help when dealing with non-R-friendly datasets
library(lubridate) # for working with dates
library(ggplot2) #plotting
library(ggpubr) #probably don't need it
library(car) #for stats stuff like crPLots
library(lmtest) #for assumptions of reg
library(vegan) #this is for PCA
library(ggrepel)
library(ggfortify)
library(MASS)



######### subset data
# With Nitrogen ------------------------------------------------------------
subset_1 <- subset(mastersheet, select = c("lake_name", "Nitrogen_Total", "DOC", "pH", "Sulphate", 
                                           "Iron_Total", "Total_Phophorus_UO", "Total_Arsenic_UO",
                                           "km_24", "kd_24", "SUVA_percent")) %>%
  na.omit(.) %>% column_to_rownames(., var="lake_name")

# Without Nitrogen ---------------------------------------------------------
subset_2 <- subset(mastersheet, select = c("lake_name", "DOC", "pH", "Sulphate", 
                                           "Iron_Total", "Total_Phophorus_UO", "Total_Arsenic_UO",
                                           "km_24", "kd_24", "SUVA_percent")) %>%
  na.omit(.) %>% column_to_rownames(., var="lake_name")


# Arrange data to work with ggplot2 ---------------------------------------
df <- subset_1
row.names(df) <- paste(df$lake_name, row.names(df), sep="_") 
df$lake_name <- NULL

head(df)

df_pca <- prcomp(df, scale = T)

#ad theme
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background=element_blank(),axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(colour="black"),
             axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

# Add a theme and plot ----------------------------------------------------
figure = autoplot(df_pca, data = subset_1, colour = "lake_name",
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 6) +
  theme + geom_point(aes(size=5, color="green"))

figure
#ggsave(file="test.pdf", plot=figure, dpi=600) #save the data


# It seems that arsenic and sulfate are perpendicualr which seems  --------

plot <- ggplot(data = subset_1, aes(x=Sulphate, y=km_24, 
                                       fill=km_24)) + 
  geom_point(aes(size=4, colour = km_24)) +
  labs(x = "Sulfate in Water (ng/L)", y = 'Total Arsenic in Water (ng/L)') +
  geom_text_repel(data = subset_1, aes(label=rownames(subset_1))) +
  theme_incubation

plot

#try it out with sulfate and Km_24
plot <- ggplot(data = subset_1, aes(x=Sulphate, y=km_24)) + 
  geom_point(aes(size=4, colour = km_24)) +
  labs(x = "Sulfate in Water (ng/L)", y = 'Methylation Rate Constant (/h)') +
  geom_text_repel(data = subset_1, aes(label=rownames(subset_1))) +
  theme_incubation

plot

#Mija is very much and outlier in this system, however it also looks like the ####
#four lakes Niven, David, BCR07, Frame, Rat####
subset_2 <- subset(subset_1, !rownames(subset_1) %in% c("DAVID","RAT","BCR07","NIVEN","FRAME"))



# Plot out now without the 5 lakes ----------------------------------------

plot <- ggplot(data = subset_2, aes(x=Sulphate, y=km_24, 
                                    fill=km_24)) + 
  geom_point(aes(size=4, colour = km_24)) +
  labs(x = "Sulfate in Water (ng/L)", y = 'Methylation Rate Constant (/h)') +
  geom_text_repel(data = subset_2, aes(label=rownames(subset_2))) +
  theme_incubation

plot

##It almost seems like an exponential relationship between rate and sulfat
##However this doesn't seem to work that well...

model <- lm(km_24 ~ I(Sulphate^(0.22)), data = subset_2)
summary(model)
crPlots(model)
bptest(model,varformula = ~fitted.values(model), studentize=T, data=subset_2)
resettest(model, power = 2:3, type = "regressor", data = subset_2)
shapiro.test(model$residuals)


#Make a box cox transformation to meet normality assumptions###
b = boxcox(km_24 ~ Sulphate, data = subset_2)
b
lamda=b$x
lik=b$y
bc=cbind(lamda, lik)
bc[order(-lik),]


# I think a stepwise regression with Km and other varaibles would  --------
full.model <- lm(km_24 ~ Sulphate + SUVA_percent + Total_Arsenic_UO + pH + DOC + Iron_Total + 
              Total_Phophorus_UO, data = subset_2)
step.model <- stepAIC(full.model, direction = "both") #sulfate came out as best predictor
summary(step.model)
bptest(step.model,varformula = ~fitted.values(step.model), studentize=T, data=subset_2)
resettest(step.model, power = 2:3, type = "regressor", data = subset_2)
shapiro.test(step.model$residuals) #did not pass = p=0.068, not great




