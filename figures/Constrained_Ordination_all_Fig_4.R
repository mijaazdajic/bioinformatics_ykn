###############################################################################################
###########THIS R SCRIPT IS FOR ANALYSING THE QIIME2 DATA FROM THE PIPELINE####################
#####The importing is done with phyloseq, then other pkgs can be use for analysis##############
############This is to do contrained ordiation wtih spatial data ##############################


# Load packages -----------------------------------------------------------
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(vegan)
library(dplyr)
library(ape)
library(ggthemes)
library(ggpubr)
library(ggrepel)


#####make a theme###########
theme_incubation <- theme(panel.background = element_rect(fill = "white", linetype = "solid", colour = "black"), 
                          legend.key = element_rect(fill = "white"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
                          axis.title = element_text(size = 20),
                          axis.text=element_text(size=40),                                               #axis numbers
                          plot.title=element_text(size=30, hjust=-0.35),                                 #Main title      --> to center title, hjust=0.5, to left align -0.2
                          legend.title = element_text(size = 16, colour = "black", angle = 0),           #Legend title
                          legend.text = element_text(size = 14, colour = "black", angle = 0),            #Legend text
                          strip.text.x = element_text(size = 50, colour = "black", angle = 0, vjust = 1),        #Facet x text size
                          strip.text.y = element_text(size = 50, colour = "black", angle = 270),
                          strip.background = element_blank(),
                          panel.border = element_blank(),
                          plot.margin = unit(c(1,4,1,1), "lines"))


# load data from previous scripts -----------------------------------------


metadata <- read_tsv("QIIME2/metadata.tsv") %>% dplyr::rename(lake_name="Description")#load metadata table

suva <- read_csv("QIIME2/SUVA_Analysis_for_QIIME.csv", local = locale(encoding = "latin1"))

years <- read_csv("INFO_SHEETS/Sampling_Year_Per_Lake.csv", locale = locale(encoding = "latin1"))

methylation_table <- read_csv("methylation_table.csv", local = locale(encoding = "latin1"))

methylation_table <- methylation_table %>% merge(.,years, by = "lake_name", all=T)
methylation_table$lake_name <- gsub("NIVEN", "NIVEN_2016", methylation_table$lake_name)
methylation_table$lake_name <- gsub("BC20", "BC20_2016", methylation_table$lake_name)


env.data <- left_join(metadata, methylation_table, by = c("lake_name","Year")) %>% 
  left_join(., suva, by = "lake_name", all=T) %>% 
  dplyr::select(.,c('#SampleID', lake_name, MeHg_water_ngl, Year, THg_water_ngl, Percent_Hg, DOC, pH,
                    Sulphate, Iron_Total, Dissolved_Phosphorus_UO, Total_Phophorus_UO, 
                    Dissolved_Arsenic_UO, Total_Arsenic_UO, km_24, kd_24, SUVA, SUVA_percent)) %>% 
  subset(., lake_name != "POCKET") %>% 
  filter(., !lake_name %in% c("NIVEN_2015", "BC20_2015"))



# Load microbial data -----------------------------------------------------
SVs <-  read_qza("QIIME2/14_dada2_table_final.qza") #load the feature table
taxonomy <- read_qza("QIIME2/taxonomy.qza") #load taxonomy
taxtable<-taxonomy$data %>% as_tibble() %>% separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) #convert the table into a tabular split version
taxtable
tree <-  read_qza("QIIME2/6_rooted-tree.qza") #adding the rooted tree made in qiime2


physeq2<-phyloseq(
  otu_table(SVs$data, taxa_are_rows = T), 
  phy_tree(tree$data), 
  tax_table(as.data.frame(taxtable) %>% 
              dplyr::select(-c(Confidence)) %>% column_to_rownames("Feature.ID") %>% as.matrix()), #moving the taxonomy to the way phyloseq wants it
  sample_data(env.data %>% as.data.frame() %>% column_to_rownames("#SampleID"))
)

physeq2



# remove Pocket -----------------------------------------------------------

physeq2 = subset_samples(physeq2, lake_name != "POCKET")



#Remove taxa not seen more than 3 times in at least 50% of the samples. This protects 
#against an OTU with small mean & trivially large C.V.
GP1 = filter_taxa(physeq2, function(x) sum(x > 3) > (0.5*length(x)), TRUE)
GP1 = subset_samples(GP1, lake_name != "BC20_2015") 
GP1 = subset_samples(GP1, lake_name != "NIVEN_2015")


#Keep only the most abundant five phyla.
phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
GP1 = prune_taxa((tax_table(GP1)[, "Phylum"] %in% top5phyla), GP1)

#physeq dataset in which the abundance values have been transformed to relative abundance 
#within each sample
GP2 = transform_sample_counts(GP1, function(x) 1E2 * x/sum(x))


#check again..
p <- plot_bar(GP2,x="lake_name" ,fill="Class")
p <- p + geom_bar(aes(color=Class), stat="identity", position="stack")
p <- p + theme(legend.title=element_blank(), legend.text = element_blank())
p



# Load spatial data -------------------------------------------------------

space <- read_csv("INFO_SHEETS/Spatial_Coordinates.csv", local = locale(encoding = "latin1"))  %>% 
  left_join(., env.data, by = c("lake_name", "Year")) %>% 
  dplyr::select(., c(lake_name, Latitude, Longitude, `#SampleID`)) %>% 
  na.omit(.) %>% filter(., !lake_name %in% c("NIVEN_2015", "BC20_2015", "POCKET")) %>% 
  column_to_rownames("#SampleID") %>% dplyr::select(-c(lake_name))
  

# PCNM calculations -------------------------------------------------------

pcnm <- pcnm(dist(space, method = "euclidian"), dist.ret = T)
pcnm <- as.matrix(pcnm$dist)



space.ord <- prcomp(pcnm, scale = T)
eigs <-  space.ord$sdev^2
eigs[1] / sum(eigs) #first eigen vector explains 0.5654
eigs[2] / sum(eigs) #second eigen vector explains 0.1306

PCspace <- as.data.frame(space.ord$rotation[,1]) #pick out only the first vector since it explains >50% of variability



# Calculate distance based RDA using Weighted UniFrac distance mat --------

# add space vector to data ------------------------------------------------
env.data <-  env.data %>% column_to_rownames(., var = "#SampleID") %>% 
  merge(., PCspace, by = "row.names", all=T) %>% 
  dplyr::select(-"lake_name") %>% dplyr::select(-"Year") %>%  #take out row names that are not num so can be standardized
  column_to_rownames(var = "Row.names") %>%
  rename(space.vector= 'space.ord$rotation[, 1]')

# Transform variables to have normal distribution -------------------------

hist(env.data$DOC) #good approximate normal dist
hist(env.data$pH) #good approximate normal dist
hist(env.data$Sulphate) #skewed to the left
hist(log(env.data$Sulphate + 1)) #better
env.data$lnsul <- log(env.data$Sulphate + 1)
hist(env.data$Iron_Total) #skewed to the left
hist(log(env.data$Iron_Total)) #better
env.data$lnFe <- log(env.data$Iron_Total)
hist(env.data$Total_Phophorus_UO) #skewed to left
hist(sqrt(env.data$Total_Phophorus_UO)) #better
env.data$sqrtP <- sqrt(env.data$Total_Phophorus_UO) 
hist(env.data$km_24) #looks ok
hist(env.data$space.vector)
hist(log10(env.data$space.vector+1)) #better
env.data$space.trans <- log10(env.data$space.vector+1)
hist(env.data$SUVA) #good

# Scale Environmental variables -------------------------------------------

scaled.env <-  env.data %>% 
  decostand(., method = "standardize", na.rm = T) %>% 
  select("space.trans", "SUVA","DOC", "pH", "lnsul", "lnFe", "sqrtP")

# Make a distance matrix using UNIFRAC weighted ---------------------------

wUF.dm<-UniFrac(GP2, weighted = TRUE, fast = TRUE) 
wUF.table<-as.matrix(dist(wUF.dm))


# Calculate distance based RDA using Weighted UniFrac distance mat --------

dbRDA<-capscale(wUF.dm ~ space.trans * DOC * pH * lnsul *
                   lnFe * SUVA * sqrtP, data=scaled.env, sqrt.dist = T)



plot(dbRDA, main = "dbRDA")
summary(dbRDA)
screeplot(dbRDA)
RsquareAdj(dbRDA)##Variance partition


vegan::anova.cca(dbRDA, by="term", step = 200)
vif.cca(dbRDA) # a vif of more than 10 indicates that a variable is strongly dependent on others

# Do a stepwise regression to find variables that best explain species differeces
rda.1 <- capscale(wUF.dm~1, data = scaled.env, sqrt.dist = T)
mod <- ordistep(rda.1, scope = formula(dbRDA),
                direction = "both", permutations = 999, trace = T, sqrt.dist = T) #bidirectional selection
summary(mod)
anova(mod, permutations = 999) #the overall model is significant
anova(mod, by = "terms", permutations = 999) #space is not significant, might slighly change since permuations
anova(mod, by = "axis", permutations = 999) #the CAP3 is not significant

mod2<-capscale(wUF.dm ~ DOC + SUVA + space.trans, data=scaled.env, sqrt.dist = T)
anova(mod2, permutations = 999) #better p-value
anova(mod2, by = "terms", permutations = 999) #DOC and SUVA is significant
anova(mod2, by = "axis", permutations = 999) #CAP1 and CAP2 is significant



summary(mod2)

#check effect of variables
doc <- capscale(wUF.dm ~ Condition(SUVA) + DOC + Condition(space.trans), data=scaled.env, sqrt.dist = T)
summary(doc) #6.84%

suva <- capscale(wUF.dm ~ Condition(DOC) + SUVA + Condition(space.trans), data=scaled.env, sqrt.dist = T)
summary(suva) #6.72%

space <- capscale(wUF.dm ~ Condition(DOC) + Condition(SUVA) + space.trans, data=scaled.env, sqrt.dist = T)
summary(space) #5.86%



plot(mod)
summary(mod)
screeplot(mod)
RsquareAdj(mod)


# Make plot ---------------------------------------------------------------

smry <- summary(mod2)
scrs <- scores(mod2)
df1  <- data.frame(smry$sites[,1:2]) # site scores for RDA1 and RDA2
df1$site <- rownames(df1)  #add site names
df2  <- data.frame(smry$biplot[,1:2])  # mapping environmental variables

name_key <- read_csv("INFO_SHEETS/lake_name_key.csv", local = locale(encoding = "latin1"))
df1 <- df1 %>% left_join(., name_key, by = "site") %>% column_to_rownames(., var = "site")



xtext <- paste("CAP 1 (", round(smry$concont[[1]][2, 1] * 100, 2), " % of fitted, ", 
               round(smry$cont[[1]][2, 1] * 100, 2), " % of total variation)", sep = "")
ytext <- paste("CAP 2 (", round(smry$concont[[1]][2, 2] * 100, 2), " % of fitted, ", 
               round(smry$cont[[1]][2, 2] * 100, 2), " % of total variation)", sep = "")



rda.plot <- ggplot(df1, aes(x=CAP1, y=CAP2, colour = lake_name)) + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  geom_point(size = 4) +
  geom_text_repel(aes(label=lake_name), colour = "grey50",size=5, point.padding = 0.5, 
                  segment.color = 'transparent') +
  geom_segment(data=df2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               color="grey25", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2, aes(x=CAP1,y=CAP2,label=rownames(df2),
                          hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2))), 
            color="grey25", size=6) +
  #xlim(-1.05, 1.65) +
  #ylim(-1.5, 1.05) +
  xlab(xtext) +
  ylab(ytext) + 
  coord_fixed() +
  theme_incubation +
  theme(legend.position = "right",
        legend.title = element_blank())


rda.plot

figure_5 <- rda.plot

ggsave(figure_5, file="Figure_5.svg", width = 10, height=10, units = "in", path = "Figures/")




# There are 7 fractions ---------------------------------------------------


# fractions [DOC+SUVA+space]:
rda.all <- capscale(wUF.dm ~ Condition(DOC) + Condition(space.trans) + SUVA, data = scaled.env, sqrt.dist = T)
# fractions [DOC+SUVA]:
rda.docsuva <- capscale(wUF.dm ~ DOC + SUVA , data = scaled.env, sqrt.dist = T)
# fractions [Spcae+SUVA]:
rda.spacesuva <- capscale(wUF.dm ~ space.trans + SUVA , data = scaled.env, sqrt.dist = T)
# fractions [DOC+Space]:
rda.docspace <- capscale(wUF.dm ~ space.trans + DOC, data = scaled.env, sqrt.dist = T)
# fractions [Spcae]:
rda.space <- capscale(wUF.dm ~ space.trans, data = scaled.env, sqrt.dist = T)
# fractions [SUVA]:
rda.suva <- capscale(wUF.dm ~ SUVA , data = scaled.env, sqrt.dist = T)
# fractions [DOC]:
rda.doc <- capscale(wUF.dm ~ DOC, data = scaled.env, sqrt.dist = T)


# calculate variations explained by individual fractions ------------------

varp <- varpart (wUF.dm, ~ DOC, ~ SUVA, data = scaled.env, sqrt.dist = T)
varp

plot (varp, digits = 2, Xnames = c('DOC', 'SUVA'), bg = c('navy', 'tomato'))


# Do the same with Unweighted UNIFRAC -------------------------------------
# Make a distance matrix using UNIFRAC weighted ---------------------------

uwUF.dm<-UniFrac(GP2, weighted = FALSE, fast = TRUE) 
uwUF.table<-as.matrix(dist(uwUF.dm))


# Calculate distance based RDA using Weighted UniFrac distance mat --------

dbRDA<-capscale(uwUF.dm ~ space.trans + DOC + pH + lnsul +
                  + lnFe, data=scaled.env, sqrt.dist = T)

plot(dbRDA, main = "dbRDA")
summary(dbRDA)
screeplot(dbRDA)
RsquareAdj(dbRDA)##Variance partition


vegan::anova.cca(dbRDA, by="term", step = 200)
vif.cca(dbRDA) # a vif of more than 10 indicates that a variable is strongly dependent on others

# Do a stepwise regression to find variables that best explain species differeces
rda.1 <- capscale(uwUF.dm~1, data = scaled.env, sqrt.dist = T)
mod <- ordistep(rda.1, scope = formula(dbRDA),
                direction = "both", permutations = 999, trace = T, sqrt.dist = T) #bidirectional selection

mod
summary(mod)
anova(mod)


plot(mod)
summary(mod)
screeplot(mod)
RsquareAdj(mod)


# Make plot ---------------------------------------------------------------

summm <- summary(mod)
plot(mod)
xtext <- paste("(", round(summm$concont[[1]][2, 1] * 100, 2), " % of fitted, ", 
               round(summm$cont[[1]][2, 1] * 100, 2), " % of total variation)", sep = "")
ytext <- paste("(",
               round(summm$cont[[1]][2, 2] * 100, 2), " % of total variation)", sep = "")
mtext(xtext, 1, line = 4)
mtext(ytext, 2, line = 2)
