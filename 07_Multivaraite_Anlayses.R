###############################################################################################
###########THIS R SCRIPT IS FOR ANALYSING THE QIIME2 DATA FROM THE PIPELINE####################
#####The importing is done with phyloseq, then other pkgs can be use for analysis##############
###############################################################################################


# Loaded needed packages --------------------------------------------------
library(dplyr)
library(stringr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggord)

#

# resolve conflicts -------------------------------------------------------
conflict_prefer("rename", "dplyr")

#load data from previous script

metadata <- read_tsv("QIIME2/metadata.tsv") %>% rename(lake_name="Description")#load metadata table

suva <- read_csv("CHEM_DATA/SUVA_Analysis.csv", local = locale(encoding = "latin1"))

years <- read_csv("INFO_SHEETS/Sampling_Year_Per_Lake.csv", locale = locale(encoding = "latin1"))

methylation_table <- read_csv("methylation_table.csv", local = locale(encoding = "latin1"))

methylation_table <- methylation_table %>% merge(.,years, by = "lake_name", all=T)

env.data <- left_join(metadata, methylation_table, by = c("lake_name","Year")) %>% 
  left_join(., suva, by = "lake_name", all=T) %>% 
  select(.,c('#SampleID', lake_name, MeHg_water_ngl, Year, THg_water_ngl, Percent_Hg, DOC, pH,
              Sulphate, Iron_Total, Dissolved_Phosphorus_UO, Total_Phophorus_UO, 
              Dissolved_Arsenic_UO, Total_Arsenic_UO, km_24, kd_24, SUVA, SUVA_percent)) %>% 
  na.omit()


spec.data <- as.data.frame(t(otu_table(GP2))) %>% dplyr::filter(row.names(.) %in% env.data$`#SampleID`) #subset the data so its the same rows as metadata



#Prune the phyloseq object to only have samples with Km and Kd

lakes <-  unique(env.data$lake_name) %>% as.vector()

env.data2 <- left_join(metadata, methylation_table, by = c("lake_name","Year")) %>% 
  left_join(., suva, by = "lake_name", all=T) %>% 
  select(.,c('#SampleID', lake_name, MeHg_water_ngl, Year, THg_water_ngl, Percent_Hg, DOC, pH,
             Sulphate, Iron_Total, Dissolved_Phosphorus_UO, Total_Phophorus_UO, 
             Dissolved_Arsenic_UO, Total_Arsenic_UO, km_24, kd_24, SUVA, SUVA_percent))

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
  sample_data(env.data2 %>% as.data.frame() %>% column_to_rownames("#SampleID"))
)

physeq2

sub <- subset_samples(physeq2, lake_name %in% c("BC36", "BC18", "ICING", "PONTOON", "YKS2",
                                               "YKN1", "YKW1", "PROSPEROUS", "YKS1", "BC43",
                                               "YK40", "YK42", "BC21", "YK67"))


sub

env.data <- as.data.frame(sample_data(sub))


# filtering samples -------------------------------------------------------

# Now your phyloseq data is ready to be worked with (Hopefully!) ----------
# Check if data is normalized ---------------------------------------------
# it is not....
p <- plot_bar(sub,x="lake_name" ,fill="Class")
p <- p + geom_bar(aes(color=Class), stat="identity", position="stack")
p <- p + theme(legend.title=element_blank(), legend.text = element_blank())
p


#Remove taxa not seen more than 3 times in at least 20% of the samples. This protects 
#against an OTU with small mean & trivially large C.V.
GP1 = filter_taxa(sub, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

#physeq dataset in which the abundance values have been transformed to relative abundance 
#within each sample
GP2 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))


#check again..
p <- plot_bar(GP2,x="lake_name" ,fill="Class")
p <- p + geom_bar(aes(color=Class), stat="identity", position="stack")
p <- p + theme(legend.title=element_blank(), legend.text = element_blank())
p


###Unconstrained ordinations
# Ordinate using Principal Coordinate analysis
ykn_pcoa <- ordinate(
  physeq = GP2, 
  method = "PCoA", 
  distance = "bray"
)

# checking the ordination by making a scree plot
plot_ordination(
  physeq = GP2,
  ordination = ykn_pcoa,
  type="scree")

# Plot PCoA
plot_ordination(
  physeq = GP2,
  ordination = ykn_pcoa,
  axes=c(1,2),   # this selects which axes to plot from the ordination
  color = "Sulphate",
  title = "PCoA of all YKN samples"
) +
  geom_point(aes(color = Sulphate), alpha = 0.7, size = 4) +
  geom_point(colour = "Sulphate", size = 1.5)+
  geom_text_repel(aes(label = lake_name))


# NMDS --------------------------------------------------------------------
# running ndms on GP2
GP2_otu_nmds <- ordinate(GP2, "NMDS", "bray")

# printing stress level
cat("stress is:", GP2_otu_nmds$stress)


#ordination plot using the otus, and color by Phylum
plot_ordination(GP2, GP2_otu_nmds, type="taxa", color="Phylum", title="NMDS ordination - OTUs")


#ordination plot using the lakes, and color by Phylum
plot_ordination(GP2, GP2_otu_nmds,
                type="samples",
                color="Sulphate",
                title="NDMS ordination - Samples") +
  geom_point(aes(color = Sulphate), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) +
  geom_text_repel(aes(label = lake_name)) #NMDS does seem to show a difference in Sulphate along the first axis


#UNIFRAC analysis of community composition
# creating a distance matrix using unweigthed unifrac distances
ykn_Unifrac_unweight <- UniFrac(GP2, weighted = FALSE)

# calculation of the NMDS using the unweigthed unifract distances
ykn_unifrac_unweight_nmds <- ordinate(GP2, "NMDS", ykn_Unifrac_unweight)

cat("stress is:", ykn_unifrac_unweight_nmds$stress) #ideally it would be below 0.01

# Plotting the ordination, 
p1 <- plot_ordination(GP2, ykn_unifrac_unweight_nmds,
                      type="samples",
                      title="NDMS unifrac Unweighted",
                      label = "lake_name" # not really needed but can be useful for exploration
) +
  geom_point(aes(color = Sulphate), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5)

# plotting object p1
print(p1)


## create weighted unifrac ordination and compare with unweighted
ykn_unifrac_W <- UniFrac(GP2, weighted = TRUE, 
                         parallel = TRUE)
# calculation of the NMDS using the unweigthed unifract distances
ykn_unifrac_W_nmds <- ordinate(GP2, "NMDS", ykn_unifrac_W)

cat("stress is:", ykn_unifrac_W_nmds$stress)

# Plotting the ordination
p2 <- plot_ordination(GP2, ykn_unifrac_W_nmds,
                      type="samples",
                      title="NDMS unifrac Weighted",
                      label = "lake_name" # not really needed but can be useful for exploration
) +
  geom_point(aes(color = Sulphate), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5)

# combine plots p1 and p2 into one figure with the plot_grid command (cowplot)
print(p2)

figure = plot_grid(p1, p2, labels = c("P1", "P2"), ncol = , nrow = 1, rel_widths=c(2,2))

ggsave(figure, file="UNIFRAC_COMPARE_SULPHATE.svg", dpi = 600, path = "Figures/")


# testing of significance for the unifrac ordinations using ANOSIM --------
#this is only good for group data type (only used sulphate to try for now)
# create vector with time point labels for each samples
ykn_sulphate <- get_variable(GP2, "Sulphate")

# run Anosim on a unweighted Unifrac distance matrix of the mouse_scaled OTU table, with the factor time
ykn_anosim_UFuW <- anosim(distance(GP2, "unifrac", weighted=FALSE), ykn_sulphate)

# run Anosim on a Weighted Unifrac distance matrix of the mouse_scaled OTU table, with the factor time
ykn_anosim_UFW <- anosim(distance(GP2, "unifrac", weighted=TRUE), ykn_sulphate)

#plot results
print(ykn_anosim_UFuW)
print(ykn_anosim_UFW)


# Run a PERMANOVA for data (non-perametric analogue to ANOVA) -------------
# run a permanova test with adonis.
set.seed(1)

# Calculate unifrac
ykn_unifrac_W <- phyloseq::distance(GP2, method = "unifrac", weighted=TRUE)

# make a data frame from the scaled sample_data
sampledf <- data.frame(sample_data(GP2))

# Adonis test
adonis(ykn_unifrac_W ~ Sulphate, data = sampledf) #there is a sig with DOC


# test of Homegeneity of dispersion
beta <- betadisper(ykn_unifrac_W, sampledf$Sulphate)

# run a permutation test to get a statistic and a significance score
permutest(beta) #didnt work with DOC, but 0.3 with sulphate



# Constrained Ordinations -------------------------------------------------

# Make a distance matrix using UNIFRAC weighted ---------------------------

wUF.dm<-UniFrac(sub, weighted = TRUE, fast = TRUE) 
wUF.table<-as.matrix(dist(wUF.dm))



# Calculate distance based RDA using Weighted UniFrac distance mat --------
scaled.env <-  env.data %>% 
  select(-"lake_name") %>% #take out row names that are not num so can be standardized
  decostand(., method = "standardize")



dbRDA<-capscale(wUF.dm~ DOC + pH + Sulphate + Iron_Total + km_24 + kd_24, data=scaled.env)

dbRDA2 <- dbrda(wUF.dm ~ DOC + pH + Sulphate + Iron_Total, data=scaled.env)
plot(dbRDA, main = "dbRDA")
summary(dbRDA)
screeplot(dbRDA)
RsquareAdj(dbRDA)##Variance partition


vegan::anova.cca(dbRDA, by="term", step = 200)
vif.cca(dbRDA) # a vif of more than 10 indicates that a variable is strongly dependent on others

# Do a stepwise regression to find variables that best explain species differeces
rda.1 <- capscale(wUF.dm~1, data = scaled.env)
mod <- ordistep(rda.1, scope = formula(dbRDA),
                direction = "both", permutations = 999, trace = T) #bidirectional selection

anova(mod)


step.mod <-  ordiR2step(rda.1, scope = formula(mod), direction = "both", 
                        permutations = 999, trace = T)



plot(mod)
summary(mod)
screeplot(mod)
RsquareAdj(mod)
vegan::anova.cca(mod, step = 200)


# plotting the db-RDA ordination with DOC as vector -----------------------

figure.full.km <- autoplot(dbRDA, arrows = T) + geom_point(size=50)
figure.full.km

figure.reduced.km <- autoplot(mod, arrows = T)
figure.reduced.km


















































# constrained ordination test using Correspondence Analysis (CA)

ykn_CA <- ordinate(GP1, "CCA")

# check ordination with a scree plot
plot_scree(ykn_CA, "Scree plot of ykn scaled Correspondence analysis")



(p1_CA <- plot_ordination(GP2, ykn_CA, "samples",
                          color="lake_name",label="lake_name") +
    geom_point(aes(color = lake_name), alpha = 0.4, size = 4))


# Now doing a constrained Correspondence Analysis (CCA), using time
ykn_CCA <- ordinate(GP2, formula = GP2 ~ kd_24 + Iron_Total, "CCA")

# check ordination with a scree plot
plot_scree(ykn_CCA, "Scree plot of ykn scaled Constrained Correspondence analysis")

# CCA plot
CCA_plot <- plot_ordination(GP2, ykn_CCA, type="samples", color="lake_name", 
                            label="lake_name") +
  geom_point(aes(color = lake_name), alpha = 0.4, size = 4)

# Now add the environmental variables as arrows into a matrix
arrowmat <- vegan::scores(ykn_CCA, display = "bp")

# transform matrix into a dataframe, and add labels
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)


# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CCA1, 
                 yend = CCA2, 
                 x = 0, 
                 y = 0,
                 shape= NULL,
                 color= NULL,
                 label=labels)

label_map <- aes(x = 1.2 * CCA1, 
                 y = 1.2 * CCA2,
                 shape= NULL,
                 color= NULL,
                 label=labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
CCA_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )

# permutational anova test on constrained axes
anova(ykn_CCA, by="margin")


##http://genoweb.toulouse.inra.fr/~formation/15_FROGS/5-June2016/FROGS_phyloseq_23062016.pdf
met <- as(sample_data(GP2), "data.frame")
cap <-  capscale(ykn_unifrac_W ~ Sulphate, data = met)

anova <- anova(cap, permutations = 999)

adonis(ykn_unifrac_W~kd_24, data = met, permutations = 999)




p<- autoplot(otu.rda.f, arrows = TRUE,axes = c(1, 2), geom = "text", layers = c( "species","sites", "biplot", "centroids"), legend.position = "right", title = "db-RDA")







#Calculate Weighted UniFrac distance matrix 
wUF.dm<-UniFrac(GP2, weighted = TRUE, fast = TRUE) 
wUF.table<-as.matrix(dist(wUF.dm))


#Calculate distance based RDA using Weighted UniFrac distance matrix 
dbRDA<-capscale(wUF.dm ~ DOC + pH, data=met, comm = ) 
plot(dbRDA, main = "dbRDA")
summary(dbRDA)
screeplot(dbRDA)




































# hellinger transform the species dataset: gives low weights to rare species 
spe.hel <- decostand(spec.data, "hellinger")

# Calculate distance matrix
bc<-vegdist(spe.hel, method="bray", binary=FALSE) 

# look at an unconstrained ordination first, it is always a good idea to look at both unconstrained and constrained ordinations
# set the seed: to reproduce the same result in the fture
set.seed(100)
bci.mds<-metaMDS(spe.hel, distance = "bray", k = 2)

# extract x and y coordinates from MDS plot into new dataframe, so you can plot with ggplot 
MDS_xy <- data.frame(bci.mds$points)
bci.mds$stress # 0.1241412

# colour by island
ggplot(MDS_xy, aes(MDS1, MDS2, col=env.data$Sulphate)) + geom_point() + theme_bw() + ggtitle('stress:0.175')



# RDA

simpleRDA <- rda(spec.data ~ Sulphate + pH + DOC + Iron_Total + SUVA_percent + km_24 + kd_24 + 
                   Total_Phophorus_UO + Total_Arsenic_UO, data=env.data)
summary(simpleRDA)
screeplot(simpleRDA) #bstick not available for constrained ordinations

# canonical coefficients
coef(simpleRDA)

# unadjusted R^2 retreived from the rda result
R2 <- RsquareAdj(simpleRDA)$r.squared
R2 

# adjusted R^2
R2adj <- RsquareAdj(simpleRDA)$adj.r.squared
R2adj 

# Triplot: three different entities in the plot: sites, response variables and explanatory variables (arrowheads are on the explanatory variables)
# Scaling 1
plot(simpleRDA, scaling=1, main="Triplot RDA matrix ~ env - scaling 1 - wa scores")

# arrows for species are missing, so lets add them without heads so they look different than the explanatory variables
spe.sc <- scores(simpleRDA, choices=1:2, scaling=1, display="sp")
arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='red')


# Scaling 2
plot(simpleRDA, main="Triplot RDA matrix ~ env - scaling 2 - wa scores")
spe2.sc <- scores(simpleRDA, choices=1:2, display="sp") # scores() choices= indicates which axes are to be selected, make sure to specify the scaling if its different than 2 
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='red')


# plot the RDA using ggplot (ggord package)
ggord(simpleRDA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # looking at the raw code, this is plotting the 'wa scores', the blue dots are different species  













