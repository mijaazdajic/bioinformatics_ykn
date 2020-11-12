###############################################################################################
###########THIS R SCRIPT IS FOR IMPORTING THE QIIME2 DATA FROM THE PIPELINE####################
#####The importing is done with phyloseq, then other pkgs can be use for analysis##############
###############################################################################################
#https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121

# Loaded needed packages --------------------------------------------------
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(vegan)
library(dplyr)
library(janitor) # lots of help when dealing with non-R-friendly datasets
library(lubridate) # for working with dates
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(grid)
library(cowplot)




# The objects that import are results of the QIIME2 pipeline as -----------
# well as the metadata ----------------------------------------------------
# This is one way of doing it - if you want to use qiime2R to make --------
# graphs with ggplot2 -----------------------------------------------------

metadata <- read_tsv("QIIME2/metadata.tsv") #load metadata table
#for rarefied data
#SVs <-  read_qza("QIIME2/14_dada2_table_final.qza") #load the feature table
#for non rarefied data
SVs <-  read_qza("QIIME2/Silva_13.8/13_dada2_table_filt_contam.qza") #load the feature table

taxonomy <- read_qza("QIIME2/Silva_13.8/classification.qza") #load taxonomy
taxtable<-taxonomy$data %>% as_tibble() %>% separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) #convert the table into a tabular split version
taxtable
tree <-  read_qza("QIIME2/Silva_13.8/asvs-tree.qza") #adding the rooted tree made in qiime2
shannon <-  read_qza("QIIME2/shannon_vector.qza") #adding shannon index
pco <- read_qza("QIIME2/unweighted_unifrac_pcoa_results.qza") #adding the PCoA made in qiime2



# Use phyloseq now to build a phyloseq object -----------------------------
physeq<-phyloseq(
  otu_table(SVs$data, taxa_are_rows = T), 
  phy_tree(tree$data), 
  tax_table(as.data.frame(taxtable) %>% 
              dplyr::select(-c(Confidence)) %>% column_to_rownames("Feature.ID") %>% as.matrix()), #moving the taxonomy to the way phyloseq wants it
  sample_data(metadata %>% as.data.frame() %>% column_to_rownames("#SampleID"))
)

physeq

# Now your phyloseq data is ready to be worked with (Hopefully!) ----------
# Check if data is normalized ---------------------------------------------
# it is not....
p <- plot_bar(physeq,x="Description" ,fill="Class")
p <- p + geom_bar(aes(color=Class), stat="identity", position="stack")
p <- p + theme(legend.title=element_blank(), legend.text = element_blank())
p

#Remove taxa not seen more than 3 times in at least 50% of the samples. This protects 
#against an OTU with small mean & trivially large C.V.
GP1 = filter_taxa(physeq, function(x) sum(x > 3) > (0.5*length(x)), TRUE)


#physeq dataset in which the abundance values have been transformed to relative abundance 
#within each sample
GP2 = transform_sample_counts(GP1, function(x) 1E2 * x/sum(x))


#check again..
p <- plot_bar(GP2,x="Description" ,fill="Class")
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
  geom_text_repel(aes(label = Description))



#NMDS
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
  geom_text_repel(aes(label = Description)) #NMDS does seem to show a difference in Sulphate along the first axis

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
                      label = "Description" # not really needed but can be useful for exploration
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
                      label = "Description" # not really needed but can be useful for exploration
) +
  geom_point(aes(color = Sulphate), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5)

# combine plots p1 and p2 into one figure with the plot_grid command (cowplot)
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

# constrained ordination test using Correspondence Analysis (CA)

ykn_CA <- ordinate(GP1, "CCA")

# check ordination with a scree plot
plot_scree(ykn_CA, "Scree plot of ykn scaled Correspondence analysis")



(p1_CA <- plot_ordination(GP2, ykn_CA, "samples",
                          color="Description",label="Description") +
    geom_point(aes(color = Description), alpha = 0.4, size = 4))


# Now doing a constrained Correspondence Analysis (CCA), using time
ykn_CCA <- ordinate(GP2, formula = GP2 ~ Sulphate + DOC, "CCA")

# check ordination with a scree plot
plot_scree(ykn_CCA, "Scree plot of ykn scaled Constrained Correspondence analysis")

# CCA plot
CCA_plot <- plot_ordination(GP2, ykn_CCA, type="samples", color="Description", 
                            label="Description") +
  geom_point(aes(color = Description), alpha = 0.4, size = 4)

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


















# prune taxa to only top 5 ------------------------------------------------
top5 <- prune_taxa(names(sort(taxa_sums(physeq),TRUE)[1:5]), physeq)
plot_tree(top5)
view(map)

# prune taxa to only top 30 -----------------------------------------------
top30 <- prune_taxa(names(sort(taxa_sums(physeq),TRUE)[1:30]), physeq)
plot_tree(top30)
view(map)

# prune taxa to only top 30 -----------------------------------------------
top50 <- prune_taxa(names(sort(taxa_sums(physeq),TRUE)[1:50]), physeq)
plot_tree(top50)
plot_tree(top50, color="Description", shape="Class", size="abundance")
plot_heatmap(top50, sample.label="Description", taxa.label = "Class")


#Can also plot map with the family names on side, and change the distance matrix
p <- plot_heatmap(top30, "NMDS", "bray", "Description", "Class")
print(p)



# This section will be from a tutorial online, filtering and analy --------
#https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#using_phyloseq


# Preporcessing data ----------------------------------------------------------
#RThe most abundant taxa are kept only if they are in the most abundant 20% of taxa 
#in at least half of the samples in dataset

f1<- filterfun_sample(topp(0.1))
wh1 <- genefilter_sample(physeq, f1, A=(1/2*nsamples(physeq)))
sum(wh1)

ex2 <- prune_taxa(wh1, physeq)

#transform sample count to relative abundance
rel  = transform_sample_counts(ex2, function(x) x / sum(x) )
ex2 = filter_taxa(ex2, function(x) sum(x) > .005, TRUE)
phy = phyloseq(otu_table(ex2))




# Show available ranks in the dataset
rank_names(physeq)

# Making richness figures -------------------------------------------------
#Side by side figures of different richness values
plot_richness(physeq, measures = c("Chao1", "Shannon"))
# Create table, number of features for each phyla
table(tax_table(physeq)[, "Phylum"], exclude = NULL)

physeq <- subset_taxa(physeq, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) #The following ensures that features with ambiguous phylum annotation are also removed.


# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(physeq),
               MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(physeq),
                    tax_table(physeq))
#Compute the total and average prevalences of the features in each phylum
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#Filter out samples that have a low count from above
# Define phyla to filter
filterPhyla = c("D_1__Acetothermia", "D_1__Chlamydiae", "D_1__Margulisbacteria", "D_1__WS2")
# Filter entries with unidentified Phylum.
phy1 = subset_taxa(physeq, !Phylum %in% filterPhyla)
phy1



# Prevalence Filtering ----------------------------------------------------
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(phy1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(physeq),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#there isn't a "natural"seperation in the abundance that we could use, will use 5% in samples as threshold
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(physeq)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
phy2 = prune_taxa(keepTaxa, physeq)



# Agglomerate taxa --------------------------------------------------------
#Agglomerate the data features corresponding to closely related taxa
# How many genera would be present after filtering?
length(get_taxa_unique(phy2, taxonomic.rank = "Genus")) #362
phy3 = tax_glom(phy2, "Genus", NArm = TRUE)

#Compare the original unfiltered data and the tree after taxonoic agglomeration
multiPlotTitleTextSize = 15
p2tree = plot_tree(phy2, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(phy3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))

# group plots together
grid.arrange(nrow = 1, p2tree, p3tree)



# Abundance value transformation ------------------------------------------
#make the function
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("D_1__Cyanobacteria"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Description",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

# Transform to relative abundance. Save as new object.
ps3ra = transform_sample_counts(phy3, function(x){x / sum(x)})

#plot the abundance values before and after transformation
plotBefore = plot_abundance(phy3,"")
plotAfter = plot_abundance(ps3ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2,  plotBefore, plotAfter)


# Subset by taxonomy ------------------------------------------------------
psOrd = subset_taxa(ps3ra, Order == "D_3__Syntrophobacterales")
plot_abundance(psOrd, Facet = "Order", Color = NULL)



.cran_packages <- c( "shiny","miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest", "dplyr", "ggrepel", "nlme", "devtools",
                     "reshape2", "PMA", "structSSI", "ade4",
                     "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("genefilter", "impute")
# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)){
  install.packages(.cran_packages[!.inst],repos = "http://cran.rstudio.com/")
}
.inst <- .github_packages %in% installed.packages()
if (any(!.inst)){
  devtools::install_github(.github_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)){
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst])
}



# Preprocessing -----------------------------------------------------------
#redo this part when we have continous data in mappingfile
quickplot(sample_data(physeq)$Description, geom = "histogram",binwidth=20) + xlab("Description")



qplot((rowSums(otu_table(phy3))),binwidth=0.2) +
  xlab("Logged counts-per-sample")


#sample_data(physeq)$Description <- cut(sample_data(physeq)$Description,
#                                  breaks = c(0, 100, 200, 400))
#levels(sample_data(physeq)$Description) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
#sample_data(physeq)$Description=gsub(" ","",sample_data(physeq)$Description)
pslog <- transform_sample_counts(top50, function(x) log(1 + x))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "Description") +
  labs(col = "Lake") +
  coord_fixed(sqrt(evals[2] / evals[1]))



# Different Ordination Projections ----------------------------------------

#do some unsupervised ordinations for exploratory purposes
#we can prune some samples if we want (I wont do it now...)
outliers <- c("F5D165", "F6D165", "M3D175", "M4D175", "M5D175", "M6D175")
ps <- prune_samples(!(sample_names(ps) %in% outliers), ps)


#We are also going to remove samples with fewer than 1000 reads
#which(!rowSums(otu_table(physeq)) > 6000)


#ps <- prune_samples(rowSums(otu_table(physeq)) > 6000, physeq)
pslog <- transform_sample_counts(physeq, function(x) log(1+x))
#same as above only bray distance
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")
evals <- out.pcoa.log$values[,1]
plot_ordination(pslog, out.pcoa.log, color = "Description") +
  labs(col = "Lake")+
  coord_fixed(sqrt(evals[2] / evals[1]))



#Make a DPCoA
out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(pslog, out.dpcoa.log, color = "Description", label= "Description") +
  labs(col = "Description")+
  coord_fixed(sqrt(evals[2] / evals[1]))


-------------------------------------------------------------------------------------------------------------------------

# this sections is just extra thoughts ------------------------------------


-------------------------------------------------------------------------------------------------------------------------
  

# PCA on ranks ------------------------------------------------------------
#Microbial abundance data is often heavy-tailed, and sometimes it can be 
#hard to identify a transformation that brings the data to normality. In 
#these cases, it can be safer to ignore the raw abundances altogether, 
#and work instead with ranks.
#I DONT THINK I WILL USE THIS SECTION
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))


abund_ranks <- abund_ranks - 1000
abund_ranks[abund_ranks < 1] <- 1


abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = 1.5) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")






# Canonical correspondence ------------------------------------------------
ps_ccpna <- ordinate(pslog, "CCA", formula = pslog ~ pH + DOC + Sulphate)

ps_scores <- vegan::scores(ps_ccpna)

sites <- data.frame(ps_scores$sites)

sites$SampleID <- rownames(sites)

sites <- sites %>% merge(sample_data(pslog),by= "row.names") %>% column_to_rownames(var = "Row.names")

species <- data.frame(ps_scores$species)

species$otu_id <- seq_along(rownames(otu_table(physeq)))

species <- species %>% merge(tax_table(pslog),by="row.names") %>% column_to_rownames(var = "Row.names")

evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)



ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, col = Phylum), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                  aes(x = CCA1, y = CCA2, label = otu_id),
                  size = 1.5, segment.size = 0.1) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_fill_discrete() +
  coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45   ) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))


https://peerj.com/articles/2836/
  https://www.biorxiv.org/content/10.1101/217919v3.full


########Build an NNMD with OTU data and environemntal variables and Km/Kd
https://rpubs.com/faysmith/fung_soil_survey

# constrained ordination test using Correspondence Analysis (CA)

mouse_CA <- ordinate(top50, "CCA")

# check ordination with a scree plot
plot_scree(mouse_CA, "Scree plot of YKN microbial analysis")

(p1_CA <- plot_ordination(top50, mouse_CA, "samples",
                          color="Sulphate", label="Description") +
    geom_point(aes(color = Sulphate), alpha = 0.4, size = 4))



# Now doing a constrained Correspondence Analysis (CCA), using time
mouse_CCA <- ordinate(top50, formula = top50 ~ Sulphate + pH + DOC + Iron_Total, "CCA")

# check ordination with a scree plot
plot_scree(mouse_CCA, "Scree plot of mouse scaled Constrained Correspondence analysis")

# CCA plot
CCA_plot <- plot_ordination(top50, mouse_CA, "samples",
                          color="Sulphate", label="Description") +
    geom_point(aes(color = Sulphate), alpha = 0.4, size = 4)

# Now add the environmental variables as arrows into a matrix
arrowmat <- vegan::scores(mouse_CCA, display = "bp")

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
anova(mouse_CCA)

# run envit for testing
map <- metadata %>% select(Sulphate,pH,Percent_Hg,DOC,Iron_Total)
ef <- envfit(mouse_CA, map, permutations = 999, na.rm = TRUE)

# plot the ordination data directly using plot
plot(mouse_CCA, display="sites" )

# overly plot with fitted variables
plot(ef, p.max= 0.001)


# export a csv file of all the samples metadata ------------------------------------
glom <- tax_glom(GP2, taxrank='Genus')
dat <- psmelt(glom)
write.csv(dat, file='otus-genus-with-sample-data.csv')
