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


# qualitative analysis for personal info ----------------------------------
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


########================================================== DONE ================================#################
