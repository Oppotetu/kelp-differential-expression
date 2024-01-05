# Safekelp RNA-seq analysis 2

rm(list=ls())

library(tidyverse)
library(waldo)
library(edgeR)
library(DESeq2)
library(WGCNA)
library(parallel)
library(egg)

setwd("C:/Users/simon/Prosjekter/RNA_seq/Workshop_analysis")

dta<- read.table("counts.txt", header = T)
names(dta) 

#creating metadata for the counts
mdta<-data.frame(id=c('S1','S49','S2','S38','S7','S43','S8','S44'),
                 sample=c('1EK','2EK','1ET','2ET','1GK','2GK','1GT','2GT'),
                 rep=c(1,2,1,2,1,2,1,2),
                 group1=c('EK','EK','ET','ET','GK','GK','GT','GT'),
                 cultivar=c('Engmo','Engmo','Engmo','Engmo','Grindstad','Grindstad','Grindstad','Grindstad'),
                 treatment=c('K','K','T','T','K','K','T','T'),
                 stringsAsFactors = TRUE)

rownames(mdta)<-mdta$id

##checking the order of samples in metadata and count data
waldo::compare(colnames(dta), rownames(mdta)) # no differences will output "✔ No differences"

#preparing DEGlist
FtDEG<-DGEList(counts = dta, samples = mdta, remove.zeros = T )
dim(FtDEG) # no. of rows(transcripts) & columns(samples)

##Filtering transcripts with low read counts
keep<-filterByExpr(FtDEG,min.count=10, group = FtDEG$samples$group1 )
table(keep)

FtDEG <- FtDEG[keep, ,keep.lib.sizes=FALSE] #subset DEGlist based the transcripts that passed filtering

## load the Differentially expressed genes list

fls = list.files(pattern = "*_ctrl.txt", full.names = F) #lists all files in the current directory matching the pattern

dfs<-lapply(fls, FUN=read.table, sep=" ", header=T) # load all files in the fls vector

names(dfs)<- substr(fls,1, nchar(fls)-4) # naming the dataframes in list

DEgenes <- unique(unlist(lapply(dfs, "[",  , c("GeneID")))) # extracting the DE genes from all contrasts of interest

#normalization
FtDEG<- calcNormFactors(FtDEG, method = "TMM")
FtDEG$samples$norm.factors

## Variance Stablizing transformation
vsd<- vst(round(FtDEG$counts), blind = T) #default blind=T(for totally unsupervised clustering)

## cleaning samples and genes with too many missing values
gsg = goodSamplesGenes(vsd, verbose = 3)

gsg$allOK #TRUE, all genes have passed the cuts

## sample clustering to detect outliers
plot(hclust(dist(t(vsd))), labels=FtDEG$samples$sample, main = " ", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 1)

#calculating the variance of the logcounts (genes are rows)
var_genes <- apply(vsd, 1, var)

#sorting the genes based on variance explained and picking top 500 genes for plotting
select_var <- names(sort(var_genes, decreasing=TRUE))[1:10000]

## subset DE genes from counts matrix 
norm_DEcounts<-vsd[rownames(vsd) %in% DEgenes,]

## reshaping matrix to dataframe 
##reshaping data to long format to plot
norm_DEcounts_df <- data.frame(norm_DEcounts) %>%
  mutate(Gene_id = row.names(norm_DEcounts)) %>%
  pivot_longer(-Gene_id)


#plotting 
norm_DEcounts_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 45)
  ) + scale_x_discrete(labels= FtDEG$samples$sample)+
  ylim(0, NA) +
  labs(
    title = "Normalized and 95 quantile Expression",
    x = "treatment",
    y = "normalized expression"
  )

#WGCNA
#transpose data to feed into WGCNA
input_mat = t(norm_DEcounts)

# enabling multi-threading
cores<-detectCores()/2
enableWGCNAThreads(nThreads = cores)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

#####################
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence"))+
  text(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = 1, col = "red")+
  abline(h = 0.80, col = "red")

plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"))+
  text(sft$fitIndices[, 1],
       sft$fitIndices[, 5],
       labels = powers,
       cex = 1, col = "red")

picked_power = 7
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 1000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = F,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

table(netwk$colors)

cor <- temp_cor

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

#Relate Module (cluster) Assignments to Treatment Groups
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
) 

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes # all the genes that best represent this module

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$group = mdta$sample
MEs0$treatment<-mdta$treatment

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(name = factor(gsub("ME", "", name), levels = module_order)
  )

# module trait relationship based on treatment
'mME %>% mutate(treatment = factor(group, levels=c("1EK","2EK","1GK","2GK",
                                                      "1ET","2ET","1GT","2GT")))%>%'
mME %>% mutate(treatment = factor(treatment, levels=c("K", "T")))%>%
  ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships(Treatment)", y = "Modules", fill="corr")

## module trait relationship based on cultivar
'mME %>% mutate(treatment = factor(group, levels=c("1EK","2EK","1ET","2ET","1GK","2GK","1GT","2GT")))%>%
  ggplot(., aes(x=treatment, y=name, fill=value)) +'
mME %>% mutate(treatment = factor(treatment, levels=c("E", "G")))%>%
  ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships(cultivar)", y = "Modules", fill="corr")

# pick out modules of interest
modules_of_interest = c("turquoise","blue", "yellow")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# genes in each module
table(submod$colors)

#extracting genes from mdoules of interest
subexpr = norm_DEcounts[submod$gene_id,]

#reshaping dataframe for plotting
submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

#attaching metadata
submod_df$group<- mdta$group1[match(submod_df$name, mdta$id)]
submod_df$id<-  mdta$sample[match(submod_df$name, mdta$id)]

#plotting based on treatment
im1<-submod_df%>% 
  ggplot(., aes(x=factor(id, level = c('1EK','2EK','1GK','2GK','1ET','2ET','1GT','2GT')),y=value,group=gene_id)) +
  geom_line(aes( color=module),alpha = 0.2) +
  scale_colour_manual(values = c(turquoise="turquoise",blue="blue", yellow="yellow"))+
  theme_bw() + theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, size = 7.5, face="bold"),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  facet_grid(rows = vars(module)) +
  labs(x = NULL, y = "Normalized expression")

my_tag <- c("n=461","n=2572") # alphabetical order of colors

im11<-tag_facet(im1,
                x=7.8, y=14.9,
                fontface = 4,
                size = 3,
                family = "serif",
                open = "", close = "",
                tag_pool = my_tag,colour = "black" )
im11








