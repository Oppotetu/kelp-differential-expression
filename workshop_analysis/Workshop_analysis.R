# Safekelp RNA-seq analysis

library(tidyverse)
library(data.table)
library(waldo)
library(readxl)
library(edgeR)
library(DESeq2)
library(adegenet)
library(RColorBrewer)
library(adegenet)
library(ade4)
library(factoextra)
library(UpSetR)
library(ggvenn)
library(VennDetail)
library(affy)

setwd("C:/Users/simon/Prosjekter/RNA_seq/Workshop_analysis")

# INKSPACE

#loding counts matrix file
dta <- read.table("counts.txt", header = T)
dim(dta)
names(dta)

#creating metadata for the counts
mdta <- data.frame(
  id=c('S1','S49','S2','S38','S7','S43','S8','S44'),
  sample=c('1EK','2EK','1ET','2ET',
           '1GK','2GK','1GT','2GT'),
  rep=c(1,2,1,2,1,2,1,2),
  group1=c('EK','EK','ET','ET','GK','GK','GT','GT'),
  cultivar=c('Engmo','Engmo','Engmo','Engmo',
             'Grindstad','Grindstad','Grindstad','Grindstad'),
  treatment=c('K','K','T','T','K','K','T','T'),
  stringsAsFactors = TRUE)

rownames(mdta) <- mdta$id

##checking the order of samples in metadata and count data
waldo::compare(
  colnames(dta), 
  rownames(mdta))

#preparing DEGlist
FtDEG <- DGEList(
  counts = dta, 
  samples = mdta, 
  remove.zeros = T)
dim(FtDEG)

##plotting distribution of counts
plotDensities(
  cpm(
    FtDEG$counts, 
    log = T), 
  col = c("red", "blue", "green", "yellow", 
          "pink", "orange","black", "gray"), 
  legend = "topright")

##Filtering transcripts with low read counts
keep <- filterByExpr(
  FtDEG,
  min.count=10, 
  group = FtDEG$samples$group1)

table(keep)  # if 5 replicates : 30, if what do you want to look at? big patterns: throw away lowly expressed genes. isoforms: keep them. 

FtDEG <- FtDEG[keep, ,keep.lib.sizes=FALSE]
dim(FtDEG)

##plotting distribution of counts 2
plotDensities(
  log(FtDEG$counts),
  col = c("red", "blue", "green", "yellow", 
          "pink", "orange", "black", "gray"), 
  legend = "topright" )

#mean lib.size
lib.size <- as.character(
  round(mean(FtDEG$samples$lib.size * 1e-6), 
        2))

#plotting library sizes
barplot(
    FtDEG$samples$lib.size*1e-6,
    names=FtDEG$samples$sample,
    ylab="Library size (millions)") + 
  abline(
    h=mean(FtDEG$samples$lib.size*1e-6), 
    col="Red", 
    lty=5, lwd=1) + 
  text(
    x=8.5, 
    y=20.5, 
    paste(
      lib.size, 
      "million reads", 
      sep = " "), 
    col= "black")

#plotting transcripts counts
barplot(
  rowSums(FtDEG$counts*1e-6), 
  las=2, 
  main="Counts per transcript", 
  axisnames = FALSE, 
  ylab = "counts in millions", 
  cex.axis=0.8)

#normalization
FtDEG <- calcNormFactors(
  FtDEG, 
  method = "TMM")

FtDEG$samples

#logcpm values
lcpm <- cpm(
  FtDEG, 
  prior.count = 2, 
  log=TRUE)

boxplot(
  lcpm,
  xlab= "samples", 
  ylab= "log-cpm")

## Variance Stablizing transformation 
vsd <- vst(
  round(FtDEG$counts), 
  blind = T) #default blind=T(for totally unsupervised clustering)

## Hierarchical Clustering 
sampleTree <- hclust(
  dist(
    t(vsd)))

plot(
  sampleTree,
  labels = FtDEG$samples$group2, 
  main = "Hierarchical Clustering", 
  sub = "", 
  xlab = "", 
  cex.lab = 1.5,
  cex.axis = 1.5, 
  cex.main = 2)

# compute pairwise correlations
vsd_cor <- cor(vsd)

#assigning colors 
ann_colors = list(
  cultivar = c(Engmo = "orange2",  
               Grindstad= "darkgreen"),
  treatment = c(K = "#084594", T = "#6BAED6" ))

# run dev.off() if problem
# plotting correlations heatmap
pheatmap::pheatmap(
  vsd_cor, 
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation", 
  annotation = mdta[,c(5,6)], 
  labels_row = mdta$sample, 
  labels_col = mdta$sample,
  annotation_colors = ann_colors,
  angle_col = 45, 
  main ="Correlation")

# transpose vst count matrix
vsdt <- t(vsd)

#calculating the variance of the vst transformed counts
var_genes <- apply(vsdt, 2, var)   #kind of a for loop, calculates the variance in all columns

# sorting & picking top 500 highly variable genes for plotting
select_var <- names(
  sort(
    var_genes, 
    decreasing=TRUE))[1:500]

# extracting 500 highly variable transcripts from vst matrix
vsdpc <- vsdt[,select_var]   # selecting subset based on the content of select_var

#principle component analysis
pc <- dudi.pca(
  vsdpc, 
  center=T, 
  scale=F, 
  scannf=F, 
  nf=6)
summary(pc) #Projected inertia (%) is the variance explained by each pc axis

#visualize the % of variance explained by each principle components
fviz_eig(pc, addlabels = T)

#adding metadata for principle components
pcaData <- as.data.frame(pc$li[,1:6])

pcaData$group <- FtDEG$samples$group1[match(
  rownames(pcaData), 
  rownames(FtDEG$samples))]

pcaData$cultivar <- FtDEG$samples$cultivar[match(
  rownames(pcaData), 
  rownames(FtDEG$samples))]

pcaData$treatment <- FtDEG$samples$treatment[match(
  rownames(pcaData), 
  rownames(FtDEG$samples))]

# plotting 
ggplot(
  pcaData, 
  aes(
    x = Axis1, 
    y = Axis2, 
    color = treatment, 
    shape = cultivar )) + 
  geom_point(size =3) +
  xlab("PC1: 43.8% variance") +
  ylab("PC2: 17.3% variance") +
  coord_fixed() +
  ggtitle("PCA (Top 500 highly variable genes)")

# model matrix
design <- model.matrix(
  ~ 0+group1,
  FtDEG$samples) # using grouped variable (treatment & cultivar)

design1 <- model.matrix(
  ~ 0+cultivar+treatment+cultivar:treatment,
  FtDEG$samples) # adding 0 removes intercept

colnames(design) <- levels(FtDEG$samples$group1)

# estimating common dispersion
FtDEG <- estimateGLMCommonDisp(FtDEG, design = design)  
FtDEG <- estimateGLMTrendedDisp(FtDEG, design = design)
FtDEG <- estimateGLMTagwiseDisp(FtDEG, design = design)

#biological coefficient of variance
sqrt(FtDEG$common.dispersion) #0.59 (is acceptable because the samples are from a population. but if samples are from feks kloner eller fra samme individ its too high)

plotBCV(FtDEG)

# making contrasts
my_contrasts<- makeContrasts(
  eT_eK = ET-EK,           # DEG due to treatment in engmo
  gT_gK = GT-GK,           # DEG due to treatment in grindstad
  egT_egK=(ET-EK)-(GT-GK), # genes which responded differently to treatment(cold) in cultivars 
  eK_gK = EK-GK,
  totT_totK = (ET+GT)-(EK+GK),
  levels = design)

# glm fit
fit <- glmFit(FtDEG, design)

# differentially expressed genes due to treatment in engmo
fitLRT <- glmLRT(
  fit,
  contrast = my_contrasts[,1])

summary(
  decideTests(
    fitLRT,
    p.value = 0.05))

# differentially expressed genes due to treatment in engmo
fit1 <- glmTreat(
  fit, 
  contrast = my_contrasts[,1],
  lfc = log2(1.2)) # should be down or up regulated by about 30 %. the ones below threshold are not important, but not all of the ones above are nececarily important. 

summary(
  decideTests(
    fit1,
    p.value = 0.05))

decideTests(fit1,p.value = 0.05)

eT_eK <- topTags(
  fit1, 
  n = nrow(FtDEG), 
  adjust.method = "fdr",
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

write.table(
  eT_eK,
  file = "Engmo_Trt_crtl.txt", 
  col.names = T, 
  row.names = F, 
  quote = F)

plotMD(fit1)

# differentially expressed genes due to treatment in Grindstad
fit2 <- glmTreat(
  fit,
  contrast = my_contrasts[,2],
  lfc = log2(1.2))

summary(
  decideTests(
    fit2,
    p.value = 0.05))

gT_gK <- topTags(
  fit2, 
  n = nrow(FtDEG), 
  adjust.method = "fdr", 
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

write.table(
  gT_gK,
  file = "Grindstad_Trt_ctrl.txt", 
  col.names = T, 
  row.names = F)

plotMD(fit2)

# genes which responded differently to treatment(cold) in engmo and grindstad
fit3 <- glmTreat(
  fit,
  contrast = my_contrasts[,3],
  lfc = log2(1.2))

summary(
  decideTests(
    fit3,
    p.value = 0.05))

egT_egK <- topTags(
  fit3, 
  n = nrow(FtDEG), 
  adjust.method = "fdr", 
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

write.table(
  egT_egK,
  file = "engmo_grindsatd_Trt_crtl.txt", 
  col.names = T , 
  row.names = F)

plotMD(fit3)

# genes which responded differently between controls in engmo and grindsmo
fit4 <- glmTreat(
  fit, 
  contrast = my_contrasts[,4],
  lfc = log2(1.2))

summary(
  decideTests(
    fit4,
    p.value = 0.05))

eK_gK <- topTags(
  fit4, 
  n = nrow(FtDEG), 
  adjust.method = "fdr",
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

write.table(
  eK_gK,file = "Engmo_Grindsatd_crtl.txt", 
  col.names = T, 
  row.names = F)

plotMD(fit4)

# total effect of treatment for both species on DEGs
fit5 <- glmTreat(
  fit, 
  contrast = my_contrasts[,5],
  lfc = log2(1.2))

summary(
  decideTests(
    fit5,
    p.value = 0.05))

totT_totK <- topTags(
  fit5, 
  n=nrow(FtDEG), 
  adjust.method = "fdr",
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

write.table(
  totT_totK,file = "total_trt_ctrl.txt", 
  col.names = T, 
  row.names = F)

 plotMD(fit4)

# universal genes responding to cold treatment
sum(eT_eK$GeneID %in% gT_gK$GeneID)

# venn plot
x= list(
  "ET vs EK" = eT_eK$GeneID, 
  "GT vs GK"= gT_gK$GeneID, 
  "ET-EK vs GT-GX" = egT_egK$GeneID)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  show_percentage = F,
  fill_alpha = 0.4,
  set_name_size = 6,
  stroke_color = "black",
  stroke_linetype = 0)

# upset plot
ven <- venndetail(
  list("ET vs EK" = eT_eK$GeneID, 
       "GT vs GK" = gT_gK$GeneID,
       "ET-EK vs GT-GX" = egT_egK$GeneID, 
       "EK vs GK" = eK_gK$GeneID, 
       "T vs K" = totT_totK$GeneID))

plot(ven, type = "upset")




