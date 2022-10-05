# Safekelp RNA-seq analysis

rm(list=ls())

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
library(dplyr)
library(tibble)

setwd("C:/Users/simon/Prosjekter/RNA_seq/r_analysis")

#loding counts matrix file
dta <- read.table("counts.txt", header = T)
dim(dta)
names(dta)

#creating metadata for the counts
mdta <- data.frame(
  id=c(
    'K1','K2','K3','K4','K5','K6','K7','K8','K10','K11','K12','K13','K14','K15','K16','K17',
    'K18','K19','K20','K21','K22','K23','K24','K25','K26','K27','K28','K29','K30','K31',
    'K32','K33','K34','K35','K36','K37','K38','K39','K40','K41','K42','K43','K44'),
  rep=c(
    1,1,1,1,2,2,2,2,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,
    8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11),
  Light=c(
    'MIN','MIN','MIN','MIN','UC','UC','UC','UC','MIN','MIN','MIN','MIN','MIN','MIN','MIN','MIN',
    'MIN','MIN','MIN','MED','MED','MED','MED','MED','MED','MED','MED','MED','MED','MED','MED',
    'MAX','MAX','MAX','MAX','MAX','MAX','MAX','MAX','MAX','MAX','MAX','MAX'),
  Days=c(
    'C0','C0','C0','C0','C0','C0','C0','C0','D1','D1','D1','D3','D3','D3','D3',
    'D9','D9','D9',
    'D9','D1','D1','D1','D1','D3','D3','D3','D3','D9','D9','D9','D9','D1','D1','D1','D1',
    'D3','D3','D3','D3','D9','D9','D9','D9'),
  group2=c(
    'C0','C0','C0','C0','UC','UC','UC','UC','C1','C1','C1','C3','C3','C3',
    'C3','C9','C9','C9','C9','MED1','MED1','MED1','MED1','MED3',
    'MED3','MED3','MED3','MED9','MED9','MED9','MED9','MAX1','MAX1','MAX1',
    'MAX1','MAX3','MAX3','MAX3','MAX3','MAX9','MAX9','MAX9','MAX9'),
  treatment2=c(
    'C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C',
    'C','C','C',
    'T','T','T','T','T','T','T','T','T','T','T','T','T','T','T','T',
    'T','T','T','T','T','T','T','T'),
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
# plotDensities(
#   cpm(
#     FtDEG$counts,
#     log = T),
#   # col = c("red", "blue", "green", "yellow",
#   #         "pink", "orange","black", "gray"),
#   legend = F,
#   main = "Density of reads")

##Filtering transcripts with low read counts
keep <- filterByExpr(
  FtDEG,
  min.count=20, ## BYTTE TIL 30?
  group = FtDEG$samples$group2)

table(keep) 

FtDEG <- FtDEG[keep, ,keep.lib.sizes=FALSE]
dim(FtDEG)

##plotting distribution of counts 2
# plotDensities(
#   log(FtDEG$counts),
#   # col = c("red", "blue", "green", "yellow",
#   #         "pink", "orange", "black", "gray"),
#   legend = F, 
#   main = "Density of reads after filtering")

#mean lib.size
lib.size <- as.character(
  round(mean(FtDEG$samples$lib.size * 1e-6), 
        2))

#plotting library sizes   # INKLUDERT 
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

#plotting transcripts counts    # IKKE INKLUDERT
# barplot(
#   rowSums(FtDEG$counts*1e-6), 
#   las=2, 
#   main="Counts per transcript", 
#   axisnames = FALSE, 
#   ylab = "counts in millions", 
#   cex.axis=0.8)

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

## Variance Stablizing transformation  - VST
vsd <- vst(
  round(FtDEG$counts), 
  blind = F) #default blind=T(for totally unsupervised clustering)

## Hierarchical Clustering 
sampleTree <- hclust(
  dist(
    t(vsd)))

# plot

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
# ann_colors = list(
#   cultivar = c(Engmo = "orange2",
#                Grindstad= "darkgreen"),
#   treatment = c(C = "#084594", T1 = "#6BAED6", T2 =  ))

# run dev.off() if problem
# plotting correlations heatmap

pheatmap::pheatmap(
  vsd_cor,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  annotation = mdta[,c(3,4)],
  labels_row = mdta$Days,
  labels_col = mdta$Light,
  # annotation_colors = ann_colors,
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

# pcaData$group <- FtDEG$samples$group2[match(
#   rownames(pcaData), 
#   rownames(FtDEG$samples))]

pcaData$Light <- FtDEG$samples$Light[match(
  rownames(pcaData), 
  rownames(FtDEG$samples))]

pcaData$Days <- FtDEG$samples$Days[match(
  rownames(pcaData),
  rownames(FtDEG$samples))]

# pcaData$treatment <- FtDEG$samples$treatment[match(
#   rownames(pcaData), 
#   rownames(FtDEG$samples))]

#plotting
ggplot(
  pcaData,
  aes(
    x = Axis1,
    y = Axis2,
    color = Light,
    shape = Days)) +
  geom_point(size =6) +
  xlab("PC1: 43.8% variance") +
  ylab("PC2: 17.3% variance") +
  coord_fixed() +
  ggtitle("PCA (top 500 highly variable genes)") +
  theme(text = element_text(size = 20),
        axis.text = element_text(face="bold", color="black", size=21),
        panel.border = element_rect(size = 2, linetype = "solid", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title=element_text(hjust=0.5))

# model matrix
design <- model.matrix(
  ~ 0+group2,
  FtDEG$samples) # using grouped variable

# design_alt <- model.matrix(
#   ~ 0+Light+Days+Light:Days,
#   FtDEG$samples) # adding 0 removes intercept

colnames(design) <- levels(FtDEG$samples$group2)

# estimating common dispersion
# FtDEG <- estimateGLMCommonDisp(FtDEG, design = design)  
# FtDEG <- estimateGLMTrendedDisp(FtDEG, design = design)
# FtDEG <- estimateGLMTagwiseDisp(FtDEG, design = design)

FtDEG <- estimateDisp(FtDEG, design = design)

#biological coefficient of variance
sqrt(FtDEG$common.dispersion) 

# = 0.49

# COMMON DISP across ALL samples ? 

# plotBCV(FtDEG, main = "Biological coefficient of variation")

# making contrasts
# my_contrasts<- makeContrasts(
#   MAX_MIN_1 = D1MAX-D1MIN,
#   MAX_MIN_3 = D3MAX-D3MIN,
#   MAX_MIN_9 = D9MAX-D9MIN,
#   
#   MIN_9_0 = D9MIN-C0, 
#   MED_9_1 = D9MED-D1MED, 
#   MAX_9_1 = D9MAX-D1MAX, 
#   
#   CUT = C0-NC,
#   
#   MED_MIN_1 = D1MED-D1MIN,
#   MED_MIN_3 = D3MED-D3MIN,
#   MED_MIN_9 = D9MED-D9MIN,
#   levels = design)

my_contrasts<- makeContrasts(
  MAX_MIN_1 = MAX1-C1,
  MAX_MIN_3 = MAX3-C3,
  MAX_MIN_9 = MAX9-C9,
  
  MIN_9_0 = C9-C0, 
  MED_9_1 = MED9-MED1, 
  MAX_9_1 = MAX9-MAX1, 
  
  CUT = C0-UC,
  
  levels = design)


# glm fit
fit <- glmFit(FtDEG, design)

# MAX1-C1

fit1 <- glmTreat(
  fit, 
  contrast = my_contrasts[,1],
  lfc = log2(1.5)) # should be down or up regulated by about 30 %. the ones below threshold are not important, 
                  # but not all of the ones above are nececarily important. 

summary(
  decideTests(
    fit1,
    p.value = 0.05))

light1 <- topTags(
  fit1, 
  n = nrow(FtDEG), 
  adjust.method = "fdr",
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

# write.table(
#   light1,
#   file = "light1_outny.txt",
#   col.names = T,
#   row.names = F,
#   quote = F)


# plotMD(fit1, main = "Mean-Difference plot, MAX1 vs C1")
#        # main.cex=2, cex.main=2, cex.lab=1.5, cex.axis=1.5)


# MAX3-C3

fit2 <- glmTreat(
  fit,
  contrast = my_contrasts[,2],
  lfc = log2(1.5))

summary(
  decideTests(
    fit2,
    p.value = 0.05))

light3 <- topTags(
  fit2, 
  n = nrow(FtDEG), 
  adjust.method = "fdr", 
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

# write.table(
#   light3,
#   file = "light3_outny.txt",
#   col.names = T,
#   row.names = F)

# plotMD(fit2, main = "Mean-Difference plot, MAX3 vs C3")

# MAX9-C9

fit3 <- glmTreat(
  fit,
  contrast = my_contrasts[,3],
  lfc = log2(1.5))

summary(
  decideTests(
    fit3,
    p.value = 0.05))

light9 <- topTags(
  fit3, 
  n = nrow(FtDEG), 
  adjust.method = "fdr", 
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

# write.table(
#   light9,
#   file = "light9_outny.txt",
#   col.names = T ,
#   row.names = F)

# plotMD(fit3, main = "Mean-Difference plot, MAX9 vs C9")

# C9-C0

fit4 <- glmTreat(
  fit, 
  contrast = my_contrasts[,4],
  lfc = log2(1.5))

summary(
  decideTests(
    fit4,
    p.value = 0.05))

timeMIN <- topTags(
  fit4, 
  n = nrow(FtDEG), 
  adjust.method = "fdr",
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

# write.table(
#   timeMIN,file = "time40_outny.txt",
#   col.names = T,
#   row.names = F)

# plotMD(fit4, main = "Mean-Difference plot, C9 vs C0")

# MED9-MED1

fit5 <- glmTreat(
  fit, 
  contrast = my_contrasts[,5],
  lfc = log2(1.5))

summary(
  decideTests(
    fit5,
    p.value = 0.05))

timeMED <- topTags(
  fit5, 
  n=nrow(FtDEG), 
  adjust.method = "fdr",
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

# write.table(
#   timeMED,file = "time100_outny.txt",
#   col.names = T,
#   row.names = F)

# plotMD(fit5, main = "Mean-Difference plot, MED9 vs MED1")

# MAX9-MAX1

fit6 <- glmTreat(
  fit, 
  contrast = my_contrasts[,6],
  lfc = log2(1.5))

summary(
  decideTests(
    fit6,
    p.value = 0.05))

timeMAX <- topTags(
  fit6, 
  n=nrow(FtDEG), 
  adjust.method = "fdr",
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

# write.table(
#   timeMAX,file = "time250_outny.txt",
#   col.names = T,
#   row.names = F)

# plotMD(fit6, main = "Mean-Difference plot, MAX9 vs MAX1")

# C0-UC

fit7 <- glmTreat(
  fit, 
  contrast = my_contrasts[,7],
  lfc = log2(1.5))

summary(
  decideTests(
    fit7,
    p.value = 0.05))

cut_nocut <- topTags(
  fit7, 
  n=nrow(FtDEG), 
  adjust.method = "fdr",
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

# write.table(
#   cut_nocut,file = "cut_nocut_outny.txt",
#   col.names = T,
#   row.names = F)

# plotMD(fit7, main = "Mean-Difference plot, C0 vs UC")

# MED1-C1

fit8 <- glmTreat(
  fit, 
  contrast = my_contrasts[,8],
  lfc = log2(1.5))

summary(
  decideTests(
    fit8,
    p.value = 0.05))

lightMED1 <- topTags(
  fit8, 
  n = nrow(FtDEG), 
  adjust.method = "fdr",
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

# write.table(
#   lightMED1,
#   file = "lightMED1_outny.txt",
#   col.names = T,
#   row.names = F,
#   quote = F)

# plotMD(fit8, main = "Mean-Difference plot, MED1 vs C1")

# MED3-C3

fit9 <- glmTreat(
  fit, 
  contrast = my_contrasts[,9],
  lfc = log2(1.5)) 

summary(
  decideTests(
    fit9,
    p.value = 0.05))

lightMED3 <- topTags(
  fit9, 
  n = nrow(FtDEG), 
  adjust.method = "fdr",
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

# write.table(
#   lightMED3,
#   file = "lightMED3_outny.txt",
#   col.names = T,
#   row.names = F,
#   quote = F)

# plotMD(fit9, main = "Mean-Difference plot, MED3 vs C3")

# MED9-C9

fit10 <- glmTreat(
  fit, 
  contrast = my_contrasts[,10],
  lfc = log2(1.5))

summary(
  decideTests(
    fit10,
    p.value = 0.05))

lightMED9 <- topTags(
  fit10, 
  n = nrow(FtDEG), 
  adjust.method = "fdr",
  p.value = 0.05)$table %>% rownames_to_column("GeneID")

# write.table(
#   lightMED9,
#   file = "lightMED9_outny.txt",
#   col.names = T,
#   row.names = F,
#   quote = F)

# plotMD(fit10, main = "Mean-Difference plot, MED9 vs C9")

dev.off()

# title(ylab="log-fold-change",mgp=c(2, 1, 0))

par(mfrow = c(1,3), cex = 0.78, cex.main = 2, cex.lab = 2)
plotMD(fit1, main = "MAX1 vs C1", ylab = '')
title(ylab="log-fold-change", mgp = c(2.4,1,0))
plotMD(fit2, main = "MAX3 vs C3", ylab = '')
plotMD(fit3, main = "MAX9 vs C9", ylab = '')

# par(mfrow = c(1,3), cex = 1, cex.main = 2, cex.lab = 2)
plotMD(fit4, main = "C9 vs C0", ylab = '')
title(ylab="log-fold-change", mgp = c(2.4,1,0))
plotMD(fit5, main = "MED9 vs MED1", ylab = '')
plotMD(fit6, main = "MAX9 vs MAX1", ylab = '')

# par(mfrow = c(1,4), cex = 0.78, cex.main = 2, cex.lab = 2)
plotMD(fit7, main = "C0 vs UC", ylab = '')
title(ylab="log-fold-change", mgp = c(2.4,1,0))
plotMD(fit8, main = "MED1 vs C1", ylab = '')
plotMD(fit10, main = "MED9 vs C9", ylab = '')

# plotMD(fit9, main = "MED3 vs C3", ylab = '')

# venn plot
x= list(
  "MAX9 vs MAX1" = timeMAX$GeneID,
  "MAX1 vs C1" = light1$GeneID,
  "MED9 vs MED1" = timeMED$GeneID)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  show_percentage = F,
  fill_alpha = 0.4,
  set_name_size = 6,
  stroke_color = "black",
  stroke_linetype = 0)

# upset plot

upset <-   list(
  "MAX1 vs C1" = light1$GeneID,
  "MAX3 vs C3" = light3$GeneID,
  "MAX9 vs C9" = light9$GeneID, 
  "C9 vs C0" = timeMIN$GeneID,
  "MED9 vs MED1" = timeMED$GeneID, 
  "MAX9 vs MAX1" = timeMAX$GeneID, 
  "C0 vs UC" = cut_nocut$GeneID)

upset(fromList(upset), nsets = 7, 
      mainbar.y.label = "DE genes in intersection", 
      sets.x.label = "Number of genes \nin the dataset",
      sets.bar.color = c("red", "green", "blue", "orange", "yellow", "lightblue", "purple" ),
      point.size = 3.5,
      line.size = 1.2,
      text.scale = c(2.8,2.8,2.5,2.5,2.5,2.5))

# c(intersection size title, 
#   intersection size tick labels, 
#   set size title, 
#   set size tick labels, 
#   set names, 
#   numbers above bars)

upset2 <- list(
  "UC" = round(rowMeans(FtDEG[["counts"]][,4:8])),
  "C0" = round(rowMeans(FtDEG[["counts"]][,1:4])),
  "C1" = round(rowMeans(FtDEG[["counts"]][,9:11])),
  "C3" = round(rowMeans(FtDEG[["counts"]][,12:15])),
  "C9" = round(rowMeans(FtDEG[["counts"]][,16:19])),
  
  "MED1" = round(rowMeans(FtDEG[["counts"]][,20:23])),
  "MED3" = round(rowMeans(FtDEG[["counts"]][,24:27])),
  "MED9" = round(rowMeans(FtDEG[["counts"]][,28:31])),
  
  "MAX1" = round(rowMeans(FtDEG[["counts"]][,32:35])),
  "MAX3" = round(rowMeans(FtDEG[["counts"]][,36:39])),
  "MAX9" = round(rowMeans(FtDEG[["counts"]][,40:43])))

upset(fromList(upset2), nsets = 11, 
      mainbar.y.label = "Genes in intersection", 
      sets.x.label = "Number of genes \nin the dataset", 
      sets.bar.color = c("red", "green", "blue", "orange", "yellow", "lightblue", "purple", "darkgreen", 
                         "darkblue", "pink", "aquamarine"),
      point.size = 3.5,
      line.size = 1.2,
      text.scale = c(2.8,2.8,2.5,2.5,2,2.5))








