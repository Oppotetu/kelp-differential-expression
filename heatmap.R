rm(list= ls())

library(tidyverse)
library(tidyr)
library(readxl)
library(qdapRegex)
library(clusterProfiler)
library(topGO)
library(ggplot2)
library(dplyr)
library(forcats)
library(stringr)
library(readr)
library(forcats)
library(stringr)
library(edgeR)
library(stats)
library(gplots)
library(tibble)
library(RColorBrewer)
library(BBmisc)

## import
setwd("C:/Users/simon/Prosjekter/RNA_seq/r_analysis")

## counts
counts <- read.table("counts.txt", header = T)
counts <- rownames_to_column(counts,var = "GeneID")
counts <- as_tibble(counts)

counts$C0 <- rowMeans(counts[,2:5])
counts$UC <- rowMeans(counts[,6:9])
counts$C1 <- rowMeans(counts[,10:12])
counts$C3 <- rowMeans(counts[,13:16])
counts$C9 <- rowMeans(counts[,17:20])
counts$MED1 <- rowMeans(counts[,21:24])
counts$MED3 <- rowMeans(counts[,25:28])
counts$MED9 <- rowMeans(counts[,29:32])
counts$MAX1 <- rowMeans(counts[,33:36])
counts$MAX3 <- rowMeans(counts[,37:40])
counts$MAX9 <- rowMeans(counts[,41:44])

counts <- counts %>% 
  select(-c(2:44)) 

## contrast files
fls <- list.files(pattern = "*_outny.txt", full.names = F)
dfs <- lapply(fls, FUN=read.table, sep="", header=T)
dfs <- lapply(dfs, FUN=as_tibble, sep="", header=T)
dfs <- lapply(dfs[], FUN=select, -c(3,5:6))
names(dfs)<- substr(fls,1, nchar(fls)-10)
list2env(dfs,envir=.GlobalEnv)
colnames(light1) <- c('GeneID', 'logFC_1', 'logCPM_1')
colnames(light3) <- c('GeneID', 'logFC_3', 'logCPM_3')
colnames(light9) <- c('GeneID', 'logFC_9', 'logCPM_9')
colnames(time40) <- c('GeneID', 'logFC_40', 'logCPM_40')
colnames(time100) <- c('GeneID', 'logFC_100', 'logCPM_100')
colnames(time250) <- c('GeneID', 'logFC_250', 'logCPM_250')
colnames(cut_nocut) <- c('GeneID', 'logFC_cut', 'logCPM_cut')

contrasts <- full_join(light1, light3, by = c("GeneID"), keep = F)
contrasts <- full_join(contrasts, light9, by = c('GeneID'), keep = F)
contrasts <- full_join(contrasts, time40, by = c('GeneID'), keep = F)
contrasts <- full_join(contrasts, time100, by = c('GeneID'), keep = F)
contrasts <- full_join(contrasts, time250, by = c('GeneID'), keep = F)
contrasts <- full_join(contrasts, cut_nocut, by = c('GeneID'), keep = F)

rm(list= c('light1','light3','light9','time40','time100','time250','cut_nocut'))

## keys
keys <- read.table("clusters_new.txt", header = F)
colnames(keys) <- c("GeneID","TrinID")

## blast results
afls <- list.files(pattern = "*_3.tsv", full.names = F)
adfs <- lapply(afls, FUN=read.table, sep="\t", header=T, na.strings = '-')
adfs <- lapply(adfs, FUN=as_tibble, sep="", header=T)
adfs <- lapply(adfs[], FUN=select, -c(3:9))
adfs <- lapply(adfs[], FUN=distinct, .keep_all = T)
names(adfs)<- substr(afls,1, nchar(afls)-6)
list2env(adfs,envir=.GlobalEnv)
colnames(annotated1) <- c('TrinID', 'sseqid')
colnames(annotated3) <- c('TrinID', 'sseqid')
colnames(annotated9) <- c('TrinID', 'sseqid')
colnames(annotated40) <- c('TrinID', 'sseqid')
colnames(annotated100) <- c('TrinID', 'sseqid')
colnames(annotated250) <- c('TrinID', 'sseqid')
colnames(annotatedcutno) <- c('TrinID', 'sseqid')

annotations <- full_join(annotated1, annotated3, by = c("TrinID", 'sseqid'), keep = F)
annotations <- full_join(annotations, annotated9, by = c('TrinID', 'sseqid'), keep = F)
annotations <- full_join(annotations, annotated40, by = c('TrinID', 'sseqid'), keep = F)
annotations <- full_join(annotations, annotated100, by = c('TrinID', 'sseqid'), keep = F)
annotations <- full_join(annotations, annotated250, by = c('TrinID', 'sseqid'), keep = F)
annotations <- full_join(annotations, annotatedcutno, by = c('TrinID', 'sseqid'), keep = F)

rm(list= c('annotated1','annotated3','annotated9','annotated40','annotated100','annotated250','annotatedcutno'))

annotations <- annotations %>% 
  drop_na(sseqid)

annotations <- filter(annotations, !(str_detect(annotations$sseqid[], 'Uncharacterized')))

for (i in 1:length(annotations$sseqid)) {
  if (str_detect(annotations$sseqid[i], 'Uncharacterized')) {
    next
  } else { 
    annotations$sseqid[i] <- strsplit(annotations$sseqid[i], "[_]", fixed = F)[[1]][2]
    annotations$sseqid[i] <- strsplit(annotations$sseqid[i], "OS", fixed = F)[[1]][1]
  }
}

# annogrep <- annotations %>%
#   filter(grepl('peroxidase|cellular oxidant', sseqid))

# annogrep <- left_join(annogrep, keys, by = "TrinID")
# annogrep <- left_join(annogrep, contrasts, by = "GeneID")

annotations <- inner_join(annotations, keys, by = c('TrinID'), keep = F)

annotations <- annotations %>%
  distinct(GeneID, .keep_all = T)


################ annotation and counts ################




anno_counts <- inner_join(annotations, counts, by = c('GeneID'), keep = F)

anno_medians0 <- anno_counts %>%
  group_by(sseqid) %>% 
  dplyr::summarize(C0 = max(C0, na.rm=T)) %>%
  ungroup()
anno_medians1 <- anno_counts %>%
  group_by(sseqid) %>% 
  dplyr::summarize(UC = max(UC, na.rm=T)) %>%
  ungroup()
anno_medians2 <- anno_counts %>%
  group_by(sseqid) %>% 
  dplyr::summarize(C1 = max(C1, na.rm=T)) %>%
  ungroup()
anno_medians3 <- anno_counts %>%
  group_by(sseqid) %>% 
  dplyr::summarize(C3 = max(C3, na.rm=T)) %>%
  ungroup()
anno_medians4 <- anno_counts %>%
  group_by(sseqid) %>% 
  dplyr::summarize(C9 = max(C9, na.rm=T)) %>%
  ungroup()
anno_medians5 <- anno_counts %>%
  group_by(sseqid) %>% 
  dplyr::summarize(MED1 = max(MED1, na.rm=T)) %>%
  ungroup()
anno_medians6 <- anno_counts %>%
  group_by(sseqid) %>% 
  dplyr::summarize(MED3 = max(MED3, na.rm=T)) %>%
  ungroup()
anno_medians7 <- anno_counts %>%
  group_by(sseqid) %>% 
  dplyr::summarize(MED9 = max(MED9, na.rm=T)) %>%
  ungroup()
anno_medians8 <- anno_counts %>%
  group_by(sseqid) %>% 
  dplyr::summarize(MAX1 = max(MAX1, na.rm=T)) %>%
  ungroup()
anno_medians9 <- anno_counts %>%
  group_by(sseqid) %>% 
  dplyr::summarize(MAX3 = max(MAX3, na.rm=T)) %>%
  ungroup()
anno_medians10 <- anno_counts %>%
  group_by(sseqid) %>% 
  dplyr::summarize(MAX9 = max(MAX9, na.rm=T)) %>%
  ungroup()

medians <- full_join(anno_medians0, anno_medians1, by = c('sseqid'), keep = F)
medians <- full_join(medians, anno_medians2, by = c('sseqid'), keep = F)
medians <- full_join(medians, anno_medians3, by = c( 'sseqid'), keep = F)
medians <- full_join(medians, anno_medians4, by = c( 'sseqid'), keep = F)
medians <- full_join(medians, anno_medians5, by = c( 'sseqid'), keep = F)
medians <- full_join(medians, anno_medians6, by = c( 'sseqid'), keep = F)
medians <- full_join(medians, anno_medians7, by = c( 'sseqid'), keep = F)
medians <- full_join(medians, anno_medians8, by = c( 'sseqid'), keep = F)
medians <- full_join(medians, anno_medians9, by = c( 'sseqid'), keep = F)
medians <- full_join(medians, anno_medians10, by = c( 'sseqid'), keep = F)


rm(list=c('anno_medians0','anno_medians1','anno_medians2','anno_medians3','anno_medians4',
          'anno_medians5','anno_medians6','anno_medians7','anno_medians8','anno_medians9',
          'anno_medians10'))

catmat <- medians[,-c(3,7:9)]




################ annotation and contrasts ################



# anno_contrasts <- inner_join(annotations, contrasts, by = c('GeneID'), keep = F)
# 
# anno_medians0 <- anno_contrasts %>%
#   group_by(sseqid) %>% 
#   dplyr::summarize(MAX1_C1_FC = max(logFC_1, na.rm=T)) %>%
#   ungroup()
# anno_medians1 <- anno_contrasts %>%
#   group_by(sseqid) %>% 
#   dplyr::summarize(MAX3_C3_FC = max(logFC_3, na.rm=T)) %>%
#   ungroup()
# anno_medians2 <- anno_contrasts %>%
#   group_by(sseqid) %>% 
#   dplyr::summarize(MAX9_C9_FC = max(logFC_9, na.rm=T)) %>%
#   ungroup()
# anno_medians3 <- anno_contrasts %>%
#   group_by(sseqid) %>% 
#   dplyr::summarize(C9_C0_FC = max(logFC_40, na.rm=T)) %>%
#   ungroup()
# anno_medians4 <- anno_contrasts %>%
#   group_by(sseqid) %>% 
#   dplyr::summarize(MED9_MED1_FC = max(logFC_100, na.rm=T)) %>%
#   ungroup()
# anno_medians5 <- anno_contrasts %>%
#   group_by(sseqid) %>% 
#   dplyr::summarize(MAX9_MAX1_FC = max(logFC_250, na.rm=T)) %>%
#   ungroup()
# anno_medians6 <- anno_contrasts %>%
#   group_by(sseqid) %>% 
#   dplyr::summarize(C0_UC_FC = max(logFC_cut, na.rm=T)) %>%
#   ungroup()
# 
# anno_medians7 <- anno_contrasts %>%
#   group_by(sseqid) %>% 
#   dplyr::summarize(MAX1_C1_CPM = max(logCPM_1, na.rm=T)) %>%
#   ungroup()
# anno_medians8 <- anno_contrasts %>%
#   group_by(sseqid) %>% 
#   dplyr::summarize(MAX3_C3_CPM = max(logCPM_3, na.rm=T)) %>%
#   ungroup()
# anno_medians9 <- anno_contrasts %>%
#   group_by(sseqid) %>% 
#   dplyr::summarize(MAX9_C9_CPM = max(logCPM_9, na.rm=T)) %>%
#   ungroup()
# anno_medians10 <- anno_contrasts %>%
#   group_by(sseqid) %>% 
#   dplyr::summarize(C9_C0_CPM = max(logCPM_40, na.rm=T)) %>%
#   ungroup()
# anno_medians11 <- anno_contrasts %>%
#   group_by(sseqid) %>% 
#   dplyr::summarize(MED9_MED1_CPM = max(logCPM_100, na.rm=T)) %>%
#   ungroup()
# anno_medians12 <- anno_contrasts %>%
#   group_by(sseqid) %>% 
#   dplyr::summarize(MAX9_MAX1_CPM = max(logCPM_250, na.rm=T)) %>%
#   ungroup()
# anno_medians13 <- anno_contrasts %>%
#   group_by(sseqid) %>% 
#   dplyr::summarize(C0_UC_CPM = max(logCPM_cut, na.rm=T)) %>%
#   ungroup()
# 
# medians2 <- full_join(anno_medians0, anno_medians1, by = c('sseqid'), keep = F)
# medians2 <- full_join(medians2, anno_medians2, by = c('sseqid'), keep = F)
# medians2 <- full_join(medians2, anno_medians3, by = c( 'sseqid'), keep = F)
# medians2 <- full_join(medians2, anno_medians4, by = c( 'sseqid'), keep = F)
# medians2 <- full_join(medians2, anno_medians5, by = c( 'sseqid'), keep = F)
# medians2 <- full_join(medians2, anno_medians6, by = c( 'sseqid'), keep = F)
# medians2 <- full_join(medians2, anno_medians7, by = c( 'sseqid'), keep = F)
# medians2 <- full_join(medians2, anno_medians8, by = c( 'sseqid'), keep = F)
# medians2 <- full_join(medians2, anno_medians9, by = c( 'sseqid'), keep = F)
# medians2 <- full_join(medians2, anno_medians10, by = c( 'sseqid'), keep = F)
# medians2 <- full_join(medians2, anno_medians11, by = c( 'sseqid'), keep = F)
# medians2 <- full_join(medians2, anno_medians12, by = c( 'sseqid'), keep = F)
# medians2 <- full_join(medians2, anno_medians13, by = c( 'sseqid'), keep = F)
# 
# rm(list=c('anno_medians0','anno_medians1','anno_medians2','anno_medians3','anno_medians4',
#           'anno_medians5','anno_medians6','anno_medians7','anno_medians8','anno_medians9',
#           'anno_medians10','anno_medians11','anno_medians12','anno_medians13'))
# 
# medians2sel <- medians2[,c(1,2,6,7,9,13,14)]
# 
# 
# ## catmat ##
# 
# catmat <- inner_join(medians1sel, medians2sel, by = c('sseqid'), keep = F)



mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# col.cell <- c("purple","orange")[sampleinfo$CellType]


################ plotting ################

catmat1 <- column_to_rownames(catmat, var = "sseqid")
catmatcpm <- cpm(catmat1[1:7], log = T)
var_genes <- apply(catmatcpm, 1, var)

select_var <- names(sort(var_genes, decreasing=TRUE))[1:10]

highly_variable_cpm <- catmatcpm[select_var,]

#write.table(highly_variable_cpm, file = "cpm1_ut.tsv", sep = '\t')

highly_variable_cpm <- read_xlsx("cpm1_inn.xlsx")

highly_variable_cpm <- column_to_rownames(highly_variable_cpm, var = "...1")

highly_variable_cpm <- as.matrix(highly_variable_cpm)

normcpm <- normalize(
  highly_variable_cpm,
  method = "range",
  range = c(0,1),
  margin = 1,
  on.constant = "quiet")

pheatmap::pheatmap(
  normcpm,
  #scale = "column",
  # clustering_distance_rows = "correlation",
  # clustering_distance_cols = "correlation",
  # cluster_rows = T,
  # cluster_cols = T,
  # annotation = mdta[,c(4,5)],
  # labels_row = mdta$Days,
  # labels_col = mdta$Light,
  annotation_colors = rev(morecols(50)),
  angle_col = 45,
  fontsize = 20,
  main ="Top 10 most variable genes across samples")

################ chosen genes ################

catgrep <- catmat %>%
  filter(grepl('vanadium|Vanadium|xanthin|glutathione|Glutathione', sseqid))

test1 <- anno_counts %>%
  distinct(sseqid, .keep_all = T)

test2 <- inner_join(catgrep, test1, by = "sseqid", keep = F)

test3 <- test2[,c(1,10)]

test4 <- inner_join(test3, contrasts, by = "GeneID")

catcpm2 <- column_to_rownames(catgrep, var = "sseqid")

#write.table(catcpm2, file = "catgrep_ut3.tsv", sep = '\t')

catgrep_inn <- read_xlsx("catgrep_inn3.xlsx")

catgrep_inn <- column_to_rownames(catgrep_inn, var = "...1")

normcat <- cpm(catgrep_inn[1:7], log = T)

normcat <- normalize(
  catgrep_inn,
  method = "standardize",
  range = c(0,1),
  margin = 1,
  on.constant = "quiet"
)
colnames(normcat) <- colnames(catgrep_inn)


highly_variable_cpm <- as.matrix(normcat)

pheatmap::pheatmap(
  highly_variable_cpm,
  scale = "row",
  #clustering_distance_rows = "correlation",
  # clustering_method = "complete",
  # clustering_distance_cols = "correlation",
  #cluster_rows = T,
  # cluster_cols = T,
  # annotation = mdta[,c(4,5)],
  # labels_row = mdta$Days,
  # labels_col = mdta$Light,
  annotation_colors = rev(morecols(50)),
  angle_col = 45,
  fontsize = 20,
  main ="Genes associated with vHPOs and heavy metal")


