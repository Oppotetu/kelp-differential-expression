rm(list=ls())

library(tidyverse)
library(tidyr)
library(data.table)
library(readxl)
library(biomaRt)
library(qdapRegex)
library(clusterProfiler)
library(topGO)
library(WGCNA)

library(ggplot2)
library(dplyr)
library(stringi)
library(readr)
library(forcats)

setwd("C:/Users/simon/Prosjekter/RNA_seq/r_analysis")

# loading DGE files
fls = list.files(pattern = "*_out.txt", full.names = F)

dfs<-lapply(fls, FUN=read.table, sep="", header=T)

names(dfs)<- substr(fls,1, nchar(fls)-4)

## eggnog annotation file
# egg_anno<- read_xlsx("", col_names = T, sheet = 1, skip = 2 ) ##KEGG 

# blast annotation file
blastlist1<- read_tsv("./tsv_matches/matches_1.tsv", col_names = F)  ## ensemble ids

names(blastlist1) <- c('qseqid','sseqid','bitscore','evalue','pident','salltitles')

blastlist3<- read_tsv("matches_3.tsv", col_names = F)  ## ensemble ids

names(blastlist3) <- c('qseqid','sseqid','bitscore','evalue','pident','salltitles')

blastlist9<- read_tsv("matches_9.tsv", col_names = F)  ## ensemble ids

names(blastlist9) <- c('qseqid','sseqid','bitscore','evalue','pident','salltitles')

blastlist40<- read_tsv("matches_40.tsv", col_names = F)  ## ensemble ids

names(blastlist40) <- c('qseqid','sseqid','bitscore','evalue','pident','salltitles')

blastlist100<- read_tsv("matches_100.tsv", col_names = F)  ## ensemble ids

names(blastlist100) <- c('qseqid','sseqid','bitscore','evalue','pident','salltitles')

blastlist250<- read_tsv("matches_250.tsv", col_names = F)  ## ensemble ids

names(blastlist250) <- c('qseqid','sseqid','bitscore','evalue','pident','salltitles')

blastlistcut<- read_tsv("matches_cut.tsv", col_names = F)  ## ensemble ids

names(blastlistcut) <- c('qseqid','sseqid','bitscore','evalue','pident','salltitles')

#######################################

# extract ensembl_id and to a new column
blastlist1$ensembl_id<-c(ex_between(blastlist1$salltitles, "gene:", "transcript:"))

# extract gene_id and to a new column
blastlist1$geneid<-c(ex_between(blastlist1$salltitles, "GeneID=", "]"))

# corset cluster file
ccluster<- read.table("clusters.txt", header = F)
names(ccluster) <- c('Trinity_ID','Cluster_ID')

# attach corset cluster to transcripts to transfer annotations
egg_anno <- merge(egg_anno, ccluster, by.x = "query", by.y = "Trinity_ID")
blastlist1 <- merge(blastlist1, ccluster, by.x= "qseqid",by.y="Trinity_ID")
blastlist3 <- merge(blastlist3, ccluster, by.x= "qseqid",by.y="Trinity_ID")

# picking best hit based on evalue within corset cluster
kegg_anno<-setDT(egg_anno)[, .SD[which.min(evalue)], by="Cluster_ID"]
ensembl_anno<-setDT(blastlist1)[, .SD[which.min(evalue)], by="qseqid"]
ncbi_anno<-setDT(blastlist1)[, .SD[which.min(evalue)], by="qseqid"]

## Retrive Go annotations from biomart

listMarts(host="https://plants.ensembl.org") 

mart <- useMart("plants_mart", host="https://plants.ensembl.org")
#listDatasets(mart) # use this to list all data sets available in plants mart

searchDatasets(mart, pattern = c("Latissima")) # displays plants species in mart matching the query pattern

Taesmart<-useMart("plants_mart", dataset="taestivum_eg_gene",host="https://plants.ensembl.org") # extract mart

#listAttributes(mart = Taesmart) #lists all available attributes in the database
#searchAttributes(mart=Taesmart, pattern = "Entrez") # search of specific attributes

options(timeout = max(100000, getOption("timeout")))
Gomap <- getBM(attributes=c("ensembl_gene_id","go_id"),
               mart=Taesmart,uniqueRows=TRUE)

#geneID2GO <- split(Gomap$go_id, Gomap$ensembl_gene_id)

## if timeout issues persist use https://plants.ensembl.org/biomart/martview/ to manually download the dataset with necessary attributes

# manually downloaded GO annotations file
Go_anno<- read_tsv("Taesmart_GO.txt", col_names = T, na = "NA")
geneID2GO1 <- base::split(Go_anno$`GO term accession`, Go_anno$`Gene stable ID`)
geneID2GO <- lapply(geneID2GO1, unique)

#transfer annotations to differentially expressed genes

dfs2<-list()

for(i in 1:length(dfs)){
  df<- as.data.frame(dfs[[i]])
  name<- names(dfs)[i]
  df$Tri_ID<- ccluster$Trinity_ID[match(df$GeneID, ccluster$Cluster_ID)]# adds Trinity id
  df$ensembl_id<-ensembl_anno$ensembl_id[match(df$GeneID, ensembl_anno$Cluster_ID)]
  df$GID<-ncbi_anno$geneid[match(df$GeneID, ncbi_anno$Cluster_ID)]
  # df$ko<-kegg_anno$KEGG_ko[match(df$GeneID , kegg_anno$Cluster_ID)]
  rownames(df)<- df$GeneID
  print(name)
  dfs2[[name]]<- df
}

allgenes<- c(unique(ensembl_anno$ensembl_id)) ### All genes under study

go_dfs<-list()

##########    INSERTED    ###########


annotated1 <- read_tsv("./annotated1.tsv")
light1 <- read.table("./light1_out.txt", header = T)

colnames(annotated1) <- c("qseqid", "sseqid", "bitscore", "evalue", "pident", "goids", "count", "terms", "category")

anno_omit <- annotated1 %>% drop_na(goids)

allkeys <- read.table("clusters_new.txt", header = F)
colnames(allkeys) <- c("GeneID","TrinID")

temp1 <- left_join(allkeys,light1,by=c("GeneID"), keep = F)
temp1_omit <- temp1 %>% drop_na(PValue)

colnames(temp1_omit)[2] <- "qseqid"

temp2 <- left_join(temp1_omit,anno_omit,by=c("qseqid"), keep=F)
temp2_omit <- temp2 %>% drop_na(goids)
temp2_omit$log10.pval<--log10(as.numeric(temp2_omit$PValue))
temp2_omit$r_log10.pval <- round(temp2_omit$log10.pval, digit = 2)

anno_means <- temp2_omit %>%
  group_by(terms) %>% 
  dplyr::summarize(Mean = mean(log10.pval, na.rm=T)) %>%
  ungroup()
    
temp3_omit <- left_join(temp2_omit,anno_means, by=c("terms"), keep= F)
temp3_omit$rmeans <- round(temp3_omit$Mean, digit = 2)

anno_grouped <- temp3_omit %>%
  distinct(terms, .keep_all = T) %>%
  group_by(category) %>% # group your data based on the variable Rating
  slice_max(order_by = rmeans, n= 10) %>%
  ungroup() %>%
  arrange(category, desc(rmeans), terms)# order in descending order the variable One_Year_Return
  

ggplot(anno_grouped, aes(x = rmeans, y = fct_reorder(terms, sort(rmeans,decreasing = T)), fill = category)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ category, scales = "free", nrow = 3) +
  geom_text(
    aes(label = rmeans),
    color = "black",
    hjust = -0.1,
    size = 4,
    position = position_dodge(0.9)
  ) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    # axis.line = element_line(colour = "black"), 
    strip.text.x = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
    plot.title=element_text(hjust=0.5)
  )

# theme(text = element_text(size = 40), 
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        # panel.background = element_blank(), 
        # axis.line = element_line(colour = "black"), 
        # plot.title=element_text(hjust=0.5))

## kanskje
ggplot(anno_grouped, aes(x=reorder(terms, count), y=count)) +
  geom_bar(stat='identity', position=position_dodge()) +
  labs(x='Description', y='Pvalue') +
  coord_flip()

## for REVIGO(bookmarks) men bruker det nok ikke
goterms <- annotated1[,c("goids","evalue")]
newgos <- na.omit(goterms)
write.table(
  newgos, file="revigo1.txt",
  col.names = T,
  row.names = F, 
  sep=" ",
  quote = F
)

## kanskje, funka ikke helt 
anno_grouped %>%
  extract(terms, into = c("Term", "id"), "(.*)\\s\\((.*)\\)") %>%
  mutate(Term=stri_trans_totitle(Term)) %>%
  ggplot(aes(count, reorder(Term,count), fill = category)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ category, nrow = 3, , scales = "free") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank()
  ) +
  geom_text(
    aes(label = paste(count, ", Adj.P= ", bitscore)),
    color = "black",
    size = 4,
    hjust=-0.1,
    position = position_dodge(0.9)
  )

## funka heller ikke helt
p <- ggplot(data = anno_grouped, aes(x = fct_reorder(terms, sort(Position,decreasing = T)), y = count, fill = category)) +
  ## Plot "Description" in the x-axis following the order stated in the "Position" column
  ## vs normalized "Number of genes" in the second y-axis
  # geom_col(data = df6grouped, aes(x = fct_reorder(terms, desc(Position)), y = `Number of genes`/normalizer)) +
  ## Add a second y-axis based on the transformation of "Percentage of genes" to "Number of genes".
  ## Notice that the transformation undoes the normalization for the earlier geom_col.
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*normalizer, name = "Number of genes")) +
  ## Modify the aesthetic of the theme
  theme(axis.text.x = element_text(angle = 70, hjust = 1), axis.title.y = element_text(size = 8),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8),
        legend.key.size =  unit(0.2, "in"), plot.title = element_text(size = 11, hjust = 0.5)) +
  ## Add a title to the plot
  labs(x = NULL, title = "Gene Ontology (GO) Annotation")
p

##########    INSERTED    ###########

for(i in 1:length(dfs2)){
  name<-names(dfs2)[i]
  print(name)
  myInterestingGenes <- allgenes[allgenes %in% dfs2[[name]]$ensembl_id] #DE genes
  geneList <- factor(as.integer(allgenes %in% myInterestingGenes)) 
  names(geneList) <- allgenes
  #str(geneList)
  #table(geneList)
  GOdata_BP<-new("topGOdata", description=name, ontology="BP", allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO1, nodeSize = 10)
  GOdata_CC<-new("topGOdata", description=name, ontology="CC", allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO1,nodeSize = 10)
  GOdata_MF<-new("topGOdata", description=name, ontology="MF", allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO1,nodeSize = 10)
  #
  resultFisher_BP <- runTest(GOdata_BP,algorithm = "classic", statistic = "fisher")
  resultFisher_CC <- runTest(GOdata_CC,algorithm = "classic", statistic = "fisher")
  resultFisher_MF <- runTest(GOdata_MF,algorithm = "classic", statistic = "fisher")
  
  allgo_BP<-usedGO(GOdata_BP)
  allgo_CC<-usedGO(GOdata_CC)
  allgo_MF<-usedGO(GOdata_MF)
  #
  tab_BP <- GenTable(GOdata_BP, classic = resultFisher_BP, topNodes=length(allgo_BP))
  tab_BP$category<-"BP"
  tab_BP$pval<-str_remove_all(tab_BP$classic, "<")
  tab_BP$p.adj<-p.adjust(tab_BP$pval, method = "fdr")
  tab_BP$log10.pval<--log10(as.numeric(tab_BP$pval))
  #
  tab_CC <- GenTable(GOdata_CC, classic = resultFisher_CC, topNodes=length(allgo_CC))
  tab_CC$category<-"CC"
  tab_CC$pval<-str_remove_all(tab_CC$classic, "<")
  tab_CC$p.adj<-p.adjust(tab_CC$pval, method = "fdr")
  tab_CC$log10.pval<--log10(as.numeric(tab_CC$pval))
  #
  tab_MF <- GenTable(GOdata_MF, classic = resultFisher_MF, topNodes=length(allgo_MF))
  tab_MF$category<-"MF"
  tab_MF$pval<-str_remove_all(tab_MF$classic, "<")
  tab_MF$p.adj<-p.adjust(tab_MF$pval, method = "fdr")
  tab_MF$log10.pval<--log10(as.numeric(tab_MF$pval))
  #
  go<-rbind(tab_BP,tab_CC,tab_MF)
  sig_go<-subset(go, go$pval<=0.01)
  go_dfs[[name]]<-sig_go
}


ggplot(go_dfs) +
  aes()
       
       # [["eT_eK_"]] %>%
       #   group_by(Term, category) %>%
       #   arrange(dplyr::desc(category),
       #           dplyr::desc(log10.pval)) %>%
       #   dplyr::group_by(category)%>%
       #   dplyr::slice(1:8) %>%
       #   mutate(Term = factor(Term,levels = rev(Term))),
       # aes(log10.pval,Term, fill=category)) + 
       # 
  
  geom_bar(stat='identity', position='dodge') + 
  geom_text(aes(label = Annotated),
            color = "black", hjust = -0.1,
            size = 2.5, position = position_dodge(0.9)) +
  scale_x_continuous(position = "top") +
  xlab("-log10(pval)") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    #axis.text.x = element_blank(), ## remove text on x axis
    #axis.ticks = element_blank(), ## to remove ticks on x and y axis labels
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    strip.background = element_blank()) +
  theme(legend.key.size = unit(0.4, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=8)) + #change legend text font size
  ggtitle("light1")+theme(plot.title = element_text(vjust = - 8, hjust = -0.5 ))

##KEGG analysis

goi<-na.omit(gsub("-", NA,dfs2[["gT_gK_"]]$ko))
##remove ko: from ko ids
goi2<- substring(goi, regexpr(":", goi) + 1, nchar(goi))
#pathwayanalysis
xx=enrichKEGG(goi2, organism = "ko" , keyType="kegg",  pvalueCutoff = 0.05,pAdjustMethod = "fdr")

dotplot(xx, showCategory=10)+ ggtitle(paste("gT_gK","Enriched pathways",sep="    "))+theme(legend.key = element_blank())

