# GO CHARTS

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
library(stringr)
library(qdapRegex)
library(tibble)

setwd("C:/Users/simon/Prosjekter/RNA_seq/r_analysis")

###############################    1    ###############################
rm(list=ls())

allkeys <- read.table("clusters_new.txt", header = F)
colnames(allkeys) <- c("GeneID","TrinID")

annotated1 <- read_tsv("./annotated1_3.tsv")
light1 <- read.table("./light1_outny.txt", header = T)

colnames(annotated1) <- c("qseqid", "sseqid", "bitscore", "evalue", "pident",
                          "goids", "count", "terms", "category")

temp1 <- left_join(allkeys,light1,by=c("GeneID"), keep = F)

colnames(temp1)[2] <- "qseqid"

temp2 <- left_join(temp1,annotated1,by=c("qseqid"), keep=F)

str_which(temp2$terms, "peroxidase")
str_which(temp2$terms, "cellular oxidant")

write.table(temp2[c(324819, 167756, 167772, 198757, 198778, 248844, 248851,
                    248857, 324822, 335874, 335899, 335925, 
                    335950, 335975, 336000, 336025),],
            file= "MAX1vsC1_print.txt", sep = "\t")

temp2_omit <- temp2 %>% drop_na()
temp2_omit$log10.pval<--log10(as.numeric(temp2_omit$PValue))
temp2_omit$r_log10.pval <- round(temp2_omit$log10.pval, digit = 2)

for (i in 1:length(temp2_omit$terms)) {
  if (isTRUE(str_detect(temp2_omit$terms[i], ","))) {
    temp2_omit$terms[i] <- strsplit(temp2_omit$terms[i], "[,]", fixed = F)[[1]][1]
  } else {
    next
  }
}


anno_means <- temp2_omit %>%
  group_by(terms) %>% 
  dplyr::summarize(Mean = mean(log10.pval, na.rm=T)) %>%
  ungroup()

temp3_omit <- left_join(temp2_omit,anno_means, by=c("terms"), keep= F)
temp3_omit$rmeans <- round(temp3_omit$Mean, digit = 2)

anno_grouped <- temp3_omit %>%
  distinct(terms, .keep_all = T) %>%
  group_by(category) %>% 
  slice_max(order_by = rmeans, n= 10) %>%
  ungroup()



ggplot(anno_grouped, aes(x = rmeans, y = fct_reorder(terms, rmeans), fill = category)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ category, scales = "free_y", nrow = 3) +
  geom_text(
    aes(label = rmeans),
    color = "black",
    hjust = -0.1,
    size = 5.5,
    position = position_dodge(0.9)
  ) +
  ggtitle("Distribution of GO terms for MAX1 vs C1") +
  xlab("-log10(p values)") + 
  ylab("Top GO terms grouped by category") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    panel.background = element_blank(), 
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 20, face = "bold"),
    strip.background = element_blank(),
    text = element_text(size = 19), 
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 20))

###############################    3    ############################### 
rm(list=ls())

allkeys <- read.table("clusters_new.txt", header = F)
colnames(allkeys) <- c("GeneID","TrinID")

annotated1 <- read_tsv("./annotated3_3.tsv")
light1 <- read.table("./light3_outny.txt", header = T)

colnames(annotated1) <- c("qseqid", "sseqid", "bitscore", "evalue",
                          "pident", "goids", "count", "terms", "category")

temp1 <- left_join(allkeys,light1,by=c("GeneID"), keep = F)

colnames(temp1)[2] <- "qseqid"

temp2 <- left_join(temp1,annotated1,by=c("qseqid"), keep=F)

str_which(temp2$terms, "peroxidase")
str_which(temp2$terms, "cellular oxidant")

write.table(temp2[c(248083, 248092, 240454, 240458, 
                    240460, 248086, 248095),], file = "MAX3vsC3_print.txt",
            sep = "\t")

temp2_omit <- temp2 %>% drop_na()
temp2_omit$log10.pval<--log10(as.numeric(temp2_omit$PValue))
temp2_omit$r_log10.pval <- round(temp2_omit$log10.pval, digit = 2)


for (i in 1:length(temp2_omit$terms)) {
  if (isTRUE(str_detect(temp2_omit$terms[i], ","))) {
    temp2_omit$terms[i] <- strsplit(temp2_omit$terms[i], "[,]", fixed = F)[[1]][1]
  } else {
    next
  }
}



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
  ungroup() 

ggplot(anno_grouped, aes(x = rmeans, y = fct_reorder(terms, rmeans), fill = category)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ category, scales = "free_y", nrow = 3) +
  geom_text(
    aes(label = rmeans),
    color = "black",
    hjust = -0.1,
    size = 5.5,
    position = position_dodge(0.9)) +
  ggtitle("Distribution of GO terms for MAX3 vs C3") +
  xlab("-log10(p values)") + 
  ylab("Top GO terms grouped by category") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    panel.background = element_blank(), 
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 20, face = "bold"),
    strip.background = element_blank(),
    text = element_text(size = 19), 
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 20))

###############################    9    ###############################
rm(list=ls())

allkeys <- read.table("clusters_new.txt", header = F)
colnames(allkeys) <- c("GeneID","TrinID")

annotated1 <- read_tsv("./annotated9_3.tsv")
light1 <- read.table("./light9_outny.txt", header = T)

colnames(annotated1) <- c("qseqid", "sseqid", "bitscore", "evalue", 
                          "pident", "goids", "count", "terms", "category")

temp1 <- left_join(allkeys,light1,by=c("GeneID"), keep = F)

colnames(temp1)[2] <- "qseqid"

temp2 <- left_join(temp1,annotated1,by=c("qseqid"), keep=F)

str_which(temp2$terms, "peroxidase")
str_which(temp2$terms, "cellular oxidant")

write.table(temp2[c(300057, 240387, 240390, 240392, 300058),], 
            file = "MAX9vsC9_print.txt", sep = "\t")

temp2_omit <- temp2 %>% drop_na()
temp2_omit$log10.pval<--log10(as.numeric(temp2_omit$PValue))
temp2_omit$r_log10.pval <- round(temp2_omit$log10.pval, digit = 2)



for (i in 1:length(temp2_omit$terms)) {
  if (isTRUE(str_detect(temp2_omit$terms[i], ","))) {
    temp2_omit$terms[i] <- strsplit(temp2_omit$terms[i], "[,]", fixed = F)[[1]][1]
  } else {
    next
  }
}



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
  ungroup()

ggplot(anno_grouped, aes(x = rmeans, y = fct_reorder(terms, rmeans), fill = category)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ category, scales = "free_y", nrow = 3) +
  geom_text(
    aes(label = rmeans),
    color = "black",
    hjust = -0.1,
    size = 5.5,
    position = position_dodge(0.9)
  ) +
  ggtitle("Distribution of GO terms for MAX9 vs C9") +
  xlab("-log10(p values)") + 
  ylab("Top GO terms grouped by category") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    panel.background = element_blank(), 
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 20, face = "bold"),
    strip.background = element_blank(),
    text = element_text(size = 19), 
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 20))

###############################    40    ###############################
rm(list=ls())

allkeys <- read.table("clusters_new.txt", header = F)
colnames(allkeys) <- c("GeneID","TrinID")

annotated1 <- read_tsv("./annotated40_3.tsv")
light1 <- read.table("./time40_outny.txt", header = T)

colnames(annotated1) <- c("qseqid", "sseqid", "bitscore", "evalue", 
                          "pident", "goids", "count", "terms", "category")

temp1 <- left_join(allkeys,light1,by=c("GeneID"), keep = F)

colnames(temp1)[2] <- "qseqid"

temp2 <- left_join(temp1,annotated1,by=c("qseqid"), keep=F)

str_which(temp2$terms, "peroxidase")
str_which(temp2$terms, "cellular oxidant")

write.table(temp2[c(180689, 180695),], 
            file = "C9vsC0_print.txt", sep = "\t")

temp2_omit <- temp2 %>% drop_na()
temp2_omit$log10.pval<--log10(as.numeric(temp2_omit$PValue))
temp2_omit$r_log10.pval <- round(temp2_omit$log10.pval, digit = 2)

for (i in 1:length(temp2_omit$terms)) {
  if (isTRUE(str_detect(temp2_omit$terms[i], ","))) {
    temp2_omit$terms[i] <- strsplit(temp2_omit$terms[i], "[,]", fixed = F)[[1]][1]
  } else {
    next
  }
}

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
  ungroup() 

ggplot(anno_grouped, aes(x = rmeans, y = fct_reorder(terms, rmeans), fill = category)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ category, scales = "free_y", nrow = 3) +
  geom_text(
    aes(label = rmeans),
    color = "black",
    hjust = -0.1,
    size = 5.5,
    position = position_dodge(0.9)
  ) +
  ggtitle("Distribution of GO terms for C9 vs C0") +
  xlab("-log10(p values)") + 
  ylab("Top GO terms grouped by category") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    panel.background = element_blank(), 
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 20, face = "bold"),
    strip.background = element_blank(),
    text = element_text(size = 19), 
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 20))

###############################    100    ###############################
rm(list=ls())

allkeys <- read.table("clusters_new.txt", header = F)
colnames(allkeys) <- c("GeneID","TrinID")

annotated1 <- read_tsv("./annotated100_3.tsv")
light1 <- read.table("./time100_outny.txt", header = T)

colnames(annotated1) <- c("qseqid", "sseqid", "bitscore", "evalue", 
                          "pident", "goids", "count", "terms", "category")

temp1 <- left_join(allkeys,light1,by=c("GeneID"), keep = F)

colnames(temp1)[2] <- "qseqid"

temp2 <- left_join(temp1,annotated1,by=c("qseqid"), keep=F)

str_which(temp2$terms, "peroxidase")
str_which(temp2$terms, "cellular oxidant")

write.table(temp2[c(182060, 182071, 250157, 335946, 335953, 338097, 338122, 
                    338147, 338172, 338197),], 
            file = "MED9vsMED1_print.txt", sep = "\t")

temp2_omit <- temp2 %>% drop_na()
temp2_omit$log10.pval<--log10(as.numeric(temp2_omit$PValue))
temp2_omit$r_log10.pval <- round(temp2_omit$log10.pval, digit = 2)

for (i in 1:length(temp2_omit$terms)) {
  if (isTRUE(str_detect(temp2_omit$terms[i], ","))) {
    temp2_omit$terms[i] <- strsplit(temp2_omit$terms[i], "[,]", fixed = F)[[1]][1]
  } else {
    next
  }
}

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
  ungroup() 

ggplot(anno_grouped, aes(x = rmeans, y = fct_reorder(terms, rmeans), fill = category)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ category, scales = "free_y", nrow = 3) +
  geom_text(
    aes(label = rmeans),
    color = "black",
    hjust = -0.1,
    size = 5.5,
    position = position_dodge(0.9)
  ) +
  ggtitle("Distribution of GO terms for MED9 vs MED1") +
  xlab("-log10(p values)") + 
  ylab("Top GO terms grouped by category") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    panel.background = element_blank(), 
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 20, face = "bold"),
    strip.background = element_blank(),
    text = element_text(size = 19), 
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 20))

###############################    250    ###############################
rm(list=ls())

allkeys <- read.table("clusters_new.txt", header = F)
colnames(allkeys) <- c("GeneID","TrinID")

annotated1 <- read_tsv("./annotated250_3.tsv")
light1 <- read.table("./time250_outny.txt", header = T)

colnames(annotated1) <- c("qseqid", "sseqid", "bitscore", "evalue", 
                          "pident", "goids", "count", "terms", "category")

temp1 <- left_join(allkeys,light1,by=c("GeneID"), keep = F)

colnames(temp1)[2] <- "qseqid"

temp2 <- left_join(temp1,annotated1,by=c("qseqid"), keep=F)

str_which(temp2$terms, "peroxidase")
str_which(temp2$terms, "cellular oxidant")

write.table(temp2[c(258814, 258823, 365345, 365349,170931, 170947, 253110,
                    253117, 253123, 258817, 258826, 340090, 340115, 340141,
                    340166,340191, 340216, 340241, 340266, 365348, 365352),], 
            file = "MAX9vsMAX1_print.txt", sep = "\t")

temp2_omit <- temp2 %>% drop_na()
temp2_omit$log10.pval<--log10(as.numeric(temp2_omit$PValue))
temp2_omit$r_log10.pval <- round(temp2_omit$log10.pval, digit = 2)

for (i in 1:length(temp2_omit$terms)) {
  if (isTRUE(str_detect(temp2_omit$terms[i], ","))) {
    temp2_omit$terms[i] <- strsplit(temp2_omit$terms[i], "[,]", fixed = F)[[1]][1]
  } else {
    next
  }
}

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
  ungroup()

ggplot(anno_grouped, aes(x = rmeans, y = fct_reorder(terms, rmeans), fill = category)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ category, scales = "free_y", nrow = 3) +
  geom_text(
    aes(label = rmeans),
    color = "black",
    hjust = -0.1,
    size = 5.5,
    position = position_dodge(0.9)
  ) + 
  ggtitle("Distribution of GO terms for MAX9 vs MAX1") +
  xlab("-log10(p values)") + 
  ylab("Top GO terms grouped by category") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    panel.background = element_blank(), 
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 20, face = "bold"),
    strip.background = element_blank(),
    text = element_text(size = 19), 
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 20))

###############################    cut    ###############################
rm(list=ls())

allkeys <- read.table("clusters_new.txt", header = F)
colnames(allkeys) <- c("GeneID","TrinID")

annotated1 <- read_tsv("./annotatedcutno_3.tsv")
light1 <- read.table("./cut_nocut_outny.txt", header = T)

colnames(annotated1) <- c("qseqid", "sseqid", "bitscore", "evalue", "pident", "goids", "count", "terms", "category")

anno_omit <- annotated1 %>% drop_na(goids)

temp1 <- left_join(allkeys,light1,by=c("GeneID"), keep = F)
temp1_omit <- temp1 %>% drop_na(PValue)

colnames(temp1_omit)[2] <- "qseqid"

temp2 <- left_join(temp1_omit,anno_omit,by=c("qseqid"), keep=F)
temp2_omit <- temp2 %>% drop_na()
temp2_omit$log10.pval<--log10(as.numeric(temp2_omit$PValue))
temp2_omit$r_log10.pval <- round(temp2_omit$log10.pval, digit = 2)

str_which(temp2_omit$terms, "peroxidase")

for (i in 1:length(temp2_omit$terms)) {
  if (isTRUE(str_detect(temp2_omit$terms[i], ","))) {
    temp2_omit$terms[i] <- strsplit(temp2_omit$terms[i], "[,]", fixed = F)[[1]][1]
  } else {
    next
  }
}

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
  ungroup() 

ggplot(anno_grouped, aes(x = rmeans, y = fct_reorder(terms, rmeans), fill = category)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ category, scales = "free_y", nrow = 3) +
  geom_text(
    aes(label = rmeans),
    color = "black",
    hjust = -0.1,
    size = 5.5,
    position = position_dodge(0.9)
  ) +
  ggtitle("Distribution of GO terms for C0 vs UC") +
  xlab("-log10(p values)") + 
  ylab("Top GO terms grouped by category") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    panel.background = element_blank(), 
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 20, face = "bold"),
    strip.background = element_blank(),
    text = element_text(size = 19), 
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 20))



