rm(list=ls())

library(stringr)
library(Bios2cor)
library(seqinr)
library(terra)

setwd("C:/Users/simon/Prosjekter/RNA_seq/r_analysis")

allkeys <- read.table("clusters_new.txt", header = F)

outfls = list.files(pattern = "*_outny.txt", full.names = F) 
dfs<-lapply(outfls, FUN=read.table, sep=" ", header=T) 
names(dfs)<- substr(outfls,1, nchar(outfls)-4)

list2env(dfs,envir=.GlobalEnv)

light1_trins <- allkeys[,1] %in% light1_outny[,1]
light3_trins <- allkeys[,1] %in% light3_outny[,1]
light9_trins <- allkeys[,1] %in% light9_outny[,1]
time40_trins <- allkeys[,1] %in% time40_outny[,1]
time100_trins <- allkeys[,1] %in% time100_outny[,1]
time250_trins <- allkeys[,1] %in% time250_outny[,1]
cutno_trins <- allkeys[,1] %in% cut_nocut_outny[,1]

light1_fin <- allkeys[light1_trins,2]
light3_fin <- allkeys[light3_trins,2]
light9_fin <- allkeys[light9_trins,2]
time40_fin <- allkeys[time40_trins,2]
time100_fin <- allkeys[time100_trins,2]
time250_fin <- allkeys[time250_trins,2]
cutno_fin <- allkeys[cutno_trins,2]

assembly_dta <- import.fasta(
  "safekelp_trinity_assembly.Trinity.fasta", 
  aa.to.upper = TRUE, gap.to.dash = TRUE, log.file = NULL)

light1_write <- assembly_dta[light1_fin]
light3_write <- assembly_dta[light3_fin]
light9_write <- assembly_dta[light9_fin]
time40_write <- assembly_dta[time40_fin]
time100_write <- assembly_dta[time100_fin]
time250_write <- assembly_dta[time250_fin]
cutno_write <- assembly_dta[cutno_fin]

file.create("light1_extracted4.fasta")
file.create("light3_extracted4.fasta")
file.create("light9_extracted4.fasta")
file.create("time40_extracted4.fasta")
file.create("time100_extracted4.fasta")
file.create("time250_extracted4.fasta")
file.create("cutno_extracted4.fasta")

write.fasta(light1_write, names = names(light1_write), open = "a",
            file.out = "light1_extracted4.fasta", nbchar = 80, as.string = FALSE)
write.fasta(light3_write, names = names(light3_write), open = "a",
            file.out = "light3_extracted4.fasta", nbchar = 80, as.string = FALSE)
write.fasta(light9_write, names = names(light9_write), open = "a",
            file.out = "light9_extracted4.fasta", nbchar = 80, as.string = FALSE)
write.fasta(time40_write, names = names(time40_write), open = "a",
            file.out = "time40_extracted4.fasta", nbchar = 80, as.string = FALSE)
write.fasta(time100_write, names = names(time100_write), open = "a",
            file.out = "time100_extracted4.fasta", nbchar = 80, as.string = FALSE)
write.fasta(time250_write, names = names(time250_write), open = "a",
            file.out = "time250_extracted4.fasta", nbchar = 80, as.string = FALSE)
write.fasta(cutno_write, names = names(cutno_write), open = "a",
            file.out = "cutno_extracted4.fasta", nbchar = 80, as.string = FALSE)








