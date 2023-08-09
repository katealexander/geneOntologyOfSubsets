setwd("~/Documents/geneOntologyOfSubsets/")

library(ggpubr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(hrbrthemes)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)


SONdown <- read.table("786O_PT2399downRNA_SONdown.txt")
SONdown <- SONdown$V1
SONns <- read.table("786O_PT2399downRNA_SONns.txt")
SONns <- SONns$V1
allHIF <- read.table("downInPT2399_RNA.txt")
allHIF <- allHIF$V1

SONdown.df <- bitr(SONdown, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
SONns.df <- bitr(SONns, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
allHIF.df <- bitr(allHIF, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db) 

all = enrichGO(gene = allHIF.df$ENTREZID,
             ont  = "BP",
             OrgDb = org.Hs.eg.db,
             pvalueCutoff = 0.01,
             pAdjustMethod = "fdr",
             minGSSize = 5,
             maxGSSize = 300, 
             qvalueCutoff = 0.01,
             readable = FALSE)
all <- clusterProfiler::simplify(all, cutoff = 0.67)
# Convert to data frame and set BP description as rownames
all = as.data.frame(all, row.names = all$Description)

down = enrichGO(gene = SONdown.df$ENTREZID,
               ont  = "BP",
               OrgDb = org.Hs.eg.db,
               pvalueCutoff = 1,
               pAdjustMethod = "fdr",
               minGSSize = 5,
               maxGSSize = 300, 
               qvalueCutoff = 1,
               readable = FALSE)
down = as.data.frame(down, row.names = down$Description)

ns = enrichGO(gene = SONns.df$ENTREZID,
                ont  = "BP",
                OrgDb = org.Hs.eg.db,
                pvalueCutoff = 1,
                pAdjustMethod = "fdr",
                minGSSize = 5,
                maxGSSize = 300, 
                qvalueCutoff = 1,
                readable = FALSE)
ns = as.data.frame(ns, row.names = ns$Description)

down$sample <- "down"
ns$sample <- "ns"

downSub <- down[down$Description %in% row.names(all),]
nsSub <- ns[ns$Description %in% row.names(all),]

merged <- merge(downSub, nsSub, by = "Description")
sharedTerms <- merged$Description[merged$qvalue.x < 0.05 & merged$qvalue.y < 0.05]
nsTerms <- merged$Description[merged$qvalue.x > 0.05 & merged$qvalue.y < 0.005]
downTerms <- merged$Description[merged$qvalue.x < 0.01 & merged$qvalue.y > 0.05]

combo <- rbind(nsSub, downSub)
combo$GeneRatio <- combo$Count/as.numeric(str_split(combo$GeneRatio, "/")[[1]][2])
comboSubDown <- combo[combo$Description %in% downTerms,]
comboSubShared <- combo[combo$Description %in% sharedTerms,]
comboSubNs <- combo[combo$Description %in% nsTerms,]

p = comboSubDown%>%
  ggplot(aes(x=sample,y=reorder(Description, -log(qvalue)))) +
  geom_point(aes(size= GeneRatio,color= -log(qvalue)) ) +
  #next line add border
  geom_point(aes(size= GeneRatio),shape = 1,colour = "black") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 2, limits = c(0,10), na.value = "red")+
  theme_classic() +
  theme(
    legend.position="top",
    legend.title = element_text(face = "bold", size = 10),
    legend.title.align = 2,
    axis.title.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12,colour = 'black', angle = 90),
    axis.text.y = element_text(size = 12,colour = 'black'),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  #if you want to change overall size of dot
  scale_size_continuous(range = c(2, 10), breaks = seq(0, 0.2, by = .05), limits = c(0,0.2))
filename = "GO_Down.pdf"
pdf(filename, width = 4.75, height = 5, onefile=FALSE)
print(p)
dev.off()

p = comboSubNs%>%
  ggplot(aes(x=sample,y=reorder(Description, -log(qvalue)))) +
  geom_point(aes(size= GeneRatio,color= -log(qvalue)) ) +
  #next line add border
  geom_point(aes(size= GeneRatio),shape = 1,colour = "black") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 2, limits = c(0,10), na.value = "red")+
  theme_classic() +
  theme(
    legend.position="top",
    legend.title = element_text(face = "bold", size = 10),
    legend.title.align = 2,
    axis.title.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12,colour = 'black', angle = 90),
    axis.text.y = element_text(size = 12,colour = 'black'),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  #if you want to change overall size of dot
  scale_size_continuous(range = c(2, 10), breaks = seq(0, 0.2, by = .05), limits = c(0,0.2))
filename = "GO_notDown.pdf"
pdf(filename, width = 5, height = 7, onefile=FALSE)
print(p)
dev.off()


p = comboSubShared%>%
  ggplot(aes(x=sample,y=reorder(Description, -log(qvalue)))) +
  geom_point(aes(size= GeneRatio,color= -log(qvalue)) ) +
  #next line add border
  geom_point(aes(size= GeneRatio),shape = 1,colour = "black") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 2, limits = c(0,10), na.value = "red")+
  theme_classic() +
  theme(
    legend.position="top",
    legend.title = element_text(face = "bold", size = 8),
    legend.title.align = 2,
    axis.title.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12,colour = 'black', angle = 90),
    axis.text.y = element_text(size = 12,colour = 'black'),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  #if you want to change overall size of dot
  scale_size_continuous(range = c(2, 10), breaks = seq(0, 0.2, by = .05), limits = c(0,0.2))
filename = "GO_shared.pdf"
pdf(filename, width = 6, height = 6.5, onefile=FALSE)
print(p)
dev.off()

