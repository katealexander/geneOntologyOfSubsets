# Rationale
The goal of this analysis is to determine whether subsets of transcription factor target genes belong to different functional categories. For example, we found that only a subset of HIF2A target genes experience HIF2A-regulated DNA-speckle association, and we wanted to know whether these particular targets genes were functionally distinct from the target genes that do not operate via speckles. 

# Strategy
To get an idea of which trancription factor functions are differential between two subgroups of target genes, I first got the functional categories of all the transcription factor target genes, without subsetting. Then, I calculated the Gene Ontology statististics of the target gene subsets for those GO terms. 

# Implementation
The below Rscript uses clusterProfiler to get Gene Ontology terms. It uses three gene lists:
1. A full list of HIF2A target genes: "downInPT2399_RNA.txt"
2. The HIF2A target genes that have decreasing SON: "786O_PT2399downRNA_SONdown.txt"
3. The HIF2A target genes that do not have decreasing SON: "786O_PT2399downRNA_SONns.txt"
```Rscript GOsubsetCompare_HIF2Atargets.R```
The above script will output dotplots of the Gene Ontology terms that are specific to the genes that decrease SON "GO_Down.pdf", the target genes that do not decrease SON "GO_notDown.pdf", and the terms that are shared between the two subgroups "GO_shared.pdf".
