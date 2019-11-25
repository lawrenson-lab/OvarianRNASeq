---
title: "Superenhancer analysis"
author: "Kate Lawrenson, Marcos A. S. Fonseca, Felipe Segato, et. al."
date: "November, 2018"
output: 
  html_document:
    smart: false
navbar:
  left:
    - text: "Home"
      href: index.html
    - text: "Normalization and DEG analysis"
      href: transc_paper_markdown.html
    - text: "One-class analysis"
      href: prj.mkd.transcriptome.OC.html
    - text: "Superenhancer analysis"
      href: se-analysis.html
    - text: "Cistrome analysis"
      href: tf-analysis.html

---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
#setwd("/var/www/html/OvarianRNAseq/")
#library(rmarkdown)
#render("se-analysis.Rmd")

```

## Introduction

In the following workflow, we describe the main procedures to analyse the gene expression profile of associated genes to superenghancers (SE) genomic regions.

```{r }
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
```
### Auxiliary functions
Obtaining only distinct peaks from a given `Genomic Ranges` Object

```{r }
distinctGR <- function(objectGR) {
  adj_hgsonc_union = objectGR
  
  index = 1
  final = length(adj_hgsonc_union)
  
  distinct_table = NULL
  
  while (index <= final) {
    
    res = as.data.frame(findOverlaps(adj_hgsonc_union[index], adj_hgsonc_union))
    n = dim(res)[1]
    
    if (n > 1) {
      adj_hgsonc_union = adj_hgsonc_union[-res$subjectHits[2:n]]
    }
    
    distinct_table = rbind(distinct_table, as.data.frame(adj_hgsonc_union[res$subjectHits[1]]))
    
    index = index + 1
    final = length(adj_hgsonc_union)
  }
  
  return(distinct_table)
  
}
```

Obtaining only distinct peaks by comparing two `Genomic Ranges` objects

```{r }
distinctGRPairs <- function(query, subject) {
  index = 1
  final = length(query)
  
  distinct_table = NULL
  
  while (index <= final) {
    
    res = as.data.frame(findOverlaps(query[index], subject))
    n = dim(res)[1]
    
    if (n < 1) {
      distinct_table = rbind(distinct_table, as.data.frame(query[index]))
    }
    
    index = index + 1
  }
  
  return(distinct_table)
  
}

```
### Part 1
### Data loading

The data is organized in a data frame with 387 samples (190 normals and 197 TCGA) as columns and 20928 unique exprerssed genes as rows, named as `combat_alln_mean`. Sample annotation is organized in another data frame `NormalSamplesTable`. The TCGA Sample annotation is also added in Normals sample table annotation as follow.

```{r load}
load("RData/190Normal-394TCGA-combatGCnorm_all_data.rda")

```

```{r echo=FALSE}
library(knitr)
kable(dados.RNAseq[1:5,1:10])
```

```{r }
load("RData/NormalsSamplesTableV3.RData")
head(NormalSamplesTable)

nTCGAsamples = 394

TCGAsamplesTable = data.frame(Samples = colnames(dados.RNAseq)[191:584], CellType = c(rep("HGSOC", nTCGAsamples)), Colors=c(rep("green", nTCGAsamples)), SP="None", Batch="None" )

SamplesTable = rbind(NormalSamplesTable, TCGAsamplesTable)

```

### Gene genomic annotation
We selected the genomic annotation, regarding chromosome name, start, end coordinates for human hg19 build. After the annation is merged with our expressed genes.

```{r results='hide'}
dados.RNAseq <- as.data.frame(dados.RNAseq)

txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- genes(txdb_hg19)
symbols <- unlist(mapIds(org.Hs.eg.db, genes$gene_id, "SYMBOL", "ENTREZID", multiVals = "first"))
genes <- as.data.frame(genes)
genes$symbol <- as.matrix(symbols)
genes <- genes[!(duplicated(genes$symbol)),]
head(genes)

dados.RNAseq <- merge(genes,dados.RNAseq, by.x="symbol",by.y="row.names")

```
Convert gene annotation into a Genomic Ranges object `dados.RNAseq.GR`. The merged data.frame also is splited into annotation `dados.RNAseq.DF` and expression information `dados.RNAseq.f`. After, the expression values are scaled as follows.

```{r }
dados.RNAseq.GR <- makeGRangesFromDataFrame(dados.RNAseq[,1:7],TRUE)
dados.RNAseq.DF <- dados.RNAseq[,1:7]
dados.RNAseq.f <- dados.RNAseq[,c(8:591)]

dados.RNAseq.scaled <- t(scale(t(dados.RNAseq.f)))
rownames(dados.RNAseq.scaled) <- dados.RNAseq.DF$symbol
```

### Superenhancer experiments
In this step we will processing the superenhancer bed files. First, the data is loaded and converted to Genomic Ranges Object. The experiments are combined in a single object. Then, using the aulixiliary function `distinctGR`, we selected only those regions that are unique considering the entire combined set. The same procedures are applied to OSEC, FTSEC and HGSOC superenhancers.

```{r }
########### OSEC #################
ose4 = read.table("se/IOSE4_ResultCount_C30H4ACXX_1_SHE1405A185.male.hg19.fa.nodup.superEnhancers.bed")
colnames(ose4) = c("chrom", "start", "end", "name", "score", "strand")
ose4 = makeGRangesFromDataFrame(ose4)

ose11 = read.table("se/IOSE11_ResultCount_D1AHBACXX_8_SHE1405A43.male.hg19.fa.nodup.superEnhancers.bed")
colnames(ose11) = c("chrom", "start", "end", "name", "score", "strand")
ose11 = makeGRangesFromDataFrame(ose11)

ose4.df <- as.data.frame(ose4)
ose11.df <- as.data.frame(ose11)

nose_total <- rbind(ose4.df, ose11.df)
nose_total <- makeGRangesFromDataFrame(nose_total)

distinct_table = distinctGR(nose_total)
nose_union = makeGRangesFromDataFrame(distinct_table)

########### FTSEC #################

ft33 = read.table("se/FT33_ResultCount_C26CUACXX_2_SHE1405A78.male.hg19.fa.nodup.superEnhancers.bed")
colnames(ft33) = c("chrom", "start", "end", "name", "score", "strand")
ft33 = makeGRangesFromDataFrame(ft33)

ft246 = read.table("se/FT246_ResultCount_C269LACXX_5_SHE1405A86.male.hg19.fa.nodup.superEnhancers.bed")
colnames(ft246)  = c("chrom", "start", "end", "name", "score", "strand")
ft246 = makeGRangesFromDataFrame(ft246)

ft33.df <- as.data.frame(ft33)
ft246.df <- as.data.frame(ft246)

nfte_total = rbind(ft33.df, ft246.df)
nfte_total = makeGRangesFromDataFrame(nfte_total)

distinct_table = distinctGR(nfte_total)
nfte_union = makeGRangesFromDataFrame(distinct_table)

########### HGSOC #################

hgs229 = read.table("se/HGSerous_229_IP.nodup.superEnhancers.bed")
colnames(hgs229) = c("chrom", "start", "end", "name", "score", "strand")
hgs229 = makeGRangesFromDataFrame(hgs229)

hgs429 = read.table("se/HGSerous_429_IP.nodup.superEnhancers.bed")
colnames(hgs429) = c("chrom", "start", "end", "name", "score", "strand")
hgs429 = makeGRangesFromDataFrame(hgs429)

hgs561 = read.table("se/HGSerous_561_IP.nodup.superEnhancers.bed")
colnames(hgs561) = c("chrom", "start", "end", "name", "score", "strand")
hgs561 = makeGRangesFromDataFrame(hgs561)

hgs550 = read.table("se/HG_Serous_550_IP.nodup.superEnhancers.bed")
colnames(hgs550) = c("chrom", "start", "end", "name", "score", "strand")
hgs550 = makeGRangesFromDataFrame(hgs550)

hgs270 = read.table("se/HGSerous_270_IP.nodup.superEnhancers.bed")
colnames(hgs270) = c("chrom", "start", "end", "name", "score", "strand")
hgs270 = makeGRangesFromDataFrame(hgs270)

hgs229.df <- as.data.frame(hgs229)
hgs429.df <- as.data.frame(hgs429)
hgs561.df <- as.data.frame(hgs561)
hgs550.df <- as.data.frame(hgs550)
hgs270.df <- as.data.frame(hgs270)

hgsoc_total <- rbind(hgs229.df, hgs429.df, hgs561.df, hgs550.df, hgs270.df)
hgsoc_total = makeGRangesFromDataFrame(hgsoc_total)

distinct_table = distinctGR(hgsoc_total)
hgsoc_union = makeGRangesFromDataFrame(distinct_table)
```
In this step the SE are cross compared in order to obtain only distict ones per cell type. The auxiliary function `distinctGRPairs` is used.

```{r }
subj = makeGRangesFromDataFrame(rbind(as.data.frame(hgsoc_union), as.data.frame(nfte_union)))
tempN = distinctGRPairs(nose_union, subj)
nose_distinct = makeGRangesFromDataFrame(tempN)

subj = makeGRangesFromDataFrame(rbind(as.data.frame(hgsoc_union), as.data.frame(nose_union)))
tempF = distinctGRPairs(nfte_union, subj)
nfte_distinct = makeGRangesFromDataFrame(tempF)

subj = makeGRangesFromDataFrame(rbind(as.data.frame(nose_union), as.data.frame(nfte_union)))
tempR = distinctGRPairs(hgsoc_union, subj)
hgsoc_distinct = makeGRangesFromDataFrame(tempR)

```
### Preparing the parameters
Identify expression profile considering all samples into 3 groups:
OSEC - 119 samples
FTSEC - 71 samples
HGSOC - 394 samples

We also defined a step size of 5000 bp and a window size of 1000000 bp.

```{r }
OSECLines <- which(SamplesTable$CellType == "NOSE")
FTSECLines <- which(SamplesTable$CellType == "NFTE")
HGSOCLines <- which(SamplesTable$CellType == "HGSOC")

step = 5000
window = 1000000
range.neigh = seq(-window, window, step)
```
### Part 2
Obtaining the gene expression profile across the 3 cell types considering OSEC specific enhancers.

```{r }
enhancer = nose_distinct
enh_name = "OSEC_enh_"
names <- paste0("OSEC_enh_", seq(1, length(enhancer)))
```
### Part 3
Now we can begin to find the genes and their correspondent expression in all 3 cell types with the following iteraction loop.

```{r }
#####################

e.index = 1
genes.range.all = NULL

while (e.index <= length(enhancer)) {
  #print(paste0("Enhancer ", e.index))
  
  e.info = enhancer[e.index]
  
  mid.pos = round(start(e.info) + ((end(e.info) - start(e.info)) / 2))
  s.pos = mid.pos - window
  e.pos = mid.pos + window
  chr.pos = as.character(seqnames(e.info))
  new.genes.range = seq(s.pos, e.pos, step)
  
  new.enh.range.GR = GRanges(seqnames=chr.pos, ranges=IRanges(s.pos, e.pos))
  
  res = findOverlaps(new.enh.range.GR, dados.RNAseq.GR, type = "any", select = "all", ignore.strand=TRUE)
  
  if (length(res) > 0) {
    
    genes.range = dados.RNAseq.DF[subjectHits(res),]
    
    OSEC.mean =  apply(dados.RNAseq.f[subjectHits(res),OSECLines],1,mean,na.rm=T)
    FTSEC.mean =  apply(dados.RNAseq.f[subjectHits(res),FTSECLines],1,mean,na.rm=T)
    HGSOC.mean =  apply(dados.RNAseq.f[subjectHits(res),HGSOCLines],1,mean,na.rm=T)
    
    genes.range = cbind(Enh = paste0(enh_name, e.index), genes.range, Mean.OSEC=OSEC.mean, Mean.FTSEC=FTSEC.mean, Mean.HGSOC=HGSOC.mean)
    
    for (gi in 1:nrow(genes.range)) {
      
      gene = genes.range[gi,]
      index = 1
      
      for (value in new.genes.range) {
        
        if (value > gene$start) {
          break
        }
        index = index + 1
      }
      
      genes.range$pos[gi] = range.neigh[index]
      
    } # end for
    
    genes.range.all = rbind(genes.range.all, genes.range)
    
  }
 
  e.index = e.index + 1 
}

#################
```
## Part 4
Calculate the z-score for each cell type gene list.

```{r }

info.plot.nose = data.frame(Position = range.neigh, Type = "OSEC", Value=NA)
info.plot.ftsec = data.frame(Position = range.neigh, Type = "FTSEC", Value=NA)
info.plot.hgsoc = data.frame(Position = range.neigh, Type = "HGSOC", Value=NA)

index.pos = 1
for ( pos in range.neigh ) {
  
  if (pos < 0) {
    sub.range = seq(pos, 0, step)
    
  } else if (pos > 0) {
    sub.range = seq(0, pos, step)
  }
  
  rows.m <- which(genes.range.all$pos %in% sub.range)
  enh.name <- unique(genes.range.all$Enh[rows.m])
  enh.index <- which(names %in% enh.name)
  genes.name = unique(genes.range.all$symbol[rows.m])
  
  info.plot.nose$Value[index.pos] <- mean(dados.RNAseq.scaled[genes.name, OSECLines])
  info.plot.ftsec$Value[index.pos] <- mean(dados.RNAseq.scaled[genes.name, FTSECLines])
  info.plot.hgsoc$Value[index.pos] <- mean(dados.RNAseq.scaled[genes.name, HGSOCLines])
  
  index.pos = index.pos + 1
}

##########################
```

Plotting the final result
```{r }
nose.info.plot.all = rbind(info.plot.nose, info.plot.ftsec, info.plot.hgsoc)

p_osec = ggplot(data=nose.info.plot.all, aes(x=Position, y=Value, colour = Type)) +
  labs(title = paste0("OSEC-specific SEs"), y = "Expression z-score", x = "Genomic positions") +
  theme(plot.title = element_text(size=12)) +
  scale_color_manual(values=c("#3333ffff", "#ff3333ff", "#33ff33ff" )) +
  geom_line()
p_osec

```

### Conclusion

We demonstrated how to obtain the expression profile of genes associated to SE genomic regions. The same procedures can be now applied to FTSEC and HGSOC specific SE by changuing Part 2 and running again Part 3 and 4.

