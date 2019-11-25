---
title: "Cistrome analysis"
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
```

## Introduction

In the following workflow, we describe the main procedures to analyse the enrichment of differentially expressed genes (DEGs) close to experimentaly peaks, derived from Cistrome, for some interesting Trascription Factor (TF).

```{r }
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
```
### Auxiliary functions
The `buildTable` function is used to select all genes close to 50000 bp of a given peak coordinates from all Cistrome peaks file. We used the `distance` function from `Genomic Ranges` package.
```{r }
buildTable <- function(pTF, pInfo.GR, pDados.RNAseq.GR) {
  
  max_len = 50000
  table.all.genes = NULL
  
  for (i in 1: length(info.GR)) {
    
    dist = distance(pInfo.GR[i], pDados.RNAseq.GR)
    info.gene = NULL
    
    list = which(dist > 0 & dist <= max_len)
    
    if (length(list) > 0 ) {
      info = data.frame(tf = pTF, peak=i, distance=dist[list], bin=max_len, gene = pDados.RNAseq.GR$symbol[list])
      
      info.gene = rbind(info.gene, info)
    }
    
    table.all.genes = rbind(table.all.genes, info.gene)
  }
  return(table.all.genes)
}
```

The `getRandomGR` function is used to generate random peaks following the same chromosome proportion and size of the original cistrome peaks files.

```{r }
getRandomGR <- function(gr){
  seqnames <- factor(x = paste0('chr',c(1:22,'X', 'Y')), levels = paste0('chr',c(1:22,'X', 'Y')))
  hg19.len <- read.table(file = 'files/hg19.len')
  hg19.len <- hg19.len[hg19.len$V1 %in% seqnames,]
  hg19.len$V1 <- factor(x = hg19.len$V1, levels = seqnames)
  hg19.len <- hg19.len[order(hg19.len$V1),]
  
  .seqnames <- factor(as.character(seqnames(gr)), levels = levels(seqnames))
  gr.random <- foreach(i=1:length(seqnames), .combine = c) %do% {
    .n <- sum(.seqnames == seqnames[i])
    .max <- hg19.len$V2[i]
    .start <- round(x = runif(n = .n, min = 1, max = .max), digits = 0)
    .width <- width(gr[.seqnames == seqnames[i]])
    .gr <- GRanges(seqnames = seqnames[i], ranges = IRanges(start = .start, width = .width))
    return(.gr)
  }
  return(gr.random)
}

```

### Part 1
### Data loading

The differentially expressed genes are organized in a data frame. All genes are loaded and only genes up and down regulated are selected in the analysis.

```{r load}
load("RData/scenario-NOSExHGSOC-394-TCGA-all-genes.RData")
table(NOSExHGSOC$STATUS)
dados.RNAseq <- NOSExHGSOC[which(NOSExHGSOC$STATUS %in% c("Down","Up")),c(1:584)]
scenario = "OxH"

```

```{r echo=FALSE}
library(knitr)
kable(dados.RNAseq[1:5,1:10])
```


### Gene genomic annotation
We selected the genomic annotation, regarding chromosome name, start, end coordinates for human hg19 build. After the annation is merged with our differentially expressed genes.

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

Convert gene annotation into a Genomic Ranges object `dados.RNAseq.GR`. 

```{r }
dados.RNAseq.GR <- makeGRangesFromDataFrame(dados.RNAseq[,1:7],TRUE)
```
### Part 2
Define the number of interaction for random peaks.

```{r }
ninter = 10
all.rd.c.T = NULL
```
### Part 3
In this step we can select a cistrome peaks file for a given TF. We filter out non-canonical chromosomes. The information are converted to `Genomic Ranges` object.
```{r }
TF = "ASCL2"

cistrome.info <- read.delim("cistrome/all.ASCL2.peaks.bed", comment.char = "#", header = F)
colnames(cistrome.info) <- c("chr", "start", "end", "name", "score")
cistrome.info <- cistrome.info[which(cistrome.info$chr %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")),]
cistrome.info <- droplevels(cistrome.info)
cistrome.info$width <- cistrome.info$end - cistrome.info$start

info.GR <- makeGRangesFromDataFrame(cistrome.info, ignore.strand=TRUE)
```

The steps bellow call both auxiliary function.
```{r }
table.TF = buildTable(TF, info.GR, dados.RNAseq.GR)
if (is.null(table.TF)) {
  table.TF = data.frame(TF=TF, Cistrome=0)
}

for (inter in 1: ninter) {
  #print(paste0("Interaction: ", inter))
  rd.peaks = getRandomGR(info.GR)
  temp_rd = buildTable(TF, rd.peaks, dados.RNAseq.GR)
  
  if (!is.null(temp_rd)) {
    all.rd.c.T = rbind(all.rd.c.T, data.frame(Var1=TF, Freq=dim(temp_rd)[1]))
  }else {
    all.rd.c.T = rbind(all.rd.c.T, data.frame(Var1=TF, Freq=0))
  }
}

all.rd.c.T$Source = "Random"
head(all.rd.c.T)
head(table.TF)

counts.TF <- as.data.frame(table(table.TF$tf))
counts.TF$Source = "Cistrome"

table.all <- rbind(counts.TF, all.rd.c.T)
colnames(table.all) <- c("TF", "Freq", "Source")

head(table.all)

fill <- "#4271AE"
line <- "#1F3552"

p <- ggplot(table.all, aes(x=TF, y=Freq)) + 
  ggtitle(paste0("Gene list - OSEC vs. HGSOC - 243 DOWN 176 UP"))+xlab("Transcription Factors")+ylab("Number of matched peaks") +
  geom_boxplot(outlier.shape = NA, fill = fill, colour = line, alpha = 0.7) +
  scale_color_manual(values=c("#4271AE", "orange")) +
  geom_jitter(shape=16, alpha=0.7, aes(colour=Source)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
p

```

### Conslusion
We demonstrated how to analyse the enrichment of DEG close to TF using Cistrome data. The same procedures can be now applied to other DEG scenarios or different TFs by changuing Part 1 and Part 3. To increase the number of random interaction change Part 2.

