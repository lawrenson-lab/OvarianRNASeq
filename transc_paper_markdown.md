---
title: "Normalization and DEG analysis"
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

#Introduction 

In the following workflow, we walk the reader through downloading the normal ovarian(OSEC,FTSEC) data and TCGA HGSOC level 1 data, and after aligning to a genome reference, get read count using *featureCount()*, correct for GC and batch effect with *EDASeq* ,*sva*,  (... to be continued)

#Read Count

To get read counts per gene and build our working matrix we need to count mapped reads for genomic features, in this case, genes. We are using the function featureCounts() from Rsubread and an external GTF file containing the GENCODE v19 annotation, even though having in-built annotation in the function.


```{r, eval = F}

library(Rsubread)
bam = dir(pattern = "bam")
fc = featureCounts(bam, annot.ext = "gencode.v19.annotation.gtf", isGTFAnnotationFile= TRUE ,nthreads = 10)
```

This should take a while to run, especially for the HGSOC TCGA data that has more sequencing depth. ( -check- )
In this study we divided the *featureCounts* into batches, for saving the process while running and to not lose the workspace for memory issues.

After assembling your data matrix with all samples (HGSOC, FTSEC and OSEC) your data should look like this:

```{r, eval=T}
load("RData/merge-fc-new-workspace.RData")

allAnno[1:5, 1:5]

```

After getting read count and building a matrix we can start to look in how the samples are distributed.

Calculate principal components for each sample and plot a PCA-plot:


```{r}
library(ggplot2)
library(cowplot)

pca.raw = prcomp(t(allAnno))


ggplot(data= data.frame(pca.raw$x[,1:2]),aes(x = PC1, y = PC2, color= factor(SamplesTableNew$Colors), shape = as.factor(SamplesTableNew$Batch))) +
  ggtitle("FTSEC - OSEC - HGSOC - raw counts (192 Normals & 394 TCGA samples)") +
  geom_point(size = 3) +
  scale_colour_manual(name = "Batch & CellType", labels = c("OSEC","HGSOC", "FTSEC"), values = c("#377EB8","darkgreen", "red")) +
  scale_shape_manual(name = "Batch", labels = c("1", "2", "HGSOC set"), values = c(17, 19, 21)) +
  
  xlab(paste0("PC1 ", prettyNum(summary(pca.raw)$importance[2,1]*100,
                                digits = 2, decimal.mark = "."), "%")) +
  ylab(paste0("PC2 ", prettyNum(summary(pca.raw)$importance[2,2]*100,
                                digits = 2, decimal.mark = "."), "%"))

```

this is a preliminar non-supervised analysis of the data, which needs to be normalized.

# Normalization

In the following section we will be describing how to normalize the data for GC bias content and to correct for batch effect.

GC content bias can describe the dependence between read coverage and GC content found in sequencing data.

Batch effects are technical sources of variation that have been added to the samples during handling. It is important that this type of technical variation does not confound with the biology, meaning that biological treatment groups overlap with technical groups. Technical batch effects confounding with the biology can be avoided by careful planning of every step of the experiment. This is an important part of the experimental design. 

```{r}
library(EDASeq)
library(sva)

head(SamplesTableNew)

data_all = newSeqExpressionSet(counts = as.matrix(allAnno[common , ]),
                               featureData= output[common, ],
                               phenoData= data.frame(conditions=   SamplesTableNew$CellType,                                                      row.names= colnames(allAnno)))

biasPlot(data_all, "GC", log = T, ylim = c(-5,5), main = "Bias plot")


```
Within-lane normalizatiion procedures reduce GC-content bias and lead to more accurate estimates of expression FC and tests of differential expression.
```{r}
dataW = withinLaneNormalization(data_all, "GC", which = "full")

```
Bias plot after withinLaneNormalization:
```{r}
biasPlot(dataW, "GC", log = T, ylim = c(-5,5), main = "Bias plot - after GC normalization")
all_gcnorm = counts(dataW)


```

PCA plot after GC-bias normalization:

```{r}
pca_gc = prcomp(t(all_gcnorm))

ggplot(data= data.frame(pca_gc$x[,1:2]),aes(x = PC1, y = PC2, color= as.factor(SamplesTableNew$Colors), shape = as.factor(SamplesTableNew$Batch))) +
  ggtitle("FTSEC - OSEC - HGSOC - GC Normalization") +
  geom_point(size = 3) +
  scale_colour_manual(name = "Batch & CellType", labels = c("OSEC","HGSOC", "FTSEC"), values = c("#377EB8","darkgreen", "red")) +
  scale_shape_manual(name = "Batch", labels = c("1", "2", "TCGA set"), values = c(17, 19, 21)) +
  
  xlab(paste0("PC1 ", prettyNum(summary(pca_gc)$importance[2,1]*100,
                                digits = 2, decimal.mark = "."), "%")) +
  ylab(paste0("PC2 ", prettyNum(summary(pca_gc)$importance[2,2]*100,
                                digits = 2, decimal.mark = "."), "%"))
```

Next step it's to correct for batch effect:

```{r}
library(sva)
library(edgeR)

log_all_gcnorm = cpm(all_gcnorm, log =T, prior.count = 2)
mod = model.matrix(~as.factor(SamplesTableNew$CellType))
mod = mod[, 1:2]
combat_alln_mean = ComBat(dat= log_all_gcnorm,batch = SamplesTableNew$Batch, mod = mod ,  par.prior = T, mean.only = T)

```

See the data PCA after GC-bias and batch effect correction:

```{r}
pca_combat = prcomp(t(combat_alln_mean))

ggplot(data= data.frame(pca_combat$x[,1:2]),aes(x = PC1, y = PC2, color= as.factor(SamplesTableNew$Colors), shape = as.factor(SamplesTableNew$Batch))) +
  ggtitle("FTSEC - OSEC - HGSOC - GC and Combat Normalization") +
  geom_point(size = 3) +
  scale_colour_manual(name = "Batch & CellType", labels = c("OSEC","HGSOC", "FTSEC"), values = c("#377EB8","darkgreen", "red")) +
  scale_shape_manual(name = "Batch", labels = c("1", "2", "TCGA set"), values = c(17, 19, 21)) +
  
  xlab(paste0("PC1 ", prettyNum(summary(pca_combat)$importance[2,1]*100,
                                digits = 2, decimal.mark = "."), "%")) +
  ylab(paste0("PC2 ", prettyNum(summary(pca_combat)$importance[2,2]*100,
                                digits = 2, decimal.mark = "."), "%"))
```

Now that the data is normalized we can see the separation of the groups, highlighting batches with shapes. The data is now ready for further analysis.




