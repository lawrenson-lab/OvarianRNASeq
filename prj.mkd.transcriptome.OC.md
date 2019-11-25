---
title: "Machine learning analysis"
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
#render("prj.mkd.transcriptome.OC.Rmd")

```

## Introduction

In the following workflow, we walk the reader through training a OC FTSEC scores signature using FTSEC samples as training set and applying it to score TCGA tumor samples (testing set).

### Data loading

The the data is organized in a data frame with 387 samples (190 normals and 197 TCGA) as columns and 21071 exprerssed genes as rows, named as `combat_alln_mean`. Sample annotation is organized in another data frame `NormalSamplesTable`.

```{r load}
load("RData/Normal-197-TestSet-TCGA-combatGCnorm_all_data.rda")
```

```{r echo=FALSE}
library(knitr)
kable(combat_alln_mean[1:5,1:5])
```

```{r }
load("RData/NormalsSamplesTableV3.RData")
head(NormalSamplesTable)
```
Add TCGA Sample annotation in Normals sample table annotation. After the data is splited into 3 groups (NFTE, NOSE, HGSOC).

```{r }
nTCGAsamples = 197
TCGAsamplesTable = data.frame(Samples = colnames(combat_alln_mean)[191:387], CellType = c(rep("HGSOC", nTCGAsamples)), Colors=c(rep("green", nTCGAsamples)), SP="None", Batch="None" )

SamplesTable = rbind(NormalSamplesTable, TCGAsamplesTable)

all_expr = combat_alln_mean
normals_expr = all_expr[,c(1:190)]

FTSEC_expr = normals_expr[,which(SamplesTable$CellType == "NFTE")]
dim(FTSEC_expr)

OSEC_expr = normals_expr[,which(SamplesTable$CellType == "NOSE")]
dim(OSEC_expr)

TCGA_expr = all_expr[, which(SamplesTable$CellType == "HGSOC")]
dim(TCGA_expr)

```
### Mean-center the data

Find the mean center by subtracting the mean of each gene from the entire data. The mean of each gene just be in a numeric vector the same size as the number of probes, in this case 21071.

```{r }
m1 <- apply(normals_expr, 1, mean )
normals_expr <- normals_expr - m1
```
### Train model
Now we can begin to train the the one-class model with the `gelnet` function. The gelnet function can be used for Linear Regression, Binary Classification and One class Problems by using an iterative method called coordinated descent (Sokolov et al. 2016). It has four main arguments described below:

X: n by p matrix => transpose(X.r)
y: NULL for one class models
l1: coefficient for the L1-norm penalty => 0
l2: coefficient for the L2-norm penalty => 1

```{r }
X1.tr <- FTSEC_expr
X1.bk <- OSEC_expr

library( gelnet )
library( dplyr )

set.seed(13)

## Train a one-class model and store the signature in mm1
mm1 <- gelnet( t(X1.tr), NULL, 0, 1 )

```
## Leave one out Cross-validation

To test how the model’s performance by using leave one out cross-validation. This process has three steps:

<p>Train model on non-left-out data</p>
<p>Score the left-out sample against the background</p>
<p>AUC = P( left-out sample is scored above the background )</p>

```{r eval=FALSE, echo=T}
## Perform leave-one-out cross-validation
auc <- c()
for( i in 1:ncol(X1.tr) )
{
  ## Train a model on non-left-out data
  cat( "Current: ", i, " of ", ncol(X1.tr), "\n")
  one <- X1.tr[,-i]
  one_res <- gelnet( t(one), NULL, 0, 1 )

  ## Score the left-out sample against the background
  score1.bk <- apply( X1.bk[rownames(X1.tr),], 2, function(z) {cor( one_res$w, z, method="sp" )} )
  score1 <- cor( one_res$w, X1.tr[,i], method="sp" )

  ## AUC = P( left-out sample is scored above the background )
  auc[i] <- sum( score1 > score1.bk ) / length(score1.bk)
  cat( "Current AUC: ", auc[i], "\n" )
  cat( "Average AUC: ", mean(auc), "\n" )
}
```

```{r echo=F, results='hide'}
## Perform leave-one-out cross-validation
auc <- c()
for( i in 1:3 )
{
  ## Train a model on non-left-out data
  cat( "Current: ", i, " of ", ncol(X1.tr), "\n")
  one <- X1.tr[,-i]
  one_res <- gelnet( t(one), NULL, 0, 1 )

  ## Score the left-out sample against the background
  score1.bk <- apply( X1.bk[rownames(X1.tr),], 2, function(z) {cor( one_res$w, z, method="sp" )} )
  score1 <- cor( one_res$w, X1.tr[,i], method="sp" )

  ## AUC = P( left-out sample is scored above the background )
  auc[i] <- sum( score1 > score1.bk ) / length(score1.bk)
  cat( "Current AUC: ", auc[i], "\n" )
  cat( "Average AUC: ", mean(auc), "\n" )
}
```


Now we can obtain the auc mean of all iteractions.

```{r }
head(mean(auc))
```
### Predict

In this step we will use the signature that was created in the previous steps and stored as `mm1` to score the TCGA data. 

```{r }
w1 = mm1$w

X1 <- TCGA_expr[names(w1),] #54 9318 
X1 <- as.matrix(X1)

## Score via Spearman correlation
s1 <- apply( X1, 2, function(z) {cor( z, w1, method="sp", use="complete.obs" )} )
```

Scale the scores into a ratio from 0 to 1 and store as data frame.
```{r }
## Scale the scores to be between 0 and 1
s1 <- s1 - min(s1)
s1 <- s1 / max(s1)
s1 <- as.data.frame(t(s1))
```

```{r echo=FALSE}
library(knitr)
kable(s1[1:5])
```


### Conclusion
We demonstrated how to derive a gene signature from FTSEC cell type samples. The robustness of the signature was estimated through leave-one-out cross-validation. The same procedures can now be applied to OSEC samples following the same steps by changuing the test and background sets.



