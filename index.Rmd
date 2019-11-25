---
title: "Home"
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

### Introduction
This document contains workflows, explaining how to perform all analysis described in `A Dualistic Model for High-Grade Serous Ovarian Cancer Origins`.

### Document authors
Marcos A. S. Fonseca & Felipe Segato.

### Install packages
To install the required packages to run the code below please execute the follwing code.
```{r eval=FALSE}
deps <- c("gelnet","dplyr","gdata","DT", "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db", "GenomicRanges")
for(pkg in deps)  if (!pkg %in% installed.packages()) install.packages(pkg, dependencies = TRUE)

```
### Available files

- [Superenhacers BED peaks](se).
- [Pre-processed RData](RData/).
- [Raw feature counts text file](190Normals_394TCGA_raw_fc.txt)
- [Cistrome BED peaks](cistrome/)
