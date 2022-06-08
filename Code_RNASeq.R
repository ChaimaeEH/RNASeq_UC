---
title: "RNASeq_UC"
author: "Chaimae"
date: "04/06/2022"
output:
  html_document: default
  pdf_document: default
editor_options: null
chunk_output_type: inline
---


# Systems Biology 2022 data analysis assignment

## Installation and loading of libraries

```{r}
#install.packages("dplyr") # to manipulate gene expression data 
#install.packages("tidyverse")
#install.packages("FactoMineR")
#install.packages("factoextra")
#install.packages("ggplot2") # to visualise and make graphs
#install.packages("pheatmap")


#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GEOquery")
#BiocManager::install("DESeq2")
#BiocManager::install("EnhancedVolcano")
#BiocManager::install('PCAtools')

#For GSEA 
#BiocManager::install("msigdbr")
#BiocManager::install("clusterProfiler")
#BiocManager::install("fgsea", force = TRUE)
#BiocManager::install("org.Hs.eg.db")

```

```{r, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
library("dplyr")
library("tidyverse")
library("GEOquery")

library("FactoMineR")
library("factoextra")
library("ggplot2")

library("DESeq2")
library("EnhancedVolcano")
library("pheatmap")

library("PCAtools")
library("airway")

#For GSEA 
library(tidyverse)         #Data manipulation
library(msigdbr)           #Broad gene set database
library(clusterProfiler)   #Hypergeo enrichment
library(fgsea)             #GSEA
library(enrichplot)
```

## Download Data 

The objective of this assignment is to explore a publicly available transcriptomics dataset (bulk RNASeq). Biopsies were collected from multiple regions of the colon from 12 ulcerative colitis (UC) patients. Four biopsies were collected from each patient, from inflamed and non-inflamed areas of the colon.

### Gene expression data 

```{r}
data <- read.delim("GSE107593_raw_reads_BCHRNAseq.txt")
head(data)
```

```{r}
nrow = nrow(data)
ncol = ncol(data)

cat("Number of rows :", nrow, "\n")

cat("Number of columns :", ncol)
```

```{r}
length(unique(data$gene_name))
```

The number of gene_names is far lower that the number of Ensembl Id. 

```{r}
glimpse(data)
```

Converting the variables from chr to factor 

```{r}
data$Row <- factor (data$Row) 
data$strand <- factor (data$strand) 
data$Type <- factor (data$Type) 
data$gene_type <- factor (data$gene_type) 
data$chr <- factor (data$chr) 
```

And see the different levels of the categorical variables. 

```{r}
levels(data$strand)

levels(data$Type)

levels(data$gene_type)

levels(data$chr)

```

Analysis : 
 - 9 descriptive variables : Row, Type, gene_name, gene_type, chr, start, end, strand and length 
 - 48 samples : Four biopsies were collected from each of the 12 ulcerate colitis (UC) patients (4*12 = 48)

 
**Preparation and selection of the columns we need** 

The rows correspond to the genes Ensemble Id and the columns to the samples 

```{r}
counts_data <- data[10:length(data)] # select only the 48 biopsies 
rownames(counts_data) <- data$Row # rename the row lines, we can't use the column gene_name because it is not composed  of unique values. 

head(counts_data)

```

### Metadata 

```{r}
gse <- getGEO(GEO = 'GSE107593', GSEMatrix = TRUE )
```

Look at the object : 

```{r}
gse
```

Get the pheno data of the first element of that list 

```{r}
metadata <- pData(phenoData(gse[[1]]))
dim(metadata)
```
 
The metadata has 48 rows and 45 variables. 

```{r}
head(metadata, 3)
```

In the metadata we have information about each of the sample, their status (inflamed or not inflamed), how they are processed, the description, ... For our analysis, we do not require all the 45 columns but only some columns of interest. 

Therefore we are going to select only some of the columns and modify the column source_name_ch1 in order to fit with the name of the columns in the gene expression dataframe. 

```{r}
# select the columns you want and modify the samples column to fit the one of the gene expression
metadata.modify <- metadata %>%
                      dplyr::select('source_name_ch1','location:ch1','site:ch1','status:ch1', 'subject:ch1') %>%
                      mutate(source_name_ch1 = gsub("Colon_", "X",source_name_ch1)) %>%
                      mutate(source_name_ch1 = gsub(" ", ".",source_name_ch1))  %>%  # because of 3 samples
                      mutate(source_name_ch1 = gsub("+", ".",source_name_ch1, fixed = TRUE))  %>%
                      mutate(source_name_ch1 = gsub("-", ".",source_name_ch1, fixed = TRUE))

# modify the 3 values for the samples that had the same number positive and negative
metadata.modify$source_name_ch1[4] <- paste0(metadata.modify$source_name_ch1[4], '.1')     # for sample 157+3
metadata.modify$source_name_ch1[15] <- paste0(metadata.modify$source_name_ch1[15], '.1')   # for sample 1077+1
metadata.modify$source_name_ch1[23] <- paste0(metadata.modify$source_name_ch1[23], '.1')   # for sample 1214+4

# Renaming columns
colnames(metadata.modify) <- c('samples','location', 'site', 'status', 'subject')

#show the first lines of the metadata 
head(metadata.modify)

```

Then we should order the 48 lines by a certain value in order it to be the same as the order of the columns in the gene expression data. 

```{r}
# same ordering between the columns of gene expression data 
metadata.modify <- metadata.modify[match(colnames(counts_data), metadata.modify$samples),]
# rename the rows 
rownames(metadata.modify) <- metadata.modify$samples
# remove the samples column 
metadata.modify <- metadata.modify[,-c(1)]
# show the data 
head(metadata.modify)
```

Check that the ordering is right and that all the column values match the row names of the Metadata

```{r}
all(colnames(counts_data) %in% rownames(metadata.modify))
all(colnames(counts_data) == rownames(metadata.modify))
```




## Introduction with DESeq2

We are going to use the package DESeq2 for the questions of this assignment. But before starting to really answer the questions we will introduce this package and make a differential expression analysis with our data. And especially calculate the DESeq2 objects that we will need. 

The first step is to construct a DESeqDataSet object 

```{r message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
dds <- DESeqDataSetFromMatrix(countData = counts_data,  # the data with 48 biopies and EnsemblId as rownames
                       colData = metadata.modify,       # the Metadata with four variables 
                       design = ~ status)               # Inflamed or Not Inflamed 
# Let's take a look at that 
dds
```

Let's perform some pre-processing to help to reduce the size of the dds object as well as increase the speed of computation : removing the rows with low gene counts and keeping the rows that have at least 10 reads total. 

```{r}
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep, ]
dds
```

We can see that we went from 60498 genes to 30952 genes. We filtered out many rows that had low gene counts. 

Then we want to compare the Inflamed VS the Non-Inflamed levels. We want to set the factor level to the Non-Inflamed. (Otherwise the choice will by done alphabetically). 

```{r}
dds$status <- relevel(dds$status, ref = "Non-Inflamed")
```

Let's run the DESeq function now :

```{r message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
dds <- DESeq(dds)
```



## Question 1 : In a PCA using the top 500 most variable genes, what metadata variable is best associated with the first principal component?
 
 
Since we have 30952 genes and we only want the 500 the most variable genes, we have to remove :
1- 500*100/30952 (cross-product)

```{r}
vsdata <- vst(dds, blind=FALSE) # estimate dispersion trend and apply a variance stabilizing transformation
vst <- assay(vsdata)

proportion = 1-500/30952
p <- pca(vst, metadata =  metadata.modify, removeVar = proportion )  # PCA 
length(p$xvars)
```

Indeed, only 500 genes are left. 

```{r message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
screeplot(p, axisLabSize = 8, titleLabSize = 20)
```

We can see on this plot that 59% of the information (variances) contained in the data is retained by the first principal component.

Let's make a biplot of the individuals and of the 10 most important variables. We color by the status. 

```{r fig.height=6, fig.width=12}
biplot(p, showLoadings = TRUE,
       labSize = 2, 
       pointSize = 2, 
       sizeLoadingsNames = 3,
       colby = 'status',
       colLegendTitle =  'Status',
       legendPosition = "right",
       legendLabSize = 20,
       legendTitleSize = 20)
```

The metadata variable status which gives us an information about the inflammation of different areas of the colon (Inflamed or not Inflamed) is best associated with the first principal component (on the x-axis). We can clearly distinguish two groups : on the right, biopsies collected from inflamed colon areas and on the left biopsies collected from non-inflamed colon areas. 


## Question 2 :	What are the top 10 genes contributing to principal component 1 & 2 (consider positive and negative loadings)?


**TOP 10 genes contributing to principal component 1** 

```{r}
# We order by the absolute value of the loadings. 
top10_genes_PCA1 <- p$loadings[order(abs(p$loadings[,1]), decreasing = TRUE),]
head(top10_genes_PCA1[1],10)
```
We have 2 negative loadings and 8 positive loadings. 

```{r}
top10_genes_PCA1 <- head(rownames(p$loadings[order(abs(p$loadings[,1]), decreasing = TRUE),]), 10)
top10_genes_PCA1
```

Let's see the name and information of those 10 genes contributing the most to PCA1. 

```{r}
data[which(data$Row %in% top10_genes_PCA1),][1:9]
```
We can notice that 8/10 are proteine coding genes and 2/10 are IG_C_gene. 

**TOP 10 genes contributing to principal component 2 :**

```{r }
# We order by the absolute value of the loadings. 
top10_genes_PCA2 <- p$loadings[order(abs(p$loadings[,2]), decreasing = TRUE),]
head(top10_genes_PCA2[2],10)

```

We have 9 negative loadings and 1 positive loadings. 

```{r}
top10_genes_PCA2 <- head(rownames(p$loadings[order(abs(p$loadings[,2]), decreasing = TRUE),]), 10)
top10_genes_PCA2
```

Name and information of those 10 genes contributing the most to PCA1. 

```{r}
data[which(data$Row %in% top10_genes_PCA2),][1:9]
```

We can notice that 9/10 genes best associated with PC2 are located in the Y chromosome and the 10th is on the X chromosome. (We noticed that the only one with negative loading is the one on Chr X)
And 7/10 are protein coding genes. 


## Question 3 :	Perform differential gene expression analysis between inflamed and non-inflamed biopsies - what are the top 20 genes increased or decreased in inflamed biopsies?

Let's save and explore the results of the dds. 

```{r}
# we put a threshold of padj = 0.05
res <- results(dds, alpha = 0.05, lfcThreshold = 0) # alpha is padj 
res
```

The log2 fold change is calculated between the Inflamed vs Non.Inflamed levels.
The statistical test used is the Wald test. 
This dataframe has multiple columns. 
  - baseMean is the average of the normalized count  taken over all the samples. 
  - the log2FoldChange is the fold change of the gene in the Inflamed condition and compared to the Non.Inflamed
    The >0 values are the up- regulated genes in the Inflamed condition. 
    The <0 values are the down-regulated genes in the Inflamed condition.
  - lfcSE : Standard Estimates for the Log2Fold change 
  - stat : stat values for the wald test for each gene 
  - pvalue : pvalue of the test statistics for each gene 
  - padj : corrected pvalue for multiple testing (indeed 5% of our deferentially expressed genes are not really             deferentially expressed but only due to random chances : they are false positive. 5% of 30k is a lot )

Let's see the summary of the results : 

```{r}
summary(res)
```

8076 genes are upregulated 
6652 genes are down regulated. 

Make an **MA plot** to visualise the results.  

```{r}
plotMA(res)
```


We want to select the top 20 genes increased or decreased in inflamed biopsies
Therefore, we order the genes by the lowest padj and select the TOP 20  


```{r}
res <- res[order(res$padj),] # order by lowest padj
head(res, 10)
```

```{r}
top20_lowest_padj <- head(rownames(res), 20)
top20_lowest_padj
```

Name and information of those 20 genes with the lowest padj. 

```{r}
data[which(data$Row %in% top20_lowest_padj),][1:9]
```

Those genes are the Top20 genes increased or decreased in inflamed biopsies. 
Let's make a volcano plot to show them : 

```{r fig.height=10, fig.width=10}
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                xlim=c(-10,10),
                y = 'pvalue',
                selectLab = top20_lowest_padj,
                title = "EnhancedVolcano plot ",
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = res$padj[19], # to see only the 20 genes that have been selected 
                FCcutoff = 0,
                pointSize = 1.0,
                labSize = 3,
                labCol = 'black',
                labFace = 'bold',
                #boxedLabels = TRUE,
                colAlpha = 3/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.1,
                colConnectors = 'black')
```






## Question 4 : Do these top 20 genes (found above, increased or decreased in inflamed biopsies) display a concordant profile across the patients?

In order to see if these 20 genes display a concordant profile across the patients, we will plot for each of them the counts of reads for a single gene across the groups.

```{r}
par(mfrow=c(2,3))
for (val in 1: length(top20_lowest_padj))
{
    plotCounts(dds, gene=top20_lowest_padj[val], intgroup="status")
}
```

We can see that for those 20 genes, we can each time distinguish 2 groups between the Non-Inflamed and Inflamed patients. It is more ambiguous for gene ENSG00000169035.11. 


```{r}
df <- as.data.frame(colData(dds)[,c("status","site")])

pheatmap(vst[top20_lowest_padj,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col = df[1])
```

These top 20 genes display a concordant profile across the patients. We can distinguish e   sily distinguish the two groups on the heatmap. Some genes seem to be under-expressed in the Inflamed patients whereas others seem to be overexpressed. 


## Question 5 : Using the Gene Set Enrichment Analysis method, find the significant Reactome Pathways enriched in inflamed or non-inflamed biopsies.

Gene set enrichment analysis (GSEA) (also functional enrichment analysis) is a method to identify classes of genes or proteins that are over-represented in a large set of genes or proteins, and may have an association with disease phenotype.

The list of deferentially expressed genes is sometimes so long that its interpretation becomes cumbersome and time consuming. A common downstream procedure is gene set testing. It aims to understand which pathways or gene networks the deferentially expressed genes are implicated in.


**Get gene set data**

First, we extract the lists of genes in C2 terms from the package `msigdbr`. We format it as a data frame so we can use the `tidyverse`.

```{r}
#C2
C2 <- as.data.frame(msigdbr(species = "Homo sapiens", 
                                category = "C2",
                                subcategory = "CP:REACTOME"))
head(C2,3)
```

Here, we see the start of the C2 data including :

* gs_cat: gene set category, in this case C2 for Curated Gene sets 
* gs_name: short gene set name
* gene_symbol: HGNC symbol for each each in the gene set
* gs_description: sentence describing the gene set name further
* additional variables: other gene identifiers like ENSEMBL, ENTREZ, etc


**Extract the matrix of genes and the log2FoldChange values**

GSEA compares mean gene expression fold change between two states. We extract the genes and the log2FoldChange from the results of DESeq2 and we sort the list in decreasing order of L2FC (required for clusterProfiler). 

```{r}
# we want the log2 fold change 
gene_matrix <- res$log2FoldChange 

# name the vector by removing everything that is after the "." in the rownames of the ensemblId
names(gene_matrix) <- gsub("\\..*","", rownames(res)) 

# omit any NA values 
gene_matrix<-na.omit(gene_matrix)

# sort the list in decreasing order (required for clusterProfiler)
gene_matrix = sort(gene_matrix, decreasing = TRUE)

head(gene_matrix)
```

This is the format we need for GSEA, for the package to run 

**Run GSEA**

```{r message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
gsea_results <- GSEA(
  geneList = gene_matrix, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    C2,
    gs_name,
    ensembl_gene
  )
)
```

We can access the results from our `gsea_results` object using `@result`
Let’s convert the contents of result into a data frame that we can use for further analysis.


```{r}
gsea_result_df <- data.frame(gsea_results@result)
```
**Visualize significant GSEA**

The most common visualization for GSEA is the enrichment or normalized enrichment score. Let's plot `NES` for the subset of significant GSEA at FDR < 0.05. This is actually pretty strict for GSEA and people often go up to FDR < 0.2.

```{r fig.height=5, fig.width=6.5}
gsea_result_df %>% 
  filter(p.adjust <= 0.05) %>% # filter only those we a certain value of padj
  #Beautify descriptions by removing _ and REACTOME
  mutate(ID = gsub("REACTOME_","", ID),
         ID = gsub("_"," ", ID)) %>% 
  
  ggplot(aes(x=reorder(ID, NES), #Reorder gene sets by NES values
             y=NES)) +
    geom_col() +
    theme_classic() +
    #Force equal max min
    #lims(y=c(-3.2,3.2)) +
    #Some more customization to pretty it up
    #Flip x and y so long labels can be read
    coord_flip() +
    #fix labels
    labs(y="Normalized enrichment score (NES)",
         x="Gene set",
         title = " Reactom GSEA (padj<0.05) ") +
    theme(axis.text.y = element_text(face="bold", color="black", size=2))
```

NES > 0 (positive enrichment), corresponds to higher expression in Inflamed areas (up-regulated)
NES < 0 (negative enrichment), corresponds to lower expression in Inflamed areas (down regulated)

**Most Positive NES**

Let’s look at the 3 gene sets with the most positive NES.

```{r}
gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(NES, n = 3)
```

The gene set `REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT` has the most positive NES score.

It's **GSEA plot** (Plot of the Running Enrichment Score (green line) for a gene set as the analysis walks down the ranked gene list, including the location of the maximum enrichment score (the red line). The black lines in the Running Enrichment Score show where the members of the gene set appear in the ranked list of genes, indicating the leading edge subset.)

```{r}
most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT",
  title = "REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT",
  color.line = "#0d76ff"
)
most_positive_nes_plot
```

Notice how the genes that are in the gene set, indicated by the black bars, tend to be on the left side of the graph indicating that they have positive gene-level scores. The red dashed line indicates the enrichment score, which is the maximum deviation from zero. As mentioned earlier, an enrichment is calculated by starting with the most highly ranked genes (according to the gene-level log2 fold changes values) and increasing the score when a gene is in the pathway and decreasing the score when a gene is not in the pathway.

**Most Negative NES**

Let’s look for the 3 gene sets with the most negative NES.

```{r}
gsea_result_df %>%
  # Return the 3 rows with the smallest (most negative) NES values
  dplyr::slice_min(NES, n = 3)
```
The gene set `REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT` has the most negative NES.

```{r}
most_negative_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT",
  title = "REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT",
  color.line = "#0d76ff"
)
most_negative_nes_plot
```

This gene set shows the opposite pattern – genes in the pathway tend to be on the right side of the graph. Again, the red dashed line here indicates the maximum deviation from zero, in other words, the enrichment score. A negative enrichment score will be returned when many genes are near the bottom of the ranked list.

**Dot Plots**

```{r}
require(DOSE)
dotplot(gsea_results, showCategory=3, split=".sign") + facet_grid(.~.sign)
```


## Question 6 : Which tissues or cell types are enriched among the top (for ex with log2 Fold Change <1.5 or >1.5) significantly up- or down-regulated genes in inflamed biopsies?

**We select the Top significantly genes and their L2FC**

The Top significantly up- or down-regulated genes in inflamed biopsies (log2 Fold Change <1.5 or >1.5)

```{r}
# we want the log2 fold change that are > or < to 1.5 
gene_matrix_l2FC_1.5 <-res %>%
                            as.data.frame() %>%
                            select('log2FoldChange') %>%
                            filter(abs(log2FoldChange) > 1.5)

# extract only the column L2FC 
gene_matrix_l2FC_1.5 <- gene_matrix_l2FC_1.5$log2FoldChange

# name the vector by removing everything that is after the "." in the rownames of the ensemblId
names(gene_matrix_l2FC_1.5) <- gsub("\\..*","", rownames(gene_matrix_l2FC_1.5)) 

# omit any NA values 
gene_matrix_l2FC_1.5<-na.omit(gene_matrix)

# sort the list in decreasing order (required for clusterProfiler)
gene_matrix_l2FC_1.5 = sort(gene_matrix_l2FC_1.5, decreasing = TRUE)

head(gene_matrix_l2FC_1.5)

```

**Get gene set data**

We extract the lists of genes in C8 terms from the package `msigdbr`. We format it as a data frame so we can use the `tidyverse`.

```{r}
#C8
C8 <- as.data.frame(msigdbr(species = "Homo sapiens", category = "C8"))
head(C8)
```
This is the format we need for GSEA, for the package to run 

**Run GSEA**

```{r message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
gsea_results_C8 <- GSEA(
  geneList = gene_matrix_l2FC_1.5, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    C8,
    gs_name,
    ensembl_gene
  )
)
```

We can access the results from our `gsea_results` object using `@result`
Let’s convert the contents of result into a data frame that we can use for further analysis and write to a file later.

```{r}
gsea_result_df_C8 <- data.frame(gsea_results_C8@result)
```

**Visualize significant GSEA**

The most common visualization for GSEA is the enrichment or normalized enrichment score. Let's plot `NES` for the subset of significant GSEA at FDR < 0.05. This is actually pretty strict for GSEA and people often go up to FDR < 0.2.

```{r fig.height=5, fig.width=6.5}
gsea_result_df_C8 %>% 
  filter(p.adjust <= 0.05) %>% # filter only those we a certain value of padj
  #Beautify descriptions by removing _ and REACTOME
  mutate(ID = gsub("_"," ", ID)) %>% 
  
  ggplot(aes(x=reorder(ID, NES), #Reorder gene sets by NES values
             y=NES)) +
    geom_col() +
    theme_classic() +
    #Force equal max min
    #lims(y=c(-3.2,3.2)) +
    #Some more customization to pretty it up
    #Flip x and y so long labels can be read
    coord_flip() +
    #fix labels
    labs(y="Normalized enrichment score (NES)",
         x="Gene set",
         title = " Reactom GSEA (padj<0.05) ") +
    theme(axis.text.y = element_text(face="bold", color="black", size=2))
```

**Most Positive NES**

Let’s look at the 5 cell type signature gene sets with the most positive NES.

```{r}
gsea_result_df_C8 %>%
  # This returns the 5 rows with the largest NES values
  dplyr::slice_max(NES, n = 5)
```

We can see that indeed the first cell type signature gene sets is `GAO_LARGE_INTESTINE_ADULT_CI_MESENCHYMAL_CELLS` which is concordant since the biopsies have been collected from that area of the human body.

Amélioration : pas sûr DU TOUT de ma conclusion. surtout que le numéro 2 c'est un rein qui revient quoi. 

**Most Negative NES**

Let’s look for the 5 gene sets with the most negative NES.

```{r}
gsea_result_df_C8 %>%
  # Return the 5 rows with the smallest (most negative) NES values
  dplyr::slice_min(NES, n = 5)
```

The TOP5 gene sets with the most negative NES are related to the intestinal epithelium. 

Enterocytes are one of the four main cell types of the intestinal epithelium


# Thank You 

Hope to see you soon

Chaimae EL HOUJJAJI

## References 

* https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html
* https://www.youtube.com/watch?v=Bzu4_yDcBLY
* https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C2
* https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/


