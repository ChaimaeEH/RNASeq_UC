# RNASeq_UC

The objective of this project is to explore a publicly available transcriptomics dataset (bulk RNASeq). Biopsies were collected from multiple regions of the colon from 12 ulcerative colitis (UC) patients. 

Four biopsies were collected from each patient, from inflamed and non-inflamed areas of the colon.
The entirety of the data (count matrix and associated metadata) can be downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107593. The Series Matrix TXT file should be used to extract the study metadata variables.

**READ THE HTML FILE**

To read the html file: 
- click on the html file then click on "View Raw"
- copy the link of the page 
- go to https://htmlpreview.github.io/ and paste the link and click on "Preview".

**Questions :** 

1)	In a PCA using the top 500 most variable genes, what metadata variable is best associated with the first principal component? 
2)	What are the top 10 genes contributing to principal component 1 & 2 (consider positive and negative loadings)? 
3)	Perform differential gene expression analysis between inflamed and non-inflamed biopsies - what are the top 20 genes increased or decreased in inflamed biopsies?
4)	Do these top 20 genes (found above, increased or decreased in inflamed biopsies) display a concordant profile across the patients?
5)	Using the Gene Set Enrichment Analysis method, find the significant Reactome Pathways enriched in inflamed or non-inflamed biopsies.
6)	Which tissues or cell types are enriched among the top (for ex with log2 Fold Change <1.5 or >1.5) significantly up- or down-regulated genes in inflamed biopsies? 


