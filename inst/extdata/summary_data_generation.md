*This information is also contained in the `GeDi` package vignette.*

# How to prepare the data: A quick summary

In order to use `GeDi` you will need results from a functional annotation analysis. In this vignette, we will show how to perform an enrichment analysis on differentially expressed (DE) genes of the `r BiocStyle::Biocpkg("macrophage")` dataset. 
In the first step, we will load the `macrophage` data and generate a `DESeqDataset` as the differential expression analysis will be performed with `r BiocStyle::Biocpkg("DESeq2")`.

```r
library("macrophage")
library("DESeq2")

data("gse", package = "macrophage")

dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
# changing the ids to Ensembl instead of the Gencode used in the object
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
dds_macrophage
```

Now that we have our `DESeqDataset`, we can perform the DE analysis. In this vignette, we will use the results of the comparison of two different conditions of the data set, namely `IFNg` and `naive`, while controlling for the cell line of origin (which has 6 levels, namely `r knitr::combine_words(levels(dds_macrophage$condition))`).
Before we perform the DE analysis, we filter lowly expressed features from the data set. In this example we filter all genes that do not have at least 10 counts in at least 6 samples (where 6 is the size of the smallest group in the data).
Afterwards, we pefrom the DE analysis and test against a null hypothesis of a log2FoldChange of 1 in order to ensure that we call those genes with a consistent *and* robust expression change DE.
In a last step, we add the gene symbols to the resulting `DataFrame` which will later serve as our Genes column in the input data to `GeDi`.

```r
keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
dds_macrophage <- dds_macrophage[keep, ]
dds_macrophage

dds_macrophage <- DESeq(dds_macrophage)

res_macrophage_IFNg_vs_naive <- results(dds_macrophage,
  contrast = c("condition", "IFNg", "naive"),
  lfcThreshold = 1, alpha = 0.05
)
res_macrophage_IFNg_vs_naive$SYMBOL <- rowData(dds_macrophage)$SYMBOL
```

After we performed the differential expression analysis, we can now perform the functional annotation analysis. For this purpose, we are first going to extract the DE genes from the previously generated results as well as determine the background genes to be used for the functional enrichment. 

For the enrichment analysis, we are going to use the overrepresentation anlysis method implemented in the `r BiocStyle::Biocpkg("topGO")` package. In order to facilitate the later usage of these results in `GeDi`, we use the `topGOtable` wrapper function available in the `r BiocStyle::Biocpkg("pcaExplorer")`. This function uses per default the `BP` ontology and the `elim` method to decorrelate the GO graph structure and deliver less redundant functional categories and generated a `DataFrame` object that can readily be used in `GeDi`.

As previously mentioned, also enrichment results generated with `r BiocStyle::Biocpkg("clusterProfiler")` can be used. Especially results generated with the `enrichGO` method have been tested during the development of `GeDi`, but also the results form the `enrichKEGG` and `enrichPathway` method can be used as input.

```r
library("pcaExplorer")
library("GeneTonic")
library("AnnotationDbi")
# we extract the differential expression genes from the result object
de_symbols_IFNg_vs_naive <- deseqresult2df(res_macrophage_IFNg_vs_naive, FDR = 0.05)$SYMBOL
# we determine the background genes for the ORA
bg_ids <- rowData(dds_macrophage)$SYMBOL[rowSums(counts(dds_macrophage)) > 0]

library("topGO")
topgoDE_macrophage_IFNg_vs_naive <-
  pcaExplorer::topGOtable(de_symbols_IFNg_vs_naive,
    bg_ids,
    ontology = "BP",
    mapping = "org.Hs.eg.db",
    geneID = "symbol",
    topTablerows = 500
  )

````
Now, we have a enrichment analysis result. However, before we can use this result in `GeDi` we have to transform the data to fit the requirements of `GeDi`. The package has three requirements at the input data: 

* The data has a column called "Genesets" which contain some sort of identifiers for the individual genesets.
* The data has a column called "Genes", which contains a list of genes belonging to the individual genesets.
* The genes in the "Genes" column are gene names and no other identifiers.

The last requirement is already fulfilled with the prepared enrichment result, however, we still have to rename the columns. In out case the respective columns are the "GO.ID" which will become "Genesets" and the "genes" column which will become "Genes".

```r
names(topgoDE_macrophage_IFNg_vs_naive)[names(topgoDE_macrophage_IFNg_vs_naive) == "GO.ID"] <- "Genesets"
names(topgoDE_macrophage_IFNg_vs_naive)[names(topgoDE_macrophage_IFNg_vs_naive) == "genes"] <- "Genes"
```

Now the data is prepared for the use in `GeDi`.
