## What is `GeDi`?

`GeDi` is a Bioconductor package containing a Shiny application designed for analyzing functional annotation analysis results. It offers an interactive interface to explore and gain insights into datasets, providing various graphs and tables to effectively mine the data. Its main focus is on GEneset DIstances (therefore the name GeDi), using distance measurements and clustering approaches to identify patterns and themes within enriched genesets. This reduction in data size simplifies interpretation and hypothesis generation for future analyses.

## What do I need to use `GeDi`?

In order to use `GeDi` you need the results of an functional annotation analysis like the ones calculated with the Bioconductor packages [`topGO`](https://bioconductor.org/packages/release/bioc/html/topGO.html) or [`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html).

The results can be provided to GeDi in various different formats such as a data.frame saved in an RDS file, a text file or in a comma-separated file. One important requirement is that your data contains at least two columns: 

* The data has a column called "Genesets" which contains some sort of identifiers for the individual genesets. In this application, we use the term "Genesets" to refer to collections of individual genes, which share common biological characteristics or functions. Such genesets can for example be obtained from databases such as the Gene Ontology (GO), the Kyoto Encyclopedia of Genes and Genomes (KEGG), Reactome, or the Molecular Signatures Database (MSigDB). The identifiers used in these databases can be directly used as geneset identifiers in [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi). 
* The data has a column called "Genes", which contains a list of genes belonging to the individual genesets in the "Genesets" column. In order to leverage all of the functionality available in [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi), the column has to contain gene names and no other commonly used identifiers. 

For a more detailed view of an example data set, please navigate to the `Data Input` panel and click the `Load the demo data` button followed by the inspection of the `Geneset preview`. This will show you a table of an example data set, which contains the necessary columns as well as additional information. You can also use this example data to explore and get familiar with the app. 
