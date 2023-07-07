## What is `GeDi`?

`GeDi` is a Bioconductor package containing a Shiny application designed for analyzing functional annotation analysis results. It offers an interactive interface to explore and gain insights into datasets, providing various graphs and tables to effectively mine the data. Its main focus is on GEneset DIstances (therefore the name GeDi), using distance measurements and clustering approaches to identify patterns and themes within enriched genesets. This reduction in data size simplifies interpretation and hypothesis generation for future analyses.

## What do I need to use `GeDi`?

In order to use `GeDi` you need the results of an functional annotation analysis like the ones calculated with the Bioconductor packages [`topGO`](https://bioconductor.org/packages/release/bioc/html/topGO.html) or [`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html).

The results can be provided to GeDi in various different formats such as a data.frame saved in an RDS file, a text file or in a comma-separated file. One important requirement is that your data contains at least two columns: 

* Genesets: The column Genesets should contain any kind of (ideally unique) identifiers for the individual sets.

* Genes: The column Genes should contain for each geneset a list of genes which are associated with this geneset.

For a more detailed view of an example data set, please navigate to the `Data Input` panel and click the `Load the demo data` button followed by the inspection of the `Geneset preview`. This will show you a table of an example data set, which contains the necessary columns as well as additional information. You can also use this example data to explore and get familiar with the app. 
