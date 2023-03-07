You can use any enrichment analysis method available in the Bioconductor universe as `GeDi` only has three small requirements at the input data which are independent of the enrichment method used: 

* The data has a column called "Genesets" which contain some sort of identifiers for the individual genesets.
* The data has a column called "Genes", which contains a list of genes belonging to the individual genesets.
* The genes in the "Genes" column are gene names and no other identifiers.

Otherwise, there are no specific requirements at the data. The data can be provided either as .RDS, .txt, or .xlsx file, where the .RDS file best contains a `data.frame` and the .txt file is either comma, semicolon or tab separated.

However if you want to get started with the app without having to generate data first, you can also explore `GeDi` using the demo data provided.


# `GeDi` 101: What to do when?

`GeDi` is an application with many different possibilities to interact and explore your data, which can easily become overwhelming. Hence we want to quickly explain the different elements of the app.

Whenever you see a button like this
