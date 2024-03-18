[`GeDi`](https://github.com/AnnekathrinSilvia/GeDi) supports many of the in the Bioconductor ecosystem availbale enrichment analysis methods. This can be achieved because [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi) only has a small set of requirements at the input data. These are independent of the used enrichment method and consist of the following:

* The data has a column called "Genesets" which contains some sort of identifiers for the individual genesets. In this application, we use the term "Genesets" to refer to collections of individual genes, which share common biological characteristics or functions. Such genesets can for example be obtained from databases such as the Gene Ontology (GO), the Kyoto Encyclopedia of Genes and Genomes (KEGG), Reactome, or the Molecular Signatures Database (MSigDB). The identifiers used in these databases can be directly used as geneset identifiers in [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi). 
* The data has a column called "Genes", which contains a list of genes belonging to the individual genesets in the "Genesets" column. In order to leverage all of the functionality available in [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi), the column has to contain gene names and no other commonly used identifiers. 

Otherwise, there are no specific requirements at the data. The data can be provided either as .RDS, .txt, or .xlsx file, where the .RDS file best contains a `data.frame` and the .txt file is either comma, semicolon or tab separated.

However if you want to get started with the app without having to generate data first, you can also explore [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi) using the demo data provided.

Although we do not set any further requirements for the input data, we highly recommend you to use differentially expressed genes or any other set of thoughfully selected genes to preselect the number and set of genes that is used with the app. This is not a requirement of [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi) as the app is designed to work with large lists of genes (althought you should be careful with the size of the input data as runtime increases with the amount of genes found in the input). However, [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi) is best suited to identify underlying themes and patterns in structured input data and hence will leverage the interpretation and evaluation of preprocessed data best. 


# `GeDi` 101: What to do when?

[`GeDi`](https://github.com/AnnekathrinSilvia/GeDi) is an application with many different possibilities to interact and explore your data. This can easily become overwhelming especially if you are new to Shiny applications. Hence we want to quickly explain the different, interactive elements of the app.

Whenever you see a button like this
