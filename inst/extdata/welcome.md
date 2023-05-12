# Welcome to `GeDi`!

Discover the power of [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi) - an R package that offers a Shiny application designed for analyzing functional annotation analysis results. 
With [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi), you can embark on an interactive journey to explore and gain a deeper understanding of the underlying mechanisms within your dataset. This application provides a user-friendly interface while maintaining accessibility to a comprehensive selection of graphs and tables, enabling you to effectively mine your data.

Leveraging the robust infrastructure of the Bioconductor project, [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi) ensures interactive and reproducible exploration of enrichment analysis results. Exporting and sharing graphs, tables, and interactive HTML reports with collaborators becomes a seamless process. The dynamic user interface delivers a rich array of content and information, categorized into thematic tasks. Our aim is to facilitate proper exploration by bridging the gap between life scientists and experienced bioinformaticians. Through high-quality visualizations and accessible documentation, [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi) promotes effective communication between both sides.

The application simplifies and summarizes the often overwhelming amount of information presented by enrichment results. Rather than being confronted with an extensive list of enriched genesets, [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi) summarizes the list through various visualization and analysis methods. The overarching theme of [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi) revolves around Geneset-Distances. By utilizing the input of enriched terms, [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi) employs diverse distance measurements to determine the distances between individual genesets. Clustering approaches are then applied to identify larger clusters, patterns, and overall themes within the data. This reduction in data size facilitates interpretation and the generation of hypotheses, ultimately influencing future analyses. 


# How to get started

To utilize [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi), you only need one input: the results of a functional annotation analysis, such as an enrichment analysis result calculated with [`topGO`](https://bioconductor.org/packages/release/bioc/html/topGO.html) or [`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html). Ideally, [`GeDi`](https://github.com/AnnekathrinSilvia/GeDi) is incorporated as an additional step in an RNA-seq analysis. It seamlessly integrates into the repertoire of Bioconductor packages developed by our group for RNA-sequencing analysis, including [`pcaExplorer`](https://bioconductor.org/packages/release/bioc/html/pcaExplorer.html), [`ideal`](https://bioconductor.org/packages/release/bioc/html/ideal.html), and [`GeneTonic`](https://bioconductor.org/packages/release/bioc/html/GeneTonic.html). In our paper titled [`Interactive and Reproducible Workflows for Exploring and Modeling RNA-seq Data with pcaExplorer, Ideal, and GeneTonic`](https://doi.org/10.1002/cpz1.411), we provide a detailed description of the usage and workflow of these packages for RNA-seq. The general workflow is visually represented in the figure below, which is also part of Figure 1 in our paper.

<center>

![](Figure1.png){width=50%}
<br>
*Figure 1: Workflow schematic of an RNA-seq data analysis*
<br>
<span style="font-size: smaller;">The different data formats to provide to each of the applications are represented as tabular or list-like elements, named by the scheme followed in the packages. Dashed arrows indicate that the provided information can be used to generate or annotate another object. Solid arrows (gray) denote that an object has been derived/computed from the other where the connector originated from. The small boxes close to each element explain in which protocol each object is used as primary (black) or secondary (gray) input. On the right side, the two main approaches delivered by our software (interactivity via web applications and reproducibility via reporting) are represented. </span>

</center>
<br>
<br>

[`GeDi`](https://github.com/AnnekathrinSilvia/GeDi) complements this workflow as an additional step alongside [`GeneTonic`](https://bioconductor.org/packages/release/bioc/html/GeneTonic.html), specifically catering to the exploration of enrichment results. For the workflow to generate and prepare enrichment results, please refer to our linked paper or explore the vignette for detailed information. A brief summary is provided below. 

