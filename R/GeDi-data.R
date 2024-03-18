#' A sample input text file
#'
#' A sample input text file taken from the GScluster package
#'
#' @details This sample input text file contains data from the GScluster package.
#'          It is identical to the sample_geneset.txt file found on the Github
#'          page of the package.
#'
#' @references Yoon, S., Kim, J., Kim, SK. et al. GScluster: network-weighted
#'             gene-set clustering analysis. BMC Genomics 20, 352 (2019).
#'             https://doi.org/10.1186/s12864-019-5738-6
#'
#' @name sample_geneset
#' @docType data
NULL


#' An empty input text file
#'
#' An empty input text file to test the application
#'
#' @details This sample input text file is empty and used for testing the
#'          application.
#'
#' @name sample_geneset_empty
#' @docType data
NULL


#' A broken input text file
#'
#' A broken input text file to test the application
#'
#' @details This sample input text file is broken and used for testing the
#'          application.
#'
#' @name sample_geneset_broken
#' @docType data
NULL


#' A small sample input text file
#'
#' A sample input text file taken from the GScluster package, which is reduced
#' to a smaller number of entries for faster testing of the application.
#'
#' @details This sample input text file contains data from the GScluster package.
#'          It was taken from the sample_geneset.txt file found on the Github
#'          page of the package and then reduced to a smaller amount of entries
#'          for faster testing of the application.
#'
#' @references Yoon, S., Kim, J., Kim, SK. et al. GScluster: network-weighted
#'             gene-set clustering analysis. BMC Genomics 20, 352 (2019).
#'             https://doi.org/10.1186/s12864-019-5738-6
#'
#' @name sample_geneset_small
#' @docType data
NULL


#' A sample input RData file
#'
#' A sample input RData file generated from the macrophage dataset.
#'
#' @details This sample input contains data from the macrophage package found on
#'          Bioconductor. The exact steps used to generated this file can be
#'          found in the package vignette.
#'
#' @references Alasoo, K., Rodrigues, J., Mukhopadhyay, S. et al. Shared
#'             genetic effects on chromatin and gene expression indicate a role
#'             for enhancer priming in immune response. Nat Genet 50, 424–431
#'             (2018). https://doi.org/10.1038/s41588-018-0046-7
#'
#' @name macrophage_topGO_example
#' @docType data
NULL

#' A small sample input RData file
#'
#' A small sample input RData file generated from the macrophage dataset.
#'
#' @details This sample input contains data from the macrophage package found on
#'          Bioconductor. It is a small version of the
#'          `macrophage_topGO_example` and only contains the first 50 rows of
#'          this example. It can be used for fast testing of the application.
#'
#' @references Alasoo, K., Rodrigues, J., Mukhopadhyay, S. et al. Shared
#'             genetic effects on chromatin and gene expression indicate a role
#'             for enhancer priming in immune response. Nat Genet 50, 424–431
#'             (2018). https://doi.org/10.1038/s41588-018-0046-7
#'
#' @name macrophage_topGO_example_small
#' @docType data
NULL

#' Sample scores
#'
#' A file containing sample distance scores for the
#' `macrophage_topGO_example_small`.
#'
#' @details This sample input contains scores  for the
#'          `macrophage_topGO_example_small`. Distance scores have been
#'          calculated using the [GeDi::getJaccardMatrix()] method.
#'
#' @references Alasoo, K., Rodrigues, J., Mukhopadhyay, S. et al. Shared
#'             genetic effects on chromatin and gene expression indicate a role
#'             for enhancer priming in immune response. Nat Genet 50, 424–431
#'             (2018). https://doi.org/10.1038/s41588-018-0046-7
#'
#' @name scores_macrophage_topGO_example_small
#' @docType data
NULL


#' PPI
#'
#' A file containing a Protein-Protein Interaction (PPI) `data.frame` for the
#' `macrophage_topGO_example_small`.
#'
#' @details This sample input contains a PPI for the
#'          `macrophage_topGO_example_small`. The PPI has been downloaded using
#'          the functions to download a PPI matrix. Please check out the
#'          vignette for further information.
#'
#' @references Alasoo, K., Rodrigues, J., Mukhopadhyay, S. et al. Shared
#'             genetic effects on chromatin and gene expression indicate a role
#'             for enhancer priming in immune response. Nat Genet 50, 424–431
#'             (2018). https://doi.org/10.1038/s41588-018-0046-7
#'
#' @name ppi_macrophage_topGO_example_small
#' @docType data
NULL


#' A sample input RData file
#'
#' A sample input RData file generated from the macrophage dataset.
#'
#' @details This sample input contains data from the macrophage package found on
#'          Bioconductor. The exact steps used to generated this file can be
#'          found in the package vignette. The used database for the enrichment
#'          was the Reactome database.
#'
#' @references Alasoo, K., Rodrigues, J., Mukhopadhyay, S. et al. Shared
#'             genetic effects on chromatin and gene expression indicate a role
#'             for enhancer priming in immune response. Nat Genet 50, 424–431
#'             (2018). https://doi.org/10.1038/s41588-018-0046-7
#'
#' @name macrophage_Reactome_example
#' @docType data
NULL

#' A sample input RData file
#'
#' A sample input RData file generated from the macrophage dataset.
#'
#' @details This sample input contains data from the macrophage package found on
#'          Bioconductor. The exact steps used to generated this file can be
#'          found in the package vignette. The used database for the enrichment
#'          was the KEGG database.
#'
#' @references Alasoo, K., Rodrigues, J., Mukhopadhyay, S. et al. Shared
#'             genetic effects on chromatin and gene expression indicate a role
#'             for enhancer priming in immune response. Nat Genet 50, 424–431
#'             (2018). https://doi.org/10.1038/s41588-018-0046-7
#'
#' @name macrophage_KEGG_example
#' @docType data
NULL
