#' Visualize the results of an enrichment analysis as word cloud
#'
#' Visualize the results of an enrichment analysis as a word cloud. The word
#' cloud highlights the most frequent terms associated with the description of
#' the genesets in the enrichment analysis.
#'
#' @param genesets_df A `data.frame` object of an enrichment analysis results.
#'                    This object should follow the input requirements of
#'                    \code{GeDi()}, check out the vignette for further details.
#'                    Besides the specified required columns, the object should
#'                    ideally include a column with a short geneset description
#'                    which is used for the word cloud. If no such column is
#'                    available, the row names of the `data.frame` are used for
#'                    the word cloud.
#'
#' @return A [wordcloud2::wordcloud2()] plot object
#' @export
#' @importFrom tm VCorpus VectorSource removeWords removePunctuation 
#' stripWhitespace stopwords TermDocumentMatrix tm_map
#' @importFrom wordcloud2 wordcloud2
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' ## Mock example showing how the data should look like
#'
#' ## If no "Term" or "Description" column is available,
#' ## the rownames of the data frame will be used.
#' geneset_df <- data.frame(
#'               Genesets = c("GO:0002503", "GO:0045087", "GO:0019886"),
#'               Genes = c("B2M, HLA-DMA, HLA-DMB",
#'                         "ACOD1, ADAM8, AIM2",
#'                         "B2M, CD74, CTSS")
#' )
#' rownames(geneset_df) <- geneset_df$Genesets
#'
#' wordcloud <- enrichmentWordcloud(geneset_df)
#'
#' ## With available "Term" column.
#' geneset_df <- data.frame(
#'               Genesets = c("GO:0002503", "GO:0045087", "GO:0019886"),
#'               Genes = c("B2M, HLA-DMA, HLA-DMB",
#'                         "ACOD1, ADAM8, AIM2",
#'                         "B2M, CD74, CTSS"),
#'               Term = c(
#'               "peptide antigen assembly with MHC class II protein complex",
#'               "innate immune response",
#'               "antigen processing and presentation of exogenous
#'                peptide antigen via MHC class II")
#' )
#'
#' wordcloud <- enrichmentWordcloud(geneset_df)
#'
#'
#' ## Example using the data available in the package
#'
#' data(macrophage_topGO_example,
#'      package = "GeDi",
#'      envir = environment())
#' wordcloud <- enrichmentWordcloud(macrophage_topGO_example)
enrichmentWordcloud <- function(genesets_df) {
  # Check if genesets are provided
  stopifnot(!is.null(genesets_df))

  # Check if there is a column named Term to use as terms for the wordcloud
  # If 'Term' column does not exist, check for a column named 'Description'
  # If neither 'Term' nor 'Description' columns exist, use row names as terms
  if ("Term" %in% names(genesets_df)) {
    terms <- genesets_df$Term
  } else if ("Description" %in% names(genesets_df)) {
    terms <- genesets_df$Description
  } else {
    terms <- rownames(genesets_df)
  }
  # Create a text corpus from the selected terms
  corpus <- VCorpus(VectorSource(terms))
  # Preprocess the text corpus by removing English stopwords, punctuation, 
  # and whitespace
  corpus <- tm_map(corpus, removeWords, stopwords("english"))
  corpus <- tm_map(corpus, removePunctuation)
  corpus <- tm_map(corpus, stripWhitespace)
  # Create a document-term matrix from the preprocessed corpus
  dtm <- as.matrix(TermDocumentMatrix(corpus))
  # Calculate the frequency of each term across documents
  v <- sort(rowSums(dtm), decreasing = TRUE)
  d <- data.frame(word = names(v), freq = v)
  # Generate a word cloud using the wordcloud2 package
  w <- wordcloud2(d,
                  size = 0.5,
                  color = brewer.pal(10, "PRGn"))

  # Return the generated word cloud
  return(w)
}
