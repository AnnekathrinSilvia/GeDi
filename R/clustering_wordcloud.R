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
#' @importFrom tm VCorpus VectorSource removeWords removePunctuation stripWhitespace stopwords TermDocumentMatrix tm_map
#' @importFrom wordcloud2 wordcloud2
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' data(macrophage_topGO_example,
#'      package = "GeDi",
#'      envir = environment())
#' wordcloud <- enrichmentWordcloud(macrophage_topGO_example)
enrichmentWordcloud <- function(genesets_df) {
  # Check if genesets are provided
  stopifnot(!is.null(genesets_df))

  # Check if there is a column named Term to use as terms for the wordcloud
  if ("Term" %in% names(genesets_df)) {
    terms <- genesets_df$Term
  }
  # If 'Term' column does not exist, check for a column named 'Description'
  else if ("Description" %in% names(genesets_df)) {
    terms <- genesets_df$Description
  }
  # If neither 'Term' nor 'Description' columns exist, use row names as terms
  else {
    terms <- rownames(genesets_df)
  }

  # Create a text corpus from the selected terms
  corpus <- VCorpus(VectorSource(terms))

  # Preprocess the text corpus by removing English stopwords, punctuation, and whitespace
  corpus <- tm_map(corpus, removeWords, stopwords("english"))
  corpus <- tm_map(corpus, removePunctuation)
  corpus <- tm_map(corpus, stripWhitespace)

  # Create a document-term matrix from the preprocessed corpus
  dtm <- as.matrix(TermDocumentMatrix(corpus))

  # Calculate the frequency of each term across documents
  v <- sort(rowSums(dtm), decreasing = TRUE)
  d <- data.frame(word = names(v), freq = v)

  # Set a fixed random seed for reproducibility
  set.seed(42)

  # Generate a word cloud using the wordcloud2 package
  w <- wordcloud2(d,
                  size = 0.5,
                  color = brewer.pal(10, "PRGn"))

  # Return the generated word cloud
  return(w)
}
