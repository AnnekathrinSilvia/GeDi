#' Show the enriched terms as wordcloud
#'
#' Show the terms of the genesets descriptions as wordcloud of enriched terms
#'
#' @param genesets_df `data.frame`, a `data.frame` of the input data, containing
#'                    the genesets, genes and a description of the individual
#'                    genesets that can be used for the wordcloud. If no
#'                    description or term is included in the data, the rownames
#'                    will be used.
#'
#' @return a `list` of parameters which can be used to generate a wordcloud
#' @export
#' @import tm
#' @importFrom wordcloud2 wordcloud2
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' data(macrophage_topGO_example,
#'      package = "GeDi",
#'      envir = environment())
#' wordcloud <- enrichmentWordcloud(macrophage_topGO_example)
enrichmentWordcloud <- function(genesets_df){
  if("Term" %in% names(genesets_df)){
    terms <- genesets_df$Term
  }else if("Description" %in% names(genesets_df)){
    terms <- genesets_df$Description
  }else{
    terms <- rownames(genesets_df)
  }

  corpus <- Corpus(VectorSource(terms))
  corpus <- tm_map(corpus, removeWords, stopwords("english"))
  corpus <- tm_map(corpus, removePunctuation)
  corpus <- tm_map(corpus, stripWhitespace)

  dtm <- as.matrix(TermDocumentMatrix(corpus))
  v <- sort(rowSums(dtm),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)

  set.seed(42)
  w <- wordcloud2(d,
                  size = 0.5,
            color = brewer.pal(10, "PRGn"))

  return(w)

}
