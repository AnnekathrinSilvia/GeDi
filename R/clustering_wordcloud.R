#' Title
#'
#' @param genesets_df
#' @param cluster
#'
#' @return
#' @export
#' @import tm
#' @importFrom wordcloud2 wordcloud2
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
enrichmentWordcloud <- function(genesets_df){
  if("Term" %in% names(genesets_df)){
    terms <- genesets_df$Term
  }else if("Description" %in% names(genesets_df)){
    terms <- genesets_df$Description
  }else{
    terms <- rownames(genesets_df)
  }
  print(terms)

  corpus <- Corpus(VectorSource(terms))
  corpus <- tm_map(corpus, removeWords, stopwords("english"))
  corpus <- tm_map(corpus, removePunctuation)
  corpus <- tm_map(corpus, stripWhitespace)

  dtm <- as.matrix(TermDocumentMatrix(corpus))
  v <- sort(rowSums(dtm),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  print(d)

  set.seed(42)
  w <- wordcloud2(d,
                  size = 0.5,
            color = brewer.pal(10, "PRGn"))

  return(w)

}
