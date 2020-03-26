
degenerate_ngrams <- function(x, element_groups) {
  
  deg_ngrams <- colnames(x) %>% 
    decode_ngrams %>% 
    strsplit("") %>% 
    lapply(degenerate, element_groups = element_groups) %>% 
    lapply(paste0, collapse = "") %>% 
    lapply(code_ngrams)
  
  res <- do.call(cbind, lapply(unique(deg_ngrams), function(ith_ngram) {
    row_sums(x[, ith_ngram == deg_ngrams, drop = FALSE])
  }))
  
  colnames(res) <- unique(deg_ngrams)
  
  res
}
