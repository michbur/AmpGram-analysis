#' Analyze AMPs containing nonstandard amino acids within a dataset
#' 
#' This function analyzes AMPs containing nonstandard amino acids within
#' their sequence. It returns letters corresponding to those amino acids
#' with a count of AMPs containing them.
#' 
#' @param raw_data Sequences read in using \code{\link{read_raw_data}} function
#' @return named integers indicating how many AMPs contain a given nonstandard
#' letter within their sequences
analyze_nonstandard_AMPs <- function(raw_data) {
  ns <- raw_data %>% 
    getElement("non_standard")
  
  nonstandard_aa <- unlist(ns) %>% 
    unique %>% 
    setdiff(toupper(biogram:::return_elements(seq_type = "prot"))) %>% 
    setNames(., .)
  
  lapply(nonstandard_aa, function(ith_aa) {
    ns[vapply(ns, function(seq) any(seq == ith_aa), c(FALSE))]
  }) %>% 
    lengths
}
