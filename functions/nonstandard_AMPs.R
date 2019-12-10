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
