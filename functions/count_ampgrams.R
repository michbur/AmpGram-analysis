count_ampgrams <- function(mer_df, ns, ds) {
  
  mer_df[, grep("^X", colnames(mer_df))] %>% 
    as.matrix() %>% 
    count_multigrams(ns = ns, 
                     ds = ds,
                     seq = .,
                     u = toupper(colnames(aaprop))) %>% 
    binarize
}

# count_ampgrams <- function(mer_df, ns, ds) {
#   
#   mer_df[, grep("^X", colnames(mer_df))] %>% 
#     as.matrix() %>% 
#     count_multigrams(ns = c(1, rep(2, 4), rep(3, 4)), 
#                      ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
#                      seq = .,
#                      u = toupper(colnames(aaprop))) %>% 
#     binarize
# }
