
#' Create all combinations of traits
#'
#' Creates all combinations of input traits.
#' @param vtraits a vector of trait indices in the expanded aaindex table
#' create by \code{\link{choose_properties}}.
#'
#' @return a list of traits combinations of given length (from single trait to
#' all traits).


create_traits_combination <- function(ftraits) {
  #vector of traits
  vtraits <- unlist(ftraits)
  
  #all combinations of traits
  #number of combinations: sum(sapply(all_traits_combn_list, nrow))
  pblapply(1L:length(vtraits), function(i)
    t(combn(vtraits, i)))
}



#' Create encodings
#'
#' Creates encodings (3-6 groups long) from list of traits.
#' @param ftraits a vector of trait indices in the expanded aaindex table
#' create by \code{\link{choose_properties}}.
#' @param list_duplicates if \code{TRUE} returns also a list of duplicates.
#' @return a named vector of encodings (for example 
#' \code{iknty_degpqrs_acfhlmvw})
create_encodings <- function(list_duplicates = FALSE) {
  
  paste_enc <- function(x)
    paste0(sapply(x, paste0, collapse = ""), collapse = "_")
  
  ftraits <- c("ARGP820101", "ARGP820103", 
               "BIGC670101", "KLEP840101", 
               "FAUJ880111", "FAUJ880112", 
               "FASG760104", "FASG760105", 
               "FASG760101", "KYTJ820101", 
               "BLAS910101")
  
  grouping_properties <- aaprop[ftraits, ]
  all_traits_combn_list <- create_traits_combination(ftraits)
  
  #create encodings
  all_aa_groups <- pblapply(3L:6, function(single_k) {
    res <- unlist(lapply(all_traits_combn_list, function(all_traits_combn)
      vapply(1L:nrow(all_traits_combn), function(single_trait_combn) {
        cl <- t(aaprop[unlist(all_traits_combn[single_trait_combn, , drop = FALSE]), , drop = FALSE]) %>%
          dist %>%
          hclust(method = "ward.D2")
        #cl <- hclust(dist(t(aaprop[unlist(all_traits_combn[single_trait_combn, , drop = FALSE]), , drop = FALSE])))
        gr <- cutree(cl, k = single_k)
        names(gr) <- tolower(names(gr))
        agg_gr <- lapply(unique(gr), function(single_group) names(gr[gr == single_group]))
        #inside encodings, amino acids are ordered alphabetically
        agg_gr <- lapply(agg_gr, sort)
        #groups are sorted by their length
        paste_enc(agg_gr[order(lengths(agg_gr))])
      }, "a")))
    names(res) <- paste0("ID", 1L:length(res), "K", single_k)
    res
  })
  
  #get indices of unique encodings
  aa_id <- lapply(all_aa_groups, function(i) !duplicated(i))
  
  aa_duplicates <- unlist(lapply(1L:length(aa_id), function(i) 
    lapply(all_aa_groups[[i]][aa_id[[i]]], function(j)
      names(which(j == all_aa_groups[[i]])))
  ), recursive = FALSE)
  #aa_duplicates <- aa_duplicates[lengths(aa_duplicates) > 1]
  
  #remove from aa_groups redundant encodings
  aa_groups <- unlist(lapply(1L:length(aa_id), function(i) {
    all_aa_groups[[i]][aa_id[[i]]]
  }), recursive = FALSE)
  
  #add as a benchmark two encodings from the literature
  aa1 = list(`1` = c("g", "a", "p", "v", "l", "i", "m"), 
             `2` = c("k", "r", "h"), 
             `3` = c("d", "e"), 
             `4` = c("f", "w", "y", "s", "t", "c", "n", "q"))
  
  aa2 = list(`1` = c("g", "a", "p", "v", "l", "i", "m", "f"), 
             `2` = c("k", "r", "h"), 
             `3` = c("d", "e"), 
             `4` = c("s", "t", "c", "n", "q", "y", "w"))
  
  if(list_duplicates) {
    list(aagroups = c(aa1 = paste_enc(aa1), aa2 = paste_enc(aa2), aa_groups),
         aa_duplicates = aa_duplicates)
  } else {
    c(aa1 = paste_enc(aa1), aa2 = paste_enc(aa2), aa_groups)
  }
}

create_encodings()




