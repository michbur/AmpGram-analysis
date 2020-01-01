#' Generate holdout groups 
#' 
#' This function firstly classifies each sequence, depending on its 
#' length, into one of the five quantiles. Then, for each quantile, 
#' it randomly selects 1/10 of sequences for the benchmark dataset 
#' and the remaining 9/10 of the sequences is assigned to the train-test
#' dataset. 
#' 
#' @param sequences input sequences (the negative dataset)
#' @return List of five lists, each corresponding to one quantile and
#' comprising of two lists:
#' \itemize{
#'  \item{traintest}{Names of sequences selected for training and testing}
#'  \item{benchmark}{Names of sequences selected for benchmarking}
#'  }
generate_holdout_groups <- function(sequences) {
  seq_length_groups <- cut(lengths(sequences), 
                           breaks = as.numeric(quantile(lengths(sequences), probs = seq(0, 1, 0.2))),
                           include.lowest = TRUE)
  
  names(seq_length_groups) <- names(sequences)
  
  holdout_list <- lapply(levels(seq_length_groups), function(ith_group) {
    peptides_in_group <- names(seq_length_groups)[seq_length_groups == ith_group]
    group_benchmark <- sample(peptides_in_group, round(length(peptides_in_group)*0.10, 0))
    list(benchmark = group_benchmark,
         traintest = setdiff(peptides_in_group, group_benchmark))
  }) 
  
  names(holdout_list) <- levels(seq_length_groups)
  holdout_list
}

# set.seed(3)
# randomed_seqs <- setNames(lapply(runif(120, min = 9, max = 160), function(ith_len)
#     sample(toupper(biogram:::return_elements("prot")), ith_len, replace = TRUE)
#   ), paste0("P", 1L:120))

generate_holdout_groups(randomed_seqs)
