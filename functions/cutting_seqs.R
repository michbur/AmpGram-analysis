#' Generates cutted sequences for construction of a negative dataset.
#' 
#' This function performs generating of cutted sequences to construct
#' the negative dataset. First, it combines all the input sequences 
#' into the one character vector. Then, it splits splits the vector 
#' into parts of length equal to the length of the positive dataset 
#' (AMPs after CD-HIT). Each part is further splitted into smaller 
#' subsequences of lengths the same as lengths of sequences in the 
#' positive dataset. To ensure the same number of sequences of the same 
#' length in both positive and negative datasets, cutted sequences are
#' selected by choosing each subsequence (lengths defined by AMPs) 
#' from a randomly selected part of the vector.
#'
#' @param sequences input sequences
#' @param lens lenghts of sequences in the positive dataset
#' @return List of length equal to the length of positive dataset containing
#' cutted sequences 
generate_cutted_sequences <- function(sequences, lens) {
  seq_vec <- unlist(sequences, use.names = FALSE)
  
  end_vec <- unname(cumsum(sample(lens)))
  start_vec <- c(1, end_vec[-length(end_vec)] + 1)
  
  seq_parts <- floor(1L:length(seq_vec)/max(end_vec)) + 1
  max_part <- length(seq_vec)%/%max(end_vec) 
  
  splits <- split(seq_vec, seq_parts)
  pos_df <- data.frame(start = start_vec, end = end_vec, part = sample(1L:max_part, size = length(lens), replace = TRUE))
  cutted_sequences <- pblapply(1L:nrow(pos_df), function(ith_row) {
    unname(splits[[pos_df[ith_row, "part"]]][pos_df[ith_row, "start"]:pos_df[ith_row, "end"]])
  })
  
  names(cutted_sequences) <- paste0("CUTTED", 1L:length(cutted_sequences))
  cutted_sequences
}

# set.seed(3)
# randomed_seqs <- lapply(runif(16, min = 9, max = 160), function(ith_len)
#   sample(toupper(biogram:::return_elements("prot")), ith_len, replace = TRUE)
# )
# randomed_lens <- floor(runif(5, min = 9, max = 50))
# 
# generate_cutted_sequences(randomed_seqs, randomed_lens)

#' Generate the negative dataset
#' 
#' Reads in the input sequences and uses \code{\link{generate_cutted_sequences}}
#' function for construction of the negative dataset.
read_and_cut <- function(path, lens) {
  seq_path <- paste0(path, "data/input-seqs.fasta")
  raw_seqs <- read_fasta(seq_path)
  generate_cutted_sequences(purify(raw_seqs)[["standard"]], lens)
}

