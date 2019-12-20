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

read_and_cut <- function(path, lens) {
  seq_path <- paste0(path, "data/input-seqs.fasta")
  generate_cutted_sequences(read_fasta(seq_path), lens)
}

