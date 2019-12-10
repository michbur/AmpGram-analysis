filter_cdhit <- function(seqs) {
  standard_seqs <- seqs[["standard"]]
  cdhit_res <- cdhit(standard_seqs, thresh = 0.90)
  print(paste0("Number of sequences after cd-hit: ", length(cdhit_res)))
  standard_seqs[cdhit_res] 
}

#' Filter sequrences using cd-hit
#'
#' Filters sequences using external software cd-hit.
#' @param input_seq list of input sequences
#' @param threshold threshold value
#' @param cdhit_path \code{character} path to cd-hit
#' @return vector of names of filtered sequences

cdhit <- function(input_seq, thresh = 0.7, word_length = 2, roi_length = Inf,
                  cdhit_path = "./third-party") {
  
  input <- tempfile(tmpdir = getwd())
  output <- tempfile(tmpdir = getwd())
  cdhit <- paste0(cdhit_path, "/cdhit -i ", input,  " -o ", output, " -c ", thresh, " -n ", word_length)
  
  write_fasta(lapply(input_seq, 
                     function(single_seq) single_seq[1L:ifelse(length(single_seq) > roi_length, 
                                                               roi_length, length(single_seq))]), 
              input)
  #browser()
  system(cdhit)
  res <- read_fasta(output)
  file.remove(input, output, paste0(output, ".clstr"))
  names(res)
}
