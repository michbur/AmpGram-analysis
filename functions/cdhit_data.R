#' Filter sequences using CD-HIT.
#' 
#' Performs CD-HIT with 90% identity threshold to remove highly similar
#' sequences from the dataset. Only sequences that do not contain nonstandard
#' amino acids are used. The threshold value corresponds to the one 
#' used by Gabere M. N. & Noble W. S. when creating the validation dataset 
#' (Bioinformatics, 33(13), 1921–1929.).
#' 
#' @param seqs list of input sequences
#' @return list of sequences after filtering using CD-HIT

filter_cdhit <- function(seqs) {
  standard_seqs <- seqs[["standard"]]
  cdhit_res <- cdhit(standard_seqs, thresh = 0.90)
  print(paste0("Number of sequences after cd-hit: ", length(cdhit_res)))
  standard_seqs[cdhit_res] 
}

#' Perform CD-HIT clustering.
#'
#' Uses external software CD-HIT to cluster sequences into clusters 
#' that meet a given sequence identity threshold. Each cluster has 
#' one representative sequence. This function returns a vector of 
#' names of sequences that represent each cluster.
#' @param input_seq list of input sequences
#' @param thresh \code{double} indicating sequence identity threshold 
#' (e.g. 0.7 means 70% identity)
#' @param word_length word size
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
