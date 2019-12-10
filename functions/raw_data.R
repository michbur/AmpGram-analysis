read_raw_data <- function() 
  read.csv("./data/dbamp_df.csv") %T>% {
  print(paste0("Number of sequences: ", nrow(.))) 
} %>%
  filter(Experimental.Evidence == "YES") %T>% {
    print(paste0("Number of sequences: ", nrow(.))) 
  } %>% 
  pull(Sequence) %>% 
  as.character() %>%
  setNames(., paste0("AMP", 1L:length(.))) %>% 
  strsplit("") %>% 
  purify() %T>% {
    print(paste0("Number of sequences with standard AA: ", length(.[["standard"]]))) 
    print(paste0("Number of sequences with non-standard AA: ", length(.[["non_standard"]]))) 
  }

# dla Kasi do poprawy
#' Removes sequences too long, too short and with nonstandard letters.
#' 
#' @param sequences list of marked sequences (each is list of character vector 
#'   of \code{sequence} and integer \code{target})
#' @return Input list of length 2 \code{sequences} excluding ones that are shorther than
#'   \code{min_length}, longer than \code{max_length} or has at least one 
#'   aminoacid not in \code{a()[-1]}.

purify <- function(sequences) {
  standard <- toupper(biogram:::return_elements(seq_type = "prot"))
  is_standard <- vapply(sequences, function(seq) all(seq %in% standard), c(FALSE))
  
  list(standard = sequences[is_standard],
       non_standard = sequences[!is_standard])
}
