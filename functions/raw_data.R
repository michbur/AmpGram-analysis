#' Reads in the raw AMP sequences for analysis
#' 
#' Reads in AMP sequences from the dbAMP database and extracts those
#' with experimental evidence for antimicrobial activity. Returns a list
#' of two lists: 
#' \itemize{
#'  \item{standard}{Sequences comprised of only standard amino acids}
#'  \item{non_standard}{Sequences containing nonstandard amino acids}
#'  }
read_raw_data <- function() 
  read.csv("./data/dbamp_df.csv") %>% 
  mutate(id = 1L:nrow(.)) %>% 
  mutate(id = sapply(id, paste_to_five)) %>% 
  mutate(id = paste0("dbAMP_", id)) %T>% {
  print(paste0("Number of sequences: ", nrow(.))) 
} %>%
  filter(Experimental.Evidence == "YES") %T>% {
    print(paste0("Number of sequences: ", nrow(.))) 
  } %>% {
    setNames(as.character(.[["Sequence"]]), .[["id"]])
  } %>% 
  strsplit("") %>% 
  purify() %T>% {
    print(paste0("Number of sequences with standard AA: ", length(.[["standard"]]))) 
    print(paste0("Number of sequences with non-standard AA: ", length(.[["non_standard"]]))) 
  }


#' Identifies sequences containing nonstandard amino acids.
#' 
#' This function checks the list of sequences for presence of nonstandard 
#' amino acids and returns a list of two lists:
#' \itemize{
#'  \item{standard}{Sequences comprised of only standard amino acids}
#'  \item{non_standard}{Sequences containing nonstandard amino acids}
#'  }
#' @param sequences list of marked sequences (each is list of character vector 
#'   of \code{sequence} and integer \code{target})
#' @return Input list of length 2 \code{sequences} excluding ones that are shorter than
#'   \code{min_length}, longer than \code{max_length} or has at least one 
#'   amino acid not in \code{a()[-1]}.

purify <- function(sequences) {
  standard <- toupper(biogram:::return_elements(seq_type = "prot"))
  is_standard <- vapply(sequences, function(seq) all(seq %in% standard), c(FALSE))
  
  list(standard = sequences[is_standard],
       non_standard = sequences[!is_standard])
}

#' function used to give dbAMP ids
paste_to_five <- function(x) {
  paste0(paste0(rep("0", 5 - nchar(x)), collapse = ""), x)
}

