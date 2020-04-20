library(dplyr)
library(drake)
library(ranger)
library(biogram)

source("./functions/get_mers.R")
source("./functions/count_ampgrams.R")
source("./functions/do_cv.R")
source("./functions/test_alphabet.R")


calc_imp_bigrams <- function(train_mer_df, train_binary_ngrams, train_groups) {
  train_dat <- filter(train_mer_df, group %in% train_groups)
  test_bis <- test_features(train_dat[["target"]],
                            train_binary_ngrams[train_mer_df[["group"]] %in% train_groups, ])
  imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
  imp_bigrams
}

train_model_mers_full_alphabet <- function(train_mer_df, train_groups, train_binary_ngrams, imp_bigrams) {
  train_dat <- filter(train_mer_df, group %in% train_groups)
  ranger_train_data <- data.frame(as.matrix(train_binary_ngrams[train_mer_df[["group"]] %in% train_groups, imp_bigrams]),
                                  tar = as.factor(train_dat[["target"]]))
  
  model_full_alphabet <- ranger(dependent.variable.name = "tar", data = ranger_train_data, 
                                write.forest = TRUE, probability = TRUE, num.trees = 2000, 
                                verbose = FALSE)
  model_full_alphabet
}



# mer_preds <- do_cv(readd(mer_df), readd(binary_ngrams)) %>%
#   filter(source_file == "[11,19]_(19,26]")
# imp_bigrams <- calc_imp_bigrams(readd(mer_df), readd(binary_ngrams), c("[11,19]", "(19,26]"))
# model_mers_full_alphabet <- train_model_mers_full_alphabet(readd(mer_df), c("[11,19]", "(19,26]"), readd(binary_ngrams), imp_bigrams)
# 
# model_peptides_full_alphabet <- train_model_peptides(calculate_statistics(mer_preds))
# 
# benchmark_mer_df <- read_fasta("./results/benchmark.fasta") %>%
#   list2matrix() %>%
#   create_mer_df()
# 
# benchmark_ngrams <- count_ampgrams(benchmark_mer_df,
#                                    ns = c(1, rep(2, 4), rep(3, 4)),
#                                    ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0), c(1, 1)))
# 
# save(list = c("mer_preds", "model_mers_full_alphabet", "imp_bigrams", "model_peptides_full_alphabet", "benchmark_mer_df", "benchmark_ngrams"),
#      file = "./data/benchmark_data.RData")

