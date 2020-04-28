library(dplyr)
library(drake)
library(ranger)
library(biogram)

source("./functions/get_mers.R")
source("./functions/count_ampgrams.R")
source("./functions/do_cv.R")
source("./functions/test_alphabet.R")

# Read in data
mer_df <- readd(mer_df)
binary_ngrams <- readd(binary_ngrams)

# Train full alphabet models
imp_bigrams <- calc_imp_bigrams(mer_df, binary_ngrams, c("[11,19]", "(19,26]"))
full_model_mers <- train_model_mers(mer_df, c("[11,19]", "(19,26]"), binary_ngrams, imp_bigrams)
mer_preds <- lapply(unique(mer_df[["group"]]), function(ith_group) {
  filter(mer_df, group == ith_group) %>% 
    mutate(pred = predict(full_model_mers, data.frame(as.matrix(binary_ngrams[mer_df[["group"]] == ith_group, imp_bigrams])))[["predictions"]][, "TRUE"])
}) %>% bind_rows()
full_model_peptides <- train_model_peptides(calculate_statistics(mer_preds))

# Train model using best alphabet
deg_binary_ngrams <- degenerate_ngrams(binary_ngrams, string2list("c_de_gw_hkr_afilmv_npqsty"), binarize = TRUE)
deg_imp_ngrams <- calc_imp_bigrams(mer_df, deg_binary_ngrams, c("[11,19]", "(19,26]"))
deg_model_mers <- train_model_mers(mer_df, c("[11,19]", "(19,26]"), deg_binary_ngrams, deg_imp_ngrams)
deg_mer_preds <- mutate(mer_df,
                        pred = predict(deg_model_mers, data.frame(as.matrix(deg_binary_ngrams[, deg_imp_ngrams])))[["predictions"]][, "TRUE"])
deg_model_peptides <- train_model_peptides(calculate_statistics(deg_mer_preds))

# Calculate ngrams for benchmark data
benchmark_mer_df <- read_fasta("./results/benchmark.fasta") %>%
  list2matrix() %>%
  create_mer_df()

benchmark_ngrams <- count_imp_ampgrams(benchmark_mer_df, imp_bigrams)
deg_benchmark_ngrams <- degenerate_ngrams(benchmark_ngrams, string2list("c_de_gw_hkr_afilmv_npqsty"), binarize = TRUE)

save(list = c("imp_bigrams", "full_model_mers", "full_model_peptides", "deg_binary_ngrams", "deg_imp_ngrams", "deg_model_mers", "deg_model_peptides", 
              "benchmark_mer_df", "benchmark_ngrams"),
     file = "./data/benchmark_data.RData")
