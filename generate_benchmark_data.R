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

train_model_mers <- function(train_mer_df, train_groups, train_binary_ngrams, imp_bigrams) {
  train_dat <- filter(train_mer_df, group %in% train_groups)
  ranger_train_data <- data.frame(as.matrix(train_binary_ngrams[train_mer_df[["group"]] %in% train_groups, imp_bigrams]),
                                  tar = as.factor(train_dat[["target"]]))
  
  model_full_alphabet <- ranger(dependent.variable.name = "tar", data = ranger_train_data, 
                                write.forest = TRUE, probability = TRUE, num.trees = 2000, 
                                verbose = FALSE)
  model_full_alphabet
}

sort_group <- function(x) {
  splitted_x <- sapply(strsplit(x, split = ","), function(i) i[1])
  x_order <- order(as.numeric(gsub(pattern = "[^0-9]", 
                                   replacement = "", x = splitted_x)),
                   na.last = FALSE)
  x[x_order]
}

count_imp_ampgrams <- function(mer_df, imp_ampgrams) {
  mer_df[, grep("^X", colnames(mer_df))] %>% 
    as.matrix() %>% 
    count_specified(imp_ampgrams) %>% 
    binarize
}


# Read in data
mer_df <- readd(mer_df)
binary_ngrams <- readd(binary_ngrams)

# Train full alphabet models
imp_bigrams <- calc_imp_bigrams(mer_df, binary_ngrams, c("[11,19]", "(19,26]"))
full_model_mers <- train_model_mers(mer_df, c("[11,19]", "(19,26]"), binary_ngrams, imp_bigrams)
mer_preds <- mutate(mer_df,
                    pred = predict(full_model_mers, data.frame(as.matrix(binary_ngrams[, imp_bigrams])))[["predictions"]][, "TRUE"])
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

