library(dplyr)
library(drake)
library(ranger)
library(biogram)

source("./functions/get_mers.R")
source("./functions/count_ampgrams.R")
source("./functions/do_cv.R")


count_longest <- function(x) {
  splitted_x <- strsplit(x = paste0(as.numeric(x > 0.5), collapse = ""),
                         split = "0")[[1]]
  len <- unname(sapply(splitted_x, nchar))
  if (length(len[len > 0]) == 0) {
    0 } else {
      len[len > 0]
    }
}



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

calculate_statistics <- function(pred_mers) {
  group_by(pred_mers, source_peptide, target) %>% 
    summarise(fraction_true = mean(pred > 0.5),
              pred_mean = mean(pred),
              pred_median = median(pred),
              n_peptide = length(pred),
              n_pos = sum(pred > 0.5),
              pred_min = min(pred),
              pred_max = max(pred), 
              longest_pos = max(count_longest(pred)),
              n_pos_10 = sum(count_longest(pred) >= 10),
              frac_0_0.2 = sum(pred <= 0.2)/n(),
              frac_0.2_0.4 = sum(pred > 0.2 & pred <= 0.4)/n(),
              frac_0.4_0.6 = sum(pred > 0.4 & pred <= 0.6)/n(),
              frac_0.6_0.8 = sum(pred > 0.6 & pred <= 0.8)/n(),
              frac_0.8_1 = sum(pred > 0.8 & pred <= 1)/n()) %>% 
    ungroup() %>% 
    mutate(target = factor(target))
}

train_model_peptides_full_alphabet <- function(mer_statistics) {
  train_dat <- mer_statistics %>% 
    select(c("target", "fraction_true", "pred_mean", "pred_median",
           "n_peptide", "n_pos", "pred_min", "pred_max", "longest_pos",
           "n_pos_10", "frac_0_0.2", "frac_0.2_0.4", "frac_0.4_0.6",
           "frac_0.6_0.8", "frac_0.8_1"))
  model_cv <- ranger(dependent.variable.name = "target", data = train_dat, 
                     write.forest = TRUE, probability = TRUE, num.trees = 500, 
                     verbose = FALSE, classification = TRUE)
  model_cv
}


# mer_preds <- do_cv(readd(mer_df), readd(binary_ngrams)) %>%
#   filter(source_file == "[11,19]_(19,26]")
# mer_preds <- all_cvs
# imp_bigrams <- calc_imp_bigrams(readd(mer_df), readd(binary_ngrams), c("[11,19]", "(19,26]"))
# model_mers_full_alphabet <- train_model_mers_full_alphabet(readd(mer_df), c("[11,19]", "(19,26]"), readd(binary_ngrams), imp_bigrams)
# 
# model_peptides_full_alphabet <- train_model_peptides_full_alphabet(calculate_statistics(mer_preds))
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

