library(dplyr)
library(drake)
library(ranger)
library(biogram)
library(pbapply)

source("./functions/get_mers.R")
source("./functions/count_ampgrams.R")
source("./functions/do_cv.R")
source("./functions/test_alphabet.R")
source("./functions/benchmark_functions.R")
  
# Read in data  
mer_df <- readd(mer_df)
binary_ngrams <- readd(binary_ngrams)

# Train full alphabet models
imp_bigrams <- calc_imp_bigrams(mer_df, binary_ngrams, c("[11,19]", "(19,26]"), 0.05)
full_model_mers <- train_model_mers(mer_df, c("[11,19]", "(19,26]"), binary_ngrams, imp_bigrams)
mer_preds <- pblapply(unique(mer_df[["group"]]), cl = 16, function(ith_group) {
  filter(mer_df, group == ith_group) %>% 
    mutate(pred = predict(full_model_mers, data.frame(as.matrix(binary_ngrams[mer_df[["group"]] == ith_group, imp_bigrams])))[["predictions"]][, "TRUE"])
}) %>% bind_rows()
full_model_peptides <- train_model_peptides(calculate_statistics(mer_preds))
AmpGram_model <- list("rf_mers" = full_model_mers, "rf_peptides" = full_model_peptides, "imp_features" = imp_bigrams)
class(AmpGram_model) <- "ag_model"
save(AmpGram_model, file = "./data/AmpGram_model.rda", compress = "xz", compression_level = 9)

# Train model using best alphabet
# deg_binary_ngrams <- degenerate_ngrams(binary_ngrams, string2list("c_de_gw_hkr_afilmv_npqsty"), binarize = TRUE)
# deg_imp_ngrams <- calc_imp_bigrams(mer_df, deg_binary_ngrams, c("[11,19]", "(19,26]"))
# deg_model_mers <- train_model_mers(mer_df, c("[11,19]", "(19,26]"), deg_binary_ngrams, deg_imp_ngrams)
# deg_mer_preds <- mutate(mer_df,
#                         pred = predict(deg_model_mers, data.frame(as.matrix(deg_binary_ngrams[, deg_imp_ngrams])))[["predictions"]][, "TRUE"])
# deg_model_peptides <- train_model_peptides(calculate_statistics(deg_mer_preds))

# Calculate ngrams for benchmark data
benchmark_mer_df <- read_fasta("./results/benchmark.fasta") %>%
  list2matrix() %>%
  create_mer_df()

benchmark_ngrams <- count_imp_ampgrams(benchmark_mer_df, imp_bigrams)
#deg_benchmark_ngrams <- degenerate_ngrams(benchmark_ngrams, string2list("c_de_gw_hkr_afilmv_npqsty"), binarize = TRUE)

# save(list = c("imp_bigrams", "full_model_mers", "full_model_peptides", "deg_binary_ngrams", "deg_imp_ngrams", "deg_model_mers", "deg_model_peptides", 
#               "benchmark_mer_df", "benchmark_ngrams"),
#      file = "./data/benchmark_data.RData")
save(list = c("imp_bigrams", "full_model_mers", "full_model_peptides", "benchmark_mer_df", "benchmark_ngrams"),
     file = "./data/benchmark_data.RData")

  
# Benchmark on WS Noble's datasets
Nobles_datasets <- preprocess_Nobles_datasets()
Nobles_datasets_mer_preds <- pblapply(Nobles_datasets[["source_peptide"]], cl = 16, function(ith_peptide) {
  seq <- filter(Nobles_datasets, source_peptide == ith_peptide)[["Peptide.sequence"]]
  target <- filter(Nobles_datasets, source_peptide == ith_peptide)[["AMP_target"]]
  mers <- strsplit(seq, "")[[1]] %>% 
    matrix(nrow = 1) %>% 
    get_single_seq_mers() %>% 
    data.frame(stringsAsFactors = FALSE) %>% 
    mutate(source_peptide = ith_peptide,
           mer_id = paste0(source_peptide, "m", 1L:nrow(.)),
           target = target)
  counted_imp_ngrams <- count_imp_ampgrams(mers, imp_bigrams)
  res <- mutate(mers, 
                pred = predict(full_model_mers, as.matrix(counted_imp_ngrams))[["predictions"]][,"TRUE"])
}) %>% bind_rows()

Nobles_datasets_stats <- calculate_statistics(Nobles_datasets_mer_preds) 

Nobles_datasets_peptide_preds <- mutate(Nobles_datasets_stats,
                                        Probability = predict(full_model_peptides, 
                                                              Nobles_datasets_stats[, 3:16])[["predictions"]][, "TRUE"],
                                        Decision = ifelse(Probability >= 0.5, TRUE, FALSE),
                                        Software = "AmpGram_full")
saveRDS(Nobles_datasets_peptide_preds, file = "./results/Nobles_datasets_benchmark_res.rds")
