library(dplyr)
library(hmeasure)
library(ranger)
library(biogram)

load("./data/benchmark_data.RData")
source("./generate_benchmark_data.R")

benchmark_mer_preds <- mutate(benchmark_mer_df,
                              pred = predict(model_mers_full_alphabet, 
                                             data.frame(as.matrix(benchmark_ngrams[, which(benchmark_ngrams[["dimnames"]][[2]] %in% imp_bigrams)])))[["predictions"]][, "TRUE"],
                              target = ifelse(str_detect(mer_id, "dbAMP"), "TRUE", "FALSE"))

benchmark_stats <- calculate_statistics(benchmark_mer_preds)

benchmark_peptide_preds <- mutate(benchmark_stats[, 1:2],
                                  pred = predict(model_peptides_full_alphabet, 
                                                 benchmark_stats[, 3:16])[["predictions"]][, "TRUE"])
HMeasure(benchmark_peptide_preds[["target"]], benchmark_peptide_preds[["pred"]])[["metrics"]]
