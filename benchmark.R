library(drake)
library(dplyr)
library(hmeasure)
library(ranger)
library(biogram)
library(ggplot2)
library(tidyr)
library(pbapply)
requireNamespace("mlr3measures")

load("./data/benchmark_data.RData")
source("./functions/do_cv.R")
source("./functions/benchmark_functions.R")
source("./functions/test_alphabet.R")

if(Sys.info()[["nodename"]] %in% c("amyloid", "phobos", "huawei")) {
  data_path <- "/home/michal/Dropbox/AMP-analysis/AmpGram-analysis/"
}
if(Sys.info()[["nodename"]] %in% c("kasia-MACH-WX9", "ryzen")) {
  data_path <- "/home/kasia/Dropbox/AmpGram-analysis/"
}


benchmark_AmpGram <- drake_plan(full_benchmark_mer_preds = mutate(benchmark_mer_df,
                                                                  pred = predict(full_model_mers, 
                                                                                 data.frame(as.matrix(benchmark_ngrams)))[["predictions"]][, "TRUE"],
                                                                  target = ifelse(grepl("dbAMP", mer_id), "TRUE", "FALSE")),
                                full_benchmark_stats = calculate_statistics(full_benchmark_mer_preds) %>% 
                                  mutate(len_group = cut(n_peptide + 9, breaks = c(11, 19, 26, 36, 60, 710),
                                                         include.lowest = TRUE)),
                                full_benchmark_peptide_preds = mutate(full_benchmark_stats[, c(1:2,17)],
                                                                      Probability = predict(full_model_peptides, 
                                                                                            full_benchmark_stats[, 3:16])[["predictions"]][, "TRUE"],
                                                                      Decision = ifelse(Probability >= 0.5, TRUE, FALSE),
                                                                      Software = "AmpGram"),
                                len_groups = select(full_benchmark_stats, c("source_peptide", "len_group")),
                                all_benchmark_res = preprocess_benchmark_data(full_benchmark_peptide_preds, len_groups),
                                benchmark_summ = calculate_benchmark_summary(all_benchmark_res)
                                # Nobles_benchmark_datasets = preprocess_Nobles_datasets(),
                                # Nobles_datasets_preds = predict_Nobles_datasets(Nobles_benchmark_datasets,
                                #                                                 full_model_mers,
                                #                                                 imp_bigrams,
                                #                                                 full_model_peptides)
                                )

make(benchmark_AmpGram)

file.copy(from = ".drake", to = paste0(data_path, "drake-cache"), recursive = TRUE, overwrite = TRUE)
