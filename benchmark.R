library(drake)
library(dplyr)
library(hmeasure)
library(ranger)
library(biogram)
library(ggplot2)
library(tidyr)
library(pbapply)
library(pROC)
library(xtable)
requireNamespace("mlr3measures")

load("./data/benchmark_data.RData")
source("./functions/do_cv.R")
source("./functions/benchmark_functions.R")
source("./functions/train_model_peptides.R")

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
                                benchmark_summ = calculate_benchmark_summary(all_benchmark_res),
                                Nobles_benchmark_datasets = preprocess_Nobles_datasets(),
                                Nobles_datasets_preds = readRDS("./results/Nobles_datasets_benchmark_res.rds"),
                                Nobles_APD_res = filter(Nobles_datasets_preds, grepl("APD", source_peptide)),
                                Nobles_DAMPD_res = filter(Nobles_datasets_preds, grepl("DAMPD", source_peptide)),
                                Nobles_APD_AUC = auc(Nobles_APD_res[["target"]], Nobles_APD_res[["Probability"]]),
                                Nobles_DAMPD_AUC = auc(Nobles_DAMPD_res[["target"]], Nobles_DAMPD_res[["Probability"]]),
                                ampscanner_res = read.csv("./data/ampscanner_noble.csv", stringsAsFactors = FALSE)[, c(1,3)] %>% 
                                  setNames(c("source_peptide", "Ampscanner")),
                                amp_only = filter(Nobles_benchmark_datasets, !is.na(AMP_target)),
                                Nobles_datasets_benchmark_res = aggregate_Nobles_datasets_benchmark_res(amp_only),
                                DAMPD_seqs_to_remove = find_apd_and_train_seqs(),
                                Nobles_datasets_benchmark_DAMPD_res = aggregate_Nobles_datasets_benchmark_res(filter(amp_only, !(`Peptide.sequence` %in% DAMPD_seqs_to_remove) & dataset == "DAMPD"))
                                )

make(benchmark_AmpGram)

file.copy(from = ".drake", to = paste0(data_path, "drake-cache"), recursive = TRUE, overwrite = TRUE)

benchmark_summ <- readd(benchmark_summ)
Nobles_datasets_benchmark_res <- readd(Nobles_datasets_benchmark_res)

write.csv(benchmark_summ, file = paste0(data_path, "publication-results/benchmark_results.csv"), row.names = FALSE)
write.csv(Nobles_datasets_benchmark_res, file = paste0(data_path, "publication-results/Nobles_benchmark_results.csv"), 
          row.names = FALSE)


# Generate tables with benchmark results on Noble's datasets
lapply(unique(Nobles_datasets_benchmark_res[["Dataset"]]), function(ith_set) {
  filter(Nobles_datasets_benchmark_res, Dataset == ith_set) %>% 
    ungroup %>% 
    select(c("Software", "AUC", "Precision", "Sensitivity", "Specificity")) %>% 
    format_table(caption = "", label = "", range=2L:5) %>% 
    writeLines(.,paste0(data_path, "publication-results/benchmark_Noble_", ith_set, ".txt"))
})


# Generate tables with benchmark results on our datasets
lapply(unique(benchmark_summ[["len_group"]]), function(ith_group) {
  print(ith_group)
  benchmark_summ[, c(7,2,1,4,5,6,8)] %>% 
    filter(len_group == ith_group) %>% 
    mutate(Software = ifelse(prob == "TRUE", paste0(Software), paste0(Software, "*"))) %>% 
    mutate(Software = relevel(as.factor(Software), "AmpGram")) %>% 
    select(-c(len_group, prob)) %>% 
    arrange(Software) %>% 
    format_table(caption = "", label = "", range=2L:5) %>% 
    writeLines(.,paste0(data_path, "publication-results/table_", ith_group, ".txt"))
})


# Generate table with benchmark on DAMPD dataset excluding sequences used for training AmpGram or AmpScanner
Nobles_datasets_benchmark_DAMPD_res %>% 
  ungroup %>% 
  select(c("Software", "AUC", "Precision", "Sensitivity", "Specificity")) %>% 
  format_table(caption = "", label = "", range=2L:5) %>% 
  writeLines(.,paste0(data_path, "publication-results/DAMPD_benchmark.txt"))
