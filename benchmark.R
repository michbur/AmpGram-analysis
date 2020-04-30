library(drake)
library(dplyr)
library(hmeasure)
library(ranger)
library(biogram)
library(ggplot2)
library(tidyr)
library(pbapply)
library(pROC)
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
                                benchmark_summ = calculate_benchmark_summary(all_benchmark_res),
                                Nobles_benchmark_datasets = preprocess_Nobles_datasets(),
                                Nobles_datasets_preds = readRDS("./results/Nobles_datasets_benchmark_res.rds"),
                                Nobles_APD_res = filter(Nobles_datasets_preds, grepl("APD", source_peptide)),
                                Nobles_DAMPD_res = filter(Nobles_datasets_preds, grepl("DAMPD", source_peptide)),
                                Nobles_APD_AUC = auc(Nobles_APD_res[["target"]], Nobles_APD_res[["Probability"]]),
                                Nobles_DAMPD_AUC = auc(Nobles_DAMPD_res[["target"]], Nobles_DAMPD_res[["Probability"]]),
                                ampscanner_res = read.csv("./results/ampscanner_noble.csv", stringsAsFactors = FALSE)[, c(1,3)] %>% 
                                  setNames(c("source_peptide", "Ampscanner")),
                                amp_only = filter(Nobles_benchmark_datasets, !is.na(AMP_target)),
                                Nobles_datasets_benchmark_res =  left_join(amp_only[,c(4:9, 15:17)], ampscanner_res) %>% 
                                  left_join(Nobles_datasets_preds[,c("source_peptide", "Probability")]) %>% 
                                  mutate(`CAMP.RF..score` = ifelse(is.infinite(`CAMP.RF..score`), 0, `CAMP.RF..score`)) %>% 
                                  mutate(`ADAM.score` = ifelse(is.infinite(`ADAM.score`), -3, `ADAM.score`)) %>% 
                                  mutate(`CAMP.SVM..score` = ifelse(is.infinite(`CAMP.SVM..score`), 0, `CAMP.SVM..score`)) %>% 
                                  mutate(Ampscanner = ifelse(is.na(Ampscanner), 0, Ampscanner)) %>% 
                                  setNames(c("ADAM", "CAMP-RF", "CAMP-SVM", "DBAASP", "MLAMP", "AMPA", "Dataset", "source_peptide", "target", "AMPScanner V2", "AmpGram")) %>% 
                                  pivot_longer(c("ADAM", "CAMP-RF", "CAMP-SVM", "DBAASP", "MLAMP", "AMPA", "AMPScanner V2", "AmpGram"), names_to = "Software",
                                               values_to = "Probability") %>% 
                                  group_by(Dataset, Software) %>% 
                                  summarise(AUC = auc(target, Probability))
                                )

make(benchmark_AmpGram)

file.copy(from = ".drake", to = paste0(data_path, "drake-cache"), recursive = TRUE, overwrite = TRUE)

