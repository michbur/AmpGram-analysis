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

# AmpGram full alphabet
full_benchmark_mer_preds <- mutate(benchmark_mer_df,
                              pred = predict(full_model_mers, 
                                             data.frame(as.matrix(benchmark_ngrams)))[["predictions"]][, "TRUE"],
                              target = ifelse(grepl("dbAMP", mer_id), "TRUE", "FALSE")) 

full_benchmark_stats <- calculate_statistics(full_benchmark_mer_preds) %>% 
  mutate(len_group = cut(n_peptide + 9, breaks = c(11, 19, 26, 36, 60, 710),
                         include.lowest = TRUE))

full_benchmark_peptide_preds <- mutate(full_benchmark_stats[, c(1:2,17)],
                                  Probability = predict(full_model_peptides, 
                                                        full_benchmark_stats[, 3:16])[["predictions"]][, "TRUE"],
                                  Decision = ifelse(Probability >= 0.5, TRUE, FALSE),
                                  Software = "AmpGram_full") 
#HMeasure(benchmark_peptide_preds[["target"]], benchmark_peptide_preds[["Probability"]])[["metrics"]]

# AmpGram simplified alphabet
deg_benchmark_ngrams <- degenerate_ngrams(benchmark_ngrams, string2list("c_de_gw_hkr_afilmv_npqsty"), binarize = TRUE)
deg_benchmark_mer_preds <- mutate(benchmark_mer_df,
                                   pred = predict(deg_model_mers, 
                                                  data.frame(as.matrix(deg_benchmark_ngrams[, deg_imp_ngrams])))[["predictions"]][, "TRUE"],
                                   target = ifelse(grepl("dbAMP", mer_id), "TRUE", "FALSE")) 

deg_benchmark_stats <- calculate_statistics(deg_benchmark_mer_preds) %>% 
  mutate(len_group = cut(n_peptide + 9, breaks = c(11, 19, 26, 36, 60, 710),
                         include.lowest = TRUE))

deg_benchmark_peptide_preds <- mutate(deg_benchmark_stats[, c(1:2,17)],
                                       Probability = predict(deg_model_peptides, 
                                                             deg_benchmark_stats[, 3:16])[["predictions"]][, "TRUE"],
                                       Decision = ifelse(Probability >= 0.5, TRUE, FALSE),
                                       Software = "AmpGram_simplified") 


len_groups <- select(full_benchmark_stats, c("source_peptide", "len_group"))

iAMPpred <- read.delim("./data/iAMPpred_benchmark.csv")[,-1] %>% 
  setNames(c("source_peptide", "iAMPpred_antibact", "iAMPpred_antivir", "iAMPpred_antifung")) %>% 
  pivot_longer(c("iAMPpred_antibact", "iAMPpred_antibact", "iAMPpred_antifung"), names_to = "Software", values_to = "Probability")

all_benchmark_res <- read.csv("./data/benchmark_all.csv") %>% 
  setNames(c("Software", "source_peptide", "Decision", "Probability")) %>% 
  bind_rows(full_benchmark_peptide_preds[, c(1,4:6)]) %>% 
  bind_rows(iAMPpred) %>% 
  filter(Software != "Amylogram") %>% 
  mutate(source_peptide = gsub("DBAMP", "dbAMP_", source_peptide),
         source_peptide = gsub("dbAMP", "dbAMP_", source_peptide),
         source_peptide = gsub("__", "_", source_peptide))

all_benchmark_res[["Decision"]][all_benchmark_res[["Software"]] %in% c("iAMPpred_antibact", "iAMPpred_antivir", "iAMPpred_antifung")] <- ifelse(
  (all_benchmark_res[["Probability"]][all_benchmark_res[["Software"]] %in% c("iAMPpred_antibact", "iAMPpred_antivir", "iAMPpred_antifung")] >= 0.5), 
  TRUE, FALSE)

all_benchmark_res <- mutate(all_benchmark_res, target = ifelse(grepl("dbAMP", source_peptide), "TRUE", "FALSE")) %>% 
  left_join(len_groups, by = "source_peptide")

benchmark_summ <- lapply(unique(all_benchmark_res[["Software"]]), function(ith_software) {
  lapply(c(as.list(unique(all_benchmark_res[["len_group"]])), 
           list(unique(all_benchmark_res[["len_group"]]))), function(ith_length) {
             ith_dat <- filter(all_benchmark_res, Software == ith_software,
                               len_group %in% ith_length) %>% 
               mutate(target = as.factor(target),
                      Decision = as.factor(Decision))
             ith_res <- if(all(is.na(ith_dat[["Probability"]]))) {
             # For softwares without probabilities AUC is calculated with hmeasure using decision.
               data.frame(AUC = tryCatch(HMeasure(ith_dat[["target"]], as.logical(ith_dat[["Decision"]]))[["metrics"]][["AUC"]], 
                                              warning = function(w) {
                                                if(w[["message"]] == "ROC curve of scores mostly lying under the diagonal. Switching scores.") {
                                                  return(1 - HMeasure(ith_dat[["target"]], as.logical(ith_dat[["Decision"]]))[["metrics"]][["AUC"]])}
                                              }),
                      prob = FALSE)
               # For softwares with probabilities AUC is calculated with mlr3measures using probabilities
             } else {
               data.frame(AUC = mlr3measures::auc(ith_dat[["target"]], ith_dat[["Probability"]], "TRUE"),
                      prob = TRUE)
             }
             # Statistics below can be calculated only if there are two levels of decisions
             if(length(levels(ith_dat[["Decision"]])) == 2) {
               mutate(ith_res,
                      MCC = mlr3measures::mcc(ith_dat[["target"]], ith_dat[["Decision"]], "TRUE"),
                      precision = mlr3measures::precision(ith_dat[["target"]], ith_dat[["Decision"]], "TRUE"),
                      sensitivity = mlr3measures::sensitivity(ith_dat[["target"]], ith_dat[["Decision"]], "TRUE"),
                      specificity = mlr3measures::specificity(ith_dat[["target"]], ith_dat[["Decision"]], "TRUE"))
             } %>% mutate(Software = ith_software,
                          len_group = ifelse(length(ith_length) > 1, "all", as.character(ith_length)))
           }) %>% bind_rows()
}) %>% bind_rows() %>% 
  mutate(len_group = factor(len_group, levels = sort_group(unique(len_group))))


pivot_longer(benchmark_summ, c(AUC, MCC, precision, 
                               sensitivity, specificity)) %>% 
  mutate(AmpGram = Software == "AmpGram_full") %>% 
  ggplot(aes(x = Software, y = value, color = AmpGram)) +
  geom_point() +
  facet_grid(len_group ~ name, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90))



### Benchmark without peptides from APD

apd <- read.csv("./data/apd_df.csv", stringsAsFactors = FALSE)
all_benchmark_sequences <- read_fasta("./results/benchmark.fasta")
benchmark_names <- names(all_benchmark_sequences)[grepl('dbAMP', names(all_benchmark_sequences))]

benchmark_seqs <- lapply(1L:(length(all_benchmark_sequences)/2), function(i){
  paste(all_benchmark_sequences[[i]], collapse = "")
}) %>% unlist()

apd_seqs <- filter(apd, Sequence %in% benchmark_seqs) 
nrow(apd_seqs)

in_apd_names <- benchmark_names[which(benchmark_seqs %in% apd_seqs[["Sequence"]])]

filter(all_benchmark_res, !is.na(Probability) & !(source_peptide %in% in_apd_names)) %>% 
  group_by(Software) %>% 
  summarise(AUC = mlr3measures::auc(as.factor(target), Probability, "TRUE")) %>% 
  ggplot(aes(x = Software, y = AUC)) +
  geom_point()



### Benchmark on WS Noble's datasets


dampd <- read.delim("./data/SuppTable1.tsv", stringsAsFactors = FALSE)
apd <- read.delim("./data/SuppTable2.tsv", stringsAsFactors = FALSE)

both_datasets <- bind_rows(preprocess_dataset(dampd, "DAMPD"),
                           preprocess_dataset(apd, "APD"))

# Selecting only AMP datasets
datasets_amp_only <- filter(both_datasets, !is.na(AMP_target))
datasets_mer_preds <- pblapply(datasets_amp_only[["source_peptide"]], cl = 8, function(ith_peptide) {
  seq <- filter(datasets_amp_only, source_peptide == ith_peptide)[["Peptide.sequence"]]
  target <- filter(datasets_amp_only, source_peptide == ith_peptide)[["AMP_target"]]
  mers <- strsplit(seq, "")[[1]] %>% 
    matrix(nrow = 1) %>% 
    get_single_seq_mers() %>% 
    data.frame(stringsAsFactors = FALSE) %>% 
    mutate(source_peptide = ith_peptide,
           mer_id = paste0(source_peptide, "m", 1L:nrow(.)),
           target = target)
  counted_imp_ngrams <- count_imp_ampgrams(mers, imp_bigrams)
  res <- mutate(mers, 
         pred = predict(model_mers_full_alphabet, as.matrix(counted_imp_ngrams))[["predictions"]][,"TRUE"])
}) %>% 
  bind_rows()

datasets_stats <- calculate_statistics(datasets_mer_preds) 

datasets_peptide_preds <- mutate(datasets_stats,
                                  Probability = predict(model_peptides_full_alphabet, 
                                                        datasets_stats[, 3:16])[["predictions"]][, "TRUE"],
                                  Decision = ifelse(Probability >= 0.5, TRUE, FALSE),
                                  Software = "AmpGram_full") 

saveRDS(datasets_peptide_preds, file = "datasets_peptide_preds.rds")
