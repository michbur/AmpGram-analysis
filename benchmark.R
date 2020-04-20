library(dplyr)
library(hmeasure)
library(ranger)
library(biogram)
library(stringr)
library(ggplot2)
library(tidyr)
requireNamespace("mlr3measures")

load("./data/benchmark_data.RData")
source("./generate_benchmark_data.R")


benchmark_mer_preds <- mutate(benchmark_mer_df,
                              pred = predict(model_mers_full_alphabet, 
                                             data.frame(as.matrix(benchmark_ngrams[, which(benchmark_ngrams[["dimnames"]][[2]] %in% imp_bigrams)])))[["predictions"]][, "TRUE"],
                              target = ifelse(str_detect(mer_id, "dbAMP"), "TRUE", "FALSE")) 

benchmark_stats <- calculate_statistics(benchmark_mer_preds) %>% 
  mutate(len_group = cut(n_peptide + 9, breaks = c(11, 19, 26, 36, 60, 710),
                         include.lowest = TRUE))


benchmark_peptide_preds <- mutate(benchmark_stats[, c(1:2,17)],
                                  Probability = predict(model_peptides_full_alphabet, 
                                                        benchmark_stats[, 3:16])[["predictions"]][, "TRUE"],
                                  Decision = ifelse(Probability >= 0.5, TRUE, FALSE),
                                  Software = "AmpGram_full") 
#HMeasure(benchmark_peptide_preds[["target"]], benchmark_peptide_preds[["Probability"]])[["metrics"]]


len_groups <- select(benchmark_stats, c("source_peptide", "len_group"))

iAMPpred <- read.delim("./data/iAMPpred_benchmark.csv")[,-1] %>% 
  setNames(c("source_peptide", "iAMPpred_antibact", "iAMPpred_antivir", "iAMPpred_antifung")) %>% 
  tidyr::gather(Software, Probability, iAMPpred_antibact:iAMPpred_antifung)

all_benchmark_res <- read.csv("./data/benchmark_all.csv") %>% 
  setNames(c("Software", "source_peptide", "Decision", "Probability")) %>% 
  bind_rows(benchmark_peptide_preds[, c(1,4:6)]) %>% 
  bind_rows(iAMPpred) %>% 
  mutate(source_peptide = gsub("DBAMP", "dbAMP_", source_peptide),
         source_peptide = gsub("dbAMP", "dbAMP_", source_peptide),
         source_peptide = gsub("__", "_", source_peptide))

all_benchmark_res[["Decision"]][all_benchmark_res[["Software"]] %in% c("Amylogram", "iAMPpred_antibact", "iAMPpred_antivir", "iAMPpred_antifung")] <- ifelse(
  (all_benchmark_res[["Probability"]][all_benchmark_res[["Software"]] %in% c("Amylogram", "iAMPpred_antibact", "iAMPpred_antivir", "iAMPpred_antifung")] >= 0.5), 
  TRUE, FALSE)

all_benchmark_res <- mutate(all_benchmark_res, target = ifelse(str_detect(source_peptide, "dbAMP"), "TRUE", "FALSE")) %>% 
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
                                  return(1 - HMeasure(ith_dat[["target"]], as.logical(ith_dat[["Decision"]]))[["metrics"]][["AUC"]])
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
             Software = ith_software,
             len_group = ifelse(length(ith_length) > 1, "all", as.character(ith_length)),
             MCC = mlr3measures::mcc(ith_dat[["target"]], ith_dat[["Decision"]], "TRUE"),
             precision = mlr3measures::precision(ith_dat[["target"]], ith_dat[["Decision"]], "TRUE"),
             sensitivity = mlr3measures::sensitivity(ith_dat[["target"]], ith_dat[["Decision"]], "TRUE"),
             specificity = mlr3measures::specificity(ith_dat[["target"]], ith_dat[["Decision"]], "TRUE"))
    } else {
      mutate(ith_res,
             Software = ith_software,
             len_group = ifelse(length(ith_length) > 1, "all", as.character(ith_length))) 
    }
  }) %>% bind_rows()
}) %>% bind_rows()


pivot_longer(benchmark_summ, c(AUC, MCC, precision, 
                               sensitivity, specificity)) %>% 
  mutate(AmpGram = Software == "AmpGram_full") %>% 
  filter(Software != "Amylogram") %>% 
  ggplot(aes(x = Software, y = value, color = AmpGram)) +
  geom_point() +
  facet_grid(len_group ~ name, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90))

# AUC without len groups
filter(all_benchmark_res, !is.na(Probability)) %>% 
  group_by(Software) %>% 
  summarise(AUC = mlr3measures::auc(as.factor(target), Probability, "TRUE")) %>% 
  ggplot(aes(x = Software, y = AUC)) +
  geom_point()




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

