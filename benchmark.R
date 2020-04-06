library(dplyr)
library(hmeasure)
library(ranger)
library(biogram)
library(stringr)
library(ggplot2)
library(tidyr)

load("./data/benchmark_data.RData")
source("./generate_benchmark_data.R")

calc_mcc <- function(TP, TN, FP, FN)
  (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))

benchmark_mer_preds <- mutate(benchmark_mer_df,
                              pred = predict(model_mers_full_alphabet, 
                                             data.frame(as.matrix(benchmark_ngrams[, which(benchmark_ngrams[["dimnames"]][[2]] %in% imp_bigrams)])))[["predictions"]][, "TRUE"],
                              target = ifelse(str_detect(mer_id, "dbAMP"), "TRUE", "FALSE")) 

benchmark_stats <- calculate_statistics(benchmark_mer_preds) %>% 
  mutate(len_group = cut(n_peptide + 9, breaks = c(0, 20, 50, 100, 200, 800),
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
  gather(Software, Probability, iAMPpred_antibact:iAMPpred_antifung)

all_benchmark_res <- read.csv("./data/benchmark_all.csv") %>% 
  setNames(c("Software", "source_peptide", "Decision", "Probability")) %>% 
  bind_rows(benchmark_peptide_preds[, c(1,4:6)]) %>% 
  bind_rows(iAMPpred)

all_benchmark_res[["source_peptide"]] <- gsub("DBAMP", "dbAMP_", all_benchmark_res[["source_peptide"]])
all_benchmark_res[["Decision"]][all_benchmark_res[["Software"]] %in% c("Amylogram", "iAMPpred_antibact", "iAMPpred_antivir", "iAMPpred_antifung")] <- ifelse(
  (all_benchmark_res[["Probability"]][all_benchmark_res[["Software"]] %in% c("Amylogram", "iAMPpred_antibact", "iAMPpred_antivir", "iAMPpred_antifung")] >= 0.5), 
  TRUE, FALSE)

all_benchmark_res <- mutate(all_benchmark_res, target = ifelse(str_detect(source_peptide, "dbAMP"), "TRUE", "FALSE")) %>% 
  left_join(len_groups, by = "source_peptide")


MCC <- group_by(all_benchmark_res, Software) %>% 
  summarise(TP = as.numeric(sum(Decision == TRUE & target == "TRUE")),
            TN = as.numeric(sum(Decision == FALSE & target == "FALSE")),
            FP = as.numeric(sum(Decision == TRUE & target == "FALSE")),
            FN = as.numeric(sum(Decision == FALSE & target == "TRUE"))) %>% 
  mutate(MCC = (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)),
         Sensitivity = TP/(TP+FN),
         Specificity = TN/(TN+FP))

AUC <- filter(all_benchmark_res, !is.na(Probability)) %>% 
  group_by(Software) %>% 
  summarise(AUC = HMeasure(target, Probability)[["metrics"]][["AUC"]])


plot_data <- left_join(MCC, AUC, by = "Software") %>% 
  filter(Software != 'Amylogram') %>% 
  gather(Measure, Value, TP:AUC) 

ggplot(filter(plot_data, Measure == "AUC" & !is.na(Value)), aes(x = Software, y = Value)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(y = "AUC")

ggplot(filter(plot_data, Measure %in% c("MCC", "Sensitivity", "Specificity")), aes(x = Software, y = Value, color = Measure)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45))


# measures by len_group
MCC_lens <- group_by(all_benchmark_res, Software, len_group) %>% 
  summarise(TP = as.numeric(sum(Decision == TRUE & target == "TRUE")),
            TN = as.numeric(sum(Decision == FALSE & target == "FALSE")),
            FP = as.numeric(sum(Decision == TRUE & target == "FALSE")),
            FN = as.numeric(sum(Decision == FALSE & target == "TRUE"))) %>% 
  mutate(MCC = (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)),
         Sensitivity = TP/(TP+FN),
         Specificity = TN/(TN+FP)) %>% 
  filter(Software != 'Amylogram') %>% 
  gather(Measure, Value, TP:Specificity) 

AUC_lens <- filter(all_benchmark_res, !is.na(Probability)) %>% 
  group_by(Software, len_group) %>% 
  summarise(AUC = HMeasure(target, Probability)[["metrics"]][["AUC"]])

ggplot(AUC_lens, aes(x = Software, y = AUC)) +
  geom_point() +
  facet_wrap(~len_group, ncol = 2) +
  theme(axis.text.x = element_text(angle = 45))

ggplot(filter(MCC_lens, Measure %in% c("MCC", "Specificity", "Sensitivity")), aes(x = Software, y = Value, color = Measure)) +
  geom_point() +
  facet_wrap(~len_group, ncol = 2) +
  theme(axis.text.x = element_text(angle = 45))
