library(dplyr)
library(hmeasure)
library(ranger)
library(biogram)
library(stringr)
library(ggplot2)

load("./data/benchmark_data.RData")
source("./generate_benchmark_data.R")

calc_mcc <- function(TP, TN, FP, FN)
  (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))

benchmark_mer_preds <- mutate(benchmark_mer_df,
                              pred = predict(model_mers_full_alphabet, 
                                             data.frame(as.matrix(benchmark_ngrams[, which(benchmark_ngrams[["dimnames"]][[2]] %in% imp_bigrams)])))[["predictions"]][, "TRUE"],
                              target = ifelse(str_detect(mer_id, "dbAMP"), "TRUE", "FALSE"))

benchmark_stats <- calculate_statistics(benchmark_mer_preds)

benchmark_peptide_preds <- mutate(benchmark_stats[, 1:2],
                                  Probability = predict(model_peptides_full_alphabet, 
                                                 benchmark_stats[, 3:16])[["predictions"]][, "TRUE"],
                                  Decision = ifelse(Probability >= 0.5, TRUE, FALSE),
                                  Software = "AmpGram_full")
HMeasure(benchmark_peptide_preds[["target"]], benchmark_peptide_preds[["Probability"]])[["metrics"]]



all_benchmark_res <- read.csv("./data/benchmark.csv") %>% 
  setNames(c("Software", "source_peptide", "Decision", "Probability")) %>% 
  bind_rows(benchmark_peptide_preds[, c(1,3:5)])

all_benchmark_res[["source_peptide"]] <- gsub("DBAMP", "dbAMP_", all_benchmark_res[["source_peptide"]])
all_benchmark_res[["Decision"]][all_benchmark_res[["Software"]] == 'Amylogram'] <- ifelse((all_benchmark_res[["Probability"]][all_benchmark_res[["Software"]] == 'Amylogram'] >= 0.5), TRUE, FALSE)

all_benchmark_res <- mutate(all_benchmark_res, target = ifelse(str_detect(source_peptide, "dbAMP"), "TRUE", "FALSE"))


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

ggplot(AUC, aes(x = Software, y = AUC)) +
  geom_point()

ggplot(MCC, aes(x = Software, y = MCC)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45))

ggplot(MCC, aes(x = Software, y = MCC)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45))