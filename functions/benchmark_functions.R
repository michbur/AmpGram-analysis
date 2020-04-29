calc_imp_bigrams <- function(train_mer_df, train_binary_ngrams, train_groups) {
  train_dat <- filter(train_mer_df, group %in% train_groups)
  test_bis <- test_features(train_dat[["target"]],
                            train_binary_ngrams[train_mer_df[["group"]] %in% train_groups, ])
  imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
  imp_bigrams
}

train_model_mers <- function(train_mer_df, train_groups, train_binary_ngrams, imp_bigrams) {
  train_dat <- filter(train_mer_df, group %in% train_groups)
  ranger_train_data <- data.frame(as.matrix(train_binary_ngrams[train_mer_df[["group"]] %in% train_groups, imp_bigrams]),
                                  tar = as.factor(train_dat[["target"]]))
  
  model_full_alphabet <- ranger(dependent.variable.name = "tar", data = ranger_train_data, 
                                write.forest = TRUE, probability = TRUE, num.trees = 2000, 
                                verbose = FALSE)
  model_full_alphabet
}

sort_group <- function(x) {
  splitted_x <- sapply(strsplit(x, split = ","), function(i) i[1])
  x_order <- order(as.numeric(gsub(pattern = "[^0-9]", 
                                   replacement = "", x = splitted_x)),
                   na.last = FALSE)
  x[x_order]
}

count_imp_ampgrams <- function(mer_df, imp_ampgrams) {
  mer_df[, grep("^X", colnames(mer_df))] %>% 
    as.matrix() %>% 
    count_specified(imp_ampgrams) %>% 
    binarize
}

preprocess_dataset <- function(dat, prefix) {
  mutate(dat,
         dataset = prefix,
         source_peptide = paste0(prefix,  sprintf('%0.5d', 1:nrow(dat))),
         AMP_target = case_when(AMP.label == 1 ~ "TRUE",
                                AMP.label == -1 ~ "FALSE"),
         Antibacterial_target = case_when(Antibacterial.label == 1 ~ "TRUE",
                                          Antibacterial.label == -1 ~ "FALSE"),
         Bacteriocin_target = case_when(Bacteriocin.label == 1 ~ "TRUE",
                                        Bacteriocin.label == -1 ~ "FALSE"))
}


get_single_seq_mers <- function(seq) {
  seq2ngrams(seq, 10, a()[-1]) %>% 
    decode_ngrams() %>% 
    unname() %>% 
    strsplit(split = "") %>% 
    do.call(rbind, .) 
}


preprocess_benchmark_data <- function(full_benchmark_peptide_preds, len_groups) {
  iAMPpred <- read.delim("./data/iAMPpred_benchmark.csv")[,-1] %>% 
    setNames(c("source_peptide", "iAMPpred_antibact", "iAMPpred_antivir", "iAMPpred_antifung")) %>% 
    pivot_longer(c("iAMPpred_antibact", "iAMPpred_antibact", "iAMPpred_antifung"), names_to = "Software", values_to = "Probability")
  
  all_benchmark_res <- read.csv("./data/benchmark_all.csv") %>% 
    setNames(c("Software", "source_peptide", "Decision", "Probability")) %>% 
    bind_rows(full_benchmark_peptide_preds[, c(1,4:6)]) %>% 
    bind_rows(iAMPpred) %>% 
    filter(!(Software %in% c("Amylogram", "ADAM-HMM"))) %>% 
    mutate(source_peptide = gsub("DBAMP", "dbAMP_", source_peptide),
           source_peptide = gsub("dbAMP", "dbAMP_", source_peptide),
           source_peptide = gsub("__", "_", source_peptide))
  
  all_benchmark_res[["Decision"]][all_benchmark_res[["Software"]] %in% c("iAMPpred_antibact", "iAMPpred_antivir", "iAMPpred_antifung")] <- ifelse(
    (all_benchmark_res[["Probability"]][all_benchmark_res[["Software"]] %in% c("iAMPpred_antibact", "iAMPpred_antivir", "iAMPpred_antifung")] >= 0.5), 
    TRUE, FALSE)
  
  all_benchmark_res <- mutate(all_benchmark_res, target = ifelse(grepl("dbAMP", source_peptide), "TRUE", "FALSE")) %>% 
    left_join(len_groups, by = "source_peptide")
}

calculate_benchmark_summary <- function(all_benchmark_res) {
  lapply(unique(all_benchmark_res[["Software"]]), function(ith_software) {
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
}


preprocess_Nobles_datasets <- function() {
  dampd <- read.delim("./data/SuppTable1.tsv", stringsAsFactors = FALSE)
  apd <- read.delim("./data/SuppTable2.tsv", stringsAsFactors = FALSE)
  both_datasets <- bind_rows(preprocess_dataset(dampd, "DAMPD"),
                             preprocess_dataset(apd, "APD"))
  filter(both_datasets, !is.na(AMP_target))
}

predict_Nobles_datasets <- function(datasets, model_mers, imp_ngrams, model_peptides) {
  datasets_mer_preds <- pblapply(datasets[["source_peptide"]], cl = 16, function(ith_peptide) {
    seq <- filter(datasets, source_peptide == ith_peptide)[["Peptide.sequence"]]
    target <- filter(datasets, source_peptide == ith_peptide)[["AMP_target"]]
    mers <- strsplit(seq, "")[[1]] %>% 
      matrix(nrow = 1) %>% 
      get_single_seq_mers() %>% 
      data.frame(stringsAsFactors = FALSE) %>% 
      mutate(source_peptide = ith_peptide,
             mer_id = paste0(source_peptide, "m", 1L:nrow(.)),
             target = target)
    counted_imp_ngrams <- count_imp_ampgrams(mers, imp_ngrams)
    res <- mutate(mers, 
           pred = predict(model_mers, as.matrix(counted_imp_ngrams))[["predictions"]][,"TRUE"])
    res
  }) %>% bind_rows()
  
  datasets_stats <- calculate_statistics(datasets_mer_preds) 
  
  datasets_peptide_preds <- mutate(datasets_stats,
                                   Probability = predict(model_peptides, 
                                                         datasets_stats[, 3:16])[["predictions"]][, "TRUE"],
                                   Decision = ifelse(Probability >= 0.5, TRUE, FALSE),
                                   Software = "AmpGram_full")
}