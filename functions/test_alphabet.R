count_longest <- function(x) {
  splitted_x <- strsplit(x = paste0(as.numeric(x > 0.5), collapse = ""),
                         split = "0")[[1]]
  len <- unname(sapply(splitted_x, nchar))
  if (length(len[len > 0]) == 0) {
    0 } else {
      len[len > 0]
    }
}

calculate_statistics <- function(pred_mers) {
  (if("fold" %in% colnames(pred_mers)) {
    group_by(pred_mers, source_peptide, target, fold)
  } else {  
    group_by(pred_mers, source_peptide, target)
  }) %>% 
    summarise(fraction_true = mean(pred > 0.5),
              pred_mean = mean(pred),
              pred_median = median(pred),
              n_peptide = length(pred),
              n_pos = sum(pred > 0.5),
              pred_min = min(pred),
              pred_max = max(pred), 
              longest_pos = max(count_longest(pred)),
              n_pos_10 = sum(count_longest(pred) >= 10),
              frac_0_0.2 = sum(pred <= 0.2)/n(),
              frac_0.2_0.4 = sum(pred > 0.2 & pred <= 0.4)/n(),
              frac_0.4_0.6 = sum(pred > 0.4 & pred <= 0.6)/n(),
              frac_0.6_0.8 = sum(pred > 0.6 & pred <= 0.8)/n(),
              frac_0.8_1 = sum(pred > 0.8 & pred <= 1)/n()) %>% 
    ungroup() %>% 
    mutate(target = factor(target))
}

train_model_peptides <- function(mer_statistics) {
  train_dat <- mer_statistics %>% 
    select(c("target", "fraction_true", "pred_mean", "pred_median",
             "n_peptide", "n_pos", "pred_min", "pred_max", "longest_pos",
             "n_pos_10", "frac_0_0.2", "frac_0.2_0.4", "frac_0.4_0.6",
             "frac_0.6_0.8", "frac_0.8_1"))
  model_cv <- ranger(dependent.variable.name = "target", data = train_dat, 
                     write.forest = TRUE, probability = TRUE, num.trees = 500, 
                     verbose = FALSE, classification = TRUE)
  model_cv
}


test_alphabet <- function(alphabet_file) {
  dat <- read.csv(alphabet_file, stringsAsFactors = FALSE)
  learner_groups <- list(`[11,19]` = 1L:175315, `(19,26]` = 175315L:350630, `[11,19]_(19,26]` = 350631L:525945)
  lapply(1L:length(learner_groups), function(i) {
    stats <- dat[learner_groups[[i]],] %>% 
      calculate_statistics() 
    lapply(1L:5, function(ith_fold) {
      test_dat <- filter(stats, fold == ith_fold)
      trained_model <- train_model_peptides(filter(stats, fold != ith_fold))
      preds <- mutate(test_dat,
                      pred = predict(trained_model, 
                                     data.frame(test_dat))[["predictions"]][, "TRUE"],
                      train_group = names(learner_groups[i]),
                      alphabet = dat[["alphabet"]][1])
    }) %>% bind_rows()
  }) %>% bind_rows()
}


test_all_alphabets <- function(data_path, alphabets) {
  lapply(alphabets, function(ith_alphabet) {
    test_alphabet(paste0(data_path, "results/", ith_alphabet, ".csv"))
  }) %>% bind_rows()
}

calc_measures_alphabets <- function(alphabets_preds) {
  alphabets_preds <- mutate(alphabets_preds, 
                            len_group = cut(n_peptide + 9, breaks = c(11, 19, 26, 36, 60, 710),
                                            include.lowest = TRUE))
  lapply(unique(alphabets_preds[["alphabet"]]), function(ith_alphabet) {
    lapply(unique(alphabets_preds[["train_group"]]), function(ith_train_group) {
      lapply(unique(alphabets_preds[["fold"]]), function(ith_fold) {
        lapply(unique(alphabets_preds[["len_group"]]), function(ith_group) {
          dat <- filter(alphabets_preds, fold == ith_fold & train_group == ith_train_group & alphabet == ith_alphabet & 
                          len_group == ith_group) %>% 
            mutate(decision = as.factor(ifelse(pred >= 0.5, "TRUE", "FALSE")))
          data.frame(alphabet = ith_alphabet,
                     train_group = ith_train_group,
                     fold = ith_fold,
                     len_group = ith_group,
                     AUC = mlr3measures::auc(dat[["target"]], dat[["pred"]], "TRUE"),
                     MCC = mlr3measures::mcc(dat[["target"]], dat[["decision"]], "TRUE"),
                     precision = mlr3measures::precision(dat[["target"]], dat[["decision"]], "TRUE"),
                     sensitivity = mlr3measures::sensitivity(dat[["target"]], dat[["decision"]], "TRUE"),
                     specificity = mlr3measures::specificity(dat[["target"]], dat[["decision"]], "TRUE"))
        }) %>% bind_rows()
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
}
