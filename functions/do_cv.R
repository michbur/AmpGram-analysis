# 
# mer_df <- readd(mer_df, cache = ampgram_cache)
# 
# binary_ngrams <- cbind(readd(ngrams12, cache = ampgram_cache),
#                        readd(ngrams3_1, cache = ampgram_cache),
#                        readd(ngrams3_2, cache = ampgram_cache))

sort_group <- function(x) {
  splitted_x <- sapply(strsplit(x, split = ","), function(i) i[1])
  x_order <- order(as.numeric(gsub(pattern = "[^0-9]", 
                                   replacement = "", x = splitted_x)))
  x[x_order]
}

string2list <- function(x) {
  pasted_group <- strsplit(x, "_", fixed = TRUE)[[1]] %>% 
    toupper()
  res <- strsplit(pasted_group, "")
  names(res) <- 1L:length(res)
  res
}



do_cv <- function(mer_df, binary_ngrams) {
  possible_groups <- sort_group(unique(mer_df[["group"]]))[1L:2]
  group_combs <- list(possible_groups[1], possible_groups[2], c(possible_groups))
  lapply(group_combs, function(ith_group) {
    lapply(unique(mer_df[["fold"]]), function(ith_fold) {
      print(paste0(ith_group, "|", ith_fold))
      train_dat <- filter(mer_df, group %in% ith_group, fold != ith_fold)
      test_dat <- filter(mer_df, fold == ith_fold)
      
      test_bis <- test_features(train_dat[["target"]],
                                binary_ngrams[mer_df[["group"]] %in% ith_group & 
                                                mer_df[["fold"]] != ith_fold, ])
      
      imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
      
      ranger_train_data <- data.frame(as.matrix(binary_ngrams[mer_df[["group"]] %in% ith_group & 
                                                                mer_df[["fold"]] != ith_fold, imp_bigrams]),
                                      tar = as.factor(train_dat[["target"]]))
      model_cv <- ranger(dependent.variable.name = "tar", data =  ranger_train_data, 
                         write.forest = TRUE, probability = TRUE, num.trees = 2000, 
                         verbose = FALSE)
      
      preds <- mutate(test_dat,
                      pred = predict(model_cv, 
                                     data.frame(as.matrix(binary_ngrams[mer_df[["fold"]] == ith_fold, imp_bigrams])))[["predictions"]][, "TRUE"])
      
      # single mer predictions
      #HMeasure(preds[["target"]], preds[["pred"]])[["metrics"]]
      
      preds[, c("source_peptide", "mer_id", "group", 
                "fold", "target", "pred")] 
    }) %>% bind_rows()
  }) %>% bind_rows()
}

