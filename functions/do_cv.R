# 
# mer_df <- readd(mer_df, cache = ampgram_cache)
# 
# binary_ngrams <- cbind(readd(ngrams12, cache = ampgram_cache),
#                        readd(ngrams3_1, cache = ampgram_cache),
#                        readd(ngrams3_2, cache = ampgram_cache))

do_cv <- function(mer_df, binary_ngrams)
  lapply(unique(mer_df[["group"]]), function(ith_group) 
    lapply(unique(mer_df[["fold"]]), function(ith_fold) {
      print(paste0(ith_group, "|", ith_fold))
      train_dat <- filter(mer_df, group == ith_group, fold != ith_fold)
      test_dat <- filter(mer_df, group == ith_group, fold == ith_fold)
      
      test_bis <- test_features(train_dat[["target"]],
                                binary_ngrams[mer_df[["group"]] == ith_group & 
                                                mer_df[["fold"]] != ith_fold, ])
      
      imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
      
      ranger_train_data <- data.frame(as.matrix(binary_ngrams[mer_df[["group"]] == ith_group & 
                                                                mer_df[["fold"]] != ith_fold, imp_bigrams]),
                                      tar = as.factor(train_dat[["target"]]))
      model_cv <- ranger(tar ~ ., ranger_train_data, write.forest = TRUE, probability = TRUE, 
                         num.trees = 2000, verbose = FALSE)
      
      preds <- mutate(test_dat,
                      pred = predict(model_cv, 
                                     data.frame(as.matrix(binary_ngrams[mer_df[["group"]] == ith_group & mer_df[["fold"]] == ith_fold, imp_bigrams])))[["predictions"]][, "TRUE"])
      
      # single mer predictions
      #HMeasure(preds[["target"]], preds[["pred"]])[["metrics"]]

      preds[, c("source_peptide", "mer_id", "group", "fold", "target", "pred")]
    }) %>% do.call(rbind, .)
  ) %>% do.call(rbind, .)
 
# peptide_preds <- group_by(preds, source_peptide, target) %>% 
#   summarise(peptide_pred = max(pred))
# rbind(HMeasure(as.numeric(peptide_preds[["target"]]), 
#                peptide_preds[["peptide_pred"]])[["metrics"]],
#       fold = ith_fold,
#       group = ith_group)
