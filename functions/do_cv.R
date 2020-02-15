
hexamer_df <- readd(mer_df, cache = ampgram_cache)

binary_ngrams <- cbind(readd(ngrams12, cache = ampgram_cache),
                       readd(ngrams3_1, cache = ampgram_cache),
                       readd(ngrams3_2, cache = ampgram_cache))

pred_df <- lapply(hexamer_df[["group"]], function(ith_group) 
  pblapply(hexamer_df[["fold"]], function(ith_fold) {
    train_dat <- filter(hexamer_df, group == ith_group, fold != ith_fold)
    test_dat <- filter(hexamer_df, group == ith_group, fold == ith_fold)
    
    test_bis <- test_features(train_dat[["target"]],
                              binary_ngrams[hexamer_df[["group"]] == ith_group & 
                                              hexamer_df[["fold"]] != ith_fold, ])
    
    imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
    
    ranger_train_data <- data.frame(as.matrix(binary_ngrams[hexamer_df[["group"]] == ith_group & 
                                                              hexamer_df[["fold"]] != ith_fold, imp_bigrams]),
                                    tar = as.factor(train_dat[["target"]]))
    model_cv <- ranger(tar ~ ., ranger_train_data, write.forest = TRUE, probability = TRUE, num.trees = 2000)
    
    preds <- mutate(test_dat,
                    pred = predict(model_cv, 
                                   data.frame(as.matrix(binary_ngrams[hexamer_df[["group"]] == ith_group & hexamer_df[["fold"]] == ith_fold, imp_bigrams])))[["predictions"]][, "TRUE"])
    
    # single mer predictions
    #HMeasure(preds[["target"]], preds[["pred"]])[["metrics"]]
    
    peptide_preds <- group_by(preds, source_peptide, target) %>% 
      summarise(peptide_pred = max(pred))
    
    rbind(HMeasure(as.numeric(peptide_preds[["target"]]), 
                   peptide_preds[["peptide_pred"]])[["metrics"]],
          fold = ith_fold,
          group = ith_group)
  }) %>% bind_rows()
) %>% bind_rows()

