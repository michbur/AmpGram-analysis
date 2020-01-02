create_hv <- function(seq) 
  do.call(rbind, lapply(1L:nrow(seq), function(i) {
    seq2ngrams(seq[i, ][!is.na(seq[i, ])], 6, a()[-1]) %>% 
      decode_ngrams() %>% 
      unname() %>% 
      strsplit(split = "") %>% 
      do.call(rbind, .) %>% 
      data.frame(stringsAsFactors = FALSE) %>% 
      mutate(source_peptide = rownames(seq)[i],
             hexapeptide_id = paste0(source_peptide, "h", 1L:nrow(.)))
  }))

pos_seq <- readd(cdhit_data)
neg_seq <- readd(negative_data)
pos_seq_ids <- readd(cdhit_data_ids)
neg_seq_ids <- readd(negative_data_ids)

seq_groups <- lapply(names(pos_seq_ids), function(i)
  c(pos_seq[pos_seq_ids[[i]][["traintest"]]], neg_seq[neg_seq_ids[[i]][["traintest"]]])) %>% 
    setNames(names(pos_seq_ids))
#lapply(pos_seq_ids, function(i) i[["traintest"]])
#lapply(neg_seq_ids, function(i) i[["traintest"]])

hexamer_df <- lapply(names(seq_groups), function(ith_group_id) {
  ith_group <- seq_groups[[ith_group_id]]
  
  folded <- cvFolds(length(ith_group), K = 5)
  fold_df <- data.frame(source_peptide = names(ith_group)[folded[["subsets"]]], 
                        which = folded[["which"]],
                        stringsAsFactors = FALSE)
  
  ith_group %>% 
    list2matrix() %>% 
    create_hv %>% 
    mutate(group = ith_group_id) %>% 
    inner_join(fold_df, by = c("source_peptide" = "source_peptide"))
}) %>% 
  do.call(rbind, .) %>% 
  mutate(target = grepl("AMP", source_peptide, fixed = TRUE))


ngrams <- hexamer_df[, 1L:6] %>% 
  as.matrix() %>% 
  count_multigrams(ns = c(1, rep(2, 2)), 
                   ds = list(0, 0, 1),
                   seq = .,
                   u = toupper(colnames(aaprop))) %>% 
  binarize

binary_ngrams <- binarize(ngrams)


ith_group <- unique(hexamer_df[["group"]])[1]
ith_fold <- sort(unique(hexamer_df[["which"]]))[1]

train_dat <- filter(hexamer_df, group == ith_group, which != ith_fold)
test_dat <- filter(hexamer_df, group == ith_group, which == ith_fold)

test_bis <- test_features(train_dat[["target"]],
                          binary_ngrams[hexamer_df[["group"]] == ith_group & hexamer_df[["which"]] != ith_fold, ])

imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]

ranger_train_data <- data.frame(as.matrix(binary_ngrams[hexamer_df[["group"]] == ith_group & hexamer_df[["which"]] != ith_fold, imp_bigrams]),
                                tar = as.factor(train_dat[["target"]]))
model_cv <- ranger(tar ~ ., ranger_train_data, write.forest = TRUE, probability = TRUE, num.trees = 1000)

preds <- mutate(test_dat,
                pred = predict(model_cv, 
                               data.frame(as.matrix(binary_ngrams[hexamer_df[["group"]] == ith_group & hexamer_df[["which"]] == ith_fold, imp_bigrams])))[["predictions"]][, "TRUE"])

HMeasure(preds[["target"]], preds[["pred"]])[["metrics"]]


hexamer_df[, 1L:6] %>% 
  as.matrix() %>% 
  count_multigrams(ns = c(1, rep(2, 4), rep(3, 3)), 
                   ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                   seq = .,
                   u = toupper(colnames(aaprop)))




create_all_folds <- function(ets, seq_lengths) {
  seq_label <- cut(seq_lengths, breaks = c(5, 6, 10, 15, 25))
  
  #set the seed to create the same folds
  set.seed(1)
  lapply(1L:3, function(constant_pos) {
    lapply(1L:3, function(constant_neg) {
      
      lapply(1L:5, function(dummy) {
        splitted_seqs <- lapply(levels(seq_label), function(single_label) {
          pos_seqs <- which(ets == 1 & seq_label == single_label)
          neg_seqs <- which(ets == 0 & seq_label == single_label)
          
          fold_list <- lapply(list(pos = pos_seqs, neg = neg_seqs), function(single_n) {
            folded <- cvFolds(length(single_n), K = 5)
            data.frame(id = single_n[folded[["subsets"]]], which = folded[["which"]])
          })
        })
        
        list(pos_train = do.call(rbind, lapply(splitted_seqs[1L:constant_pos], 
                                               function(single_split) single_split[["pos"]])),
             pos_test = do.call(rbind, lapply(splitted_seqs[-(1L:constant_pos)], 
                                              function(single_split) single_split[["pos"]])),
             neg_train = do.call(rbind, lapply(splitted_seqs[1L:constant_pos], 
                                               function(single_split) single_split[["neg"]])),
             neg_test = do.call(rbind, lapply(splitted_seqs[-(1L:constant_pos)], 
                                              function(single_split) single_split[["neg"]])))
      })
    })
  }) %>% unlist(recursive = FALSE) %>% 
    unlist(recursive = FALSE)
}



do_cv <- function(all_folds, extracted_ngrams, hv)
  lapply(extracted_ngrams, function(encoded_group)
    lapply(all_folds, function(fold_list) 
      do_single_cv(fold_list, encoded_group, hv)
    )
  )


do_single_cv <- function(fold_list, encoded_group, hv) {
  lapply(1L:5, function(fold) {
    #training data
    train_pos <- encoded_group[hv %in% fold_list[[1]][fold_list[[1]][, "which"] != fold, "id"], ]
    train_neg <- encoded_group[hv %in% fold_list[[3]][fold_list[[3]][, "which"] != fold, "id"], ]
    
    test_bis <- test_features(c(rep(1, nrow(train_pos)), rep(0, nrow(train_neg))),
                              rbind(train_pos, train_neg), adjust = NULL)
    imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
    
    train_data <- data.frame(as.matrix(rbind(train_pos, train_neg)[, imp_bigrams]),
                             tar = as.factor(c(rep(1, nrow(train_pos)), rep(0, nrow(train_neg)))))
    model_cv <- ranger(tar ~ ., train_data, write.forest = TRUE, probability = TRUE)
    
    #test data
    test_pos <- rbind(encoded_group[hv %in% fold_list[[1]][fold_list[[1]][, "which"] == fold, "id"], ],
                      encoded_group[hv %in% fold_list[[2]][fold_list[[2]][, "which"] == fold, "id"], ])
    test_neg <- rbind(encoded_group[hv %in% fold_list[[3]][fold_list[[3]][, "which"] == fold, "id"], ], 
                      encoded_group[hv %in% fold_list[[4]][fold_list[[4]][, "which"] == fold, "id"], ])
    
    #number of n-grams from protein
    ngram_prots_pos <- c(as.vector(table(hv[hv %in% fold_list[[1]][fold_list[[1]][, "which"] == fold, "id"]])),
                         as.vector(table(hv[hv %in% fold_list[[2]][fold_list[[2]][, "which"] == fold, "id"]])))
    ngram_prots_neg <- c(as.vector(table(hv[hv %in% fold_list[[3]][fold_list[[3]][, "which"] == fold, "id"]])),
                         as.vector(table(hv[hv %in% fold_list[[4]][fold_list[[4]][, "which"] == fold, "id"]])))
    
    preds <- cbind(predict(model_cv, data.frame(as.matrix(rbind(test_pos, test_neg)[, imp_bigrams])))[["predictions"]][, 2], 
                   # probability of being amyloid
                   c(rep(1, nrow(test_pos)), rep(0, nrow(test_neg))),
                   c(unlist(lapply(1L:length(ngram_prots_pos), function(prot_id)
                     rep(prot_id, ngram_prots_pos[prot_id]))), 
                     unlist(lapply(1L:length(ngram_prots_neg), function(prot_id)
                       rep(prot_id, ngram_prots_neg[prot_id])))))
    
    # is randomForest consistent with ranger?
    # model_rf <- randomForest(tar ~ ., train_data)
    # predict(model_rf, data.frame(rbind(test_pos, test_neg)[, imp_bigrams])) == 
    #   predict(model_cv, data.frame(rbind(test_pos, test_neg)[, imp_bigrams]))[["predictions"]]
    
    preds_df <- preds %>% 
      data.frame %>% 
      rename(prob = X1, tar = X2, prot = X3) %>% 
      mutate(prot = paste0(tar, "_", prot)) %>%
      group_by(prot) %>%
      # assumption - peptide is amyloid if at least one hexagram has prob > 0.5, 
      # so we take maximum probabilities for all hexagrams belonging to the peptide
      summarise(prob = max(prob), tar = unique(tar), len = 5 + length(prot)) %>%
      mutate(len_range = cut(len, include.lowest = TRUE, breaks = c(5, 6, 10, 15, 25)))
    
    perf_measures <- lapply(levels(preds_df[["len_range"]]), function(single_range) {
      dat <- preds_df[preds_df[["len_range"]] == single_range, ]
      tryCatch(HMeasure(dat[["tar"]], dat[["prob"]])[["metrics"]], 
               warning = function(w) 
                 list(HMeasure(dat[["tar"]], dat[["prob"]])[["metrics"]],w))
    })
    
    #check if scores were switched
    reverted <- !sapply(perf_measures, is.data.frame)
    
    aggregated_measures <- data.frame(len_range = levels(preds_df[["len_range"]]), 
                                      n = as.vector(table(preds_df[["len_range"]])),
                                      reverted = reverted,
                                      do.call(rbind, lapply(1L:length(reverted), function(single_range_id) {
                                        if(reverted[single_range_id]) {
                                          perf_measures[[single_range_id]][[1]]
                                        } else {
                                          perf_measures[[single_range_id]]
                                        }
                                      }))
    )
    
    list(aggregated_measures, imp_bigrams)
  })
}



train_forests <- function(pos, pos_id, neg, neg_id) {

}
