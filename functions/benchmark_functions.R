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
