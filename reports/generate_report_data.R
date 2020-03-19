library(dplyr)
library(ranger)
library(kernlab)
library(biogram)
library(ggplot2)
library(hmeasure)
library(reshape2)
library(ggbeeswarm)

if(Sys.info()[["nodename"]] %in% c("amyloid", "phobos", "huawei")) {
  data_path <- "/home/michal/Dropbox/AMP-analysis/AmpGram-analysis/"
}
if(Sys.info()[["nodename"]] %in% c("kasia-MACH-WX9", "ryzen")) {
  data_path <- "/home/kasia/Dropbox/AmpGram-analysis/"
}

sort_group <- function(x) {
  splitted_x <- sapply(strsplit(x, split = ","), function(i) i[1])
  x_order <- order(as.numeric(gsub(pattern = "[^0-9]", 
                                   replacement = "", x = splitted_x)))
  x[x_order]
}

get_penultimate <- function(x) {
  x[length(x) - 1]
}

all_cvs <- lapply(list.files(paste0(data_path, "results"), 
                             pattern = "csv", full.names = TRUE), function(ith_file) {
                               read.csv(ith_file, stringsAsFactors = FALSE) %>% 
                                 mutate(source_file = get_penultimate(strsplit(ith_file, 
                                                                               split = "[/|]")[[1]]))
                             }) %>% 
  bind_rows() 

all_cvs_pos <- group_by(all_cvs, source_peptide) %>% 
  mutate(fraction_true = mean(pred > 0.5)) %>% 
  ungroup() %>% 
  mutate(mer_pos = as.numeric(sapply(strsplit(mer_id, split = "m"), 
                                     last))) %>% 
  mutate(cfrac_true = cut(fraction_true, breaks = 0L:5/5, 
                          include.lowest = TRUE))

pred_len <- group_by(all_cvs, source_peptide, target) %>% 
  summarise(fraction_true = mean(pred > 0.5),
            len = length(pred) + 9) %>% 
  mutate(cfrac_true = cut(fraction_true, breaks = 0L:5/5, 
                          include.lowest = TRUE)) %>% 
  group_by(target, len, cfrac_true) %>% 
  summarise(n = length(fraction_true)) %>% 
  group_by(target, len) %>% 
  mutate(n_prop = n/sum(n))


layer_dat <- group_by(all_cvs, source_peptide, target, group, fold, source_file) %>% 
  summarise(fraction_true = mean(pred > 0.5),
            pred_mean = mean(pred),
            pred_median = median(pred),
            n_peptide = length(pred),
            n_pos = sum(pred > 0.5),
            pred_min = min(pred),
            pred_max = max(pred),
            n_pred_0.9 = sum(pred > 0.9),
            n_pred_0.8 = sum(pred > 0.8)) %>% 
  ungroup() %>% 
  mutate(target = factor(target))

perf_rf <- lapply(unique(all_cvs[["source_file"]]), function(ith_learner) 
  lapply(1L:5, function(ith_fold) {
    ranger_train_data <- filter(layer_dat, 
                                fold != ith_fold,
                                source_file == ith_learner)[, c("target", "fraction_true",
                                                                "pred_mean", "pred_median",
                                                                "n_peptide", "n_pos", "pred_min",
                                                                "pred_max", "n_pred_0.9", "n_pred_0.8")]
    
    ranger_test_data <- filter(layer_dat, 
                               fold == ith_fold,
                               source_file == ith_learner)
    
    model_cv <- ranger(dependent.variable.name = "target", data =  ranger_train_data, 
                       write.forest = TRUE, probability = TRUE, num.trees = 500, 
                       verbose = FALSE, classification = TRUE)
    
    pred_df <- mutate(ranger_test_data,
                      source_peptide = ranger_test_data[["source_peptide"]],
                      pred = predict(model_cv, ranger_test_data)[["predictions"]][, "TRUE"]) %>% 
      mutate(len_group = cut(n_peptide + 9, breaks = c(0, 20, 50, 100, 200, 800),
                             include.lowest = TRUE))
    
    lapply(unique(pred_df[["len_group"]]), function(ith_group)
      HMeasure(true.class = filter(pred_df, len_group == ith_group)[["target"]],
               scores = filter(pred_df, len_group == ith_group)[["pred"]])[["metrics"]] %>%
        mutate(fold = ith_fold, source_file = ith_learner, len_group = ith_group)
    ) %>% bind_rows
  }) %>% bind_rows
) %>% bind_rows %>% 
  mutate(len_group = as.character(len_group), 
         len_group = factor(len_group, levels = sort_group(unique(len_group))))



perf_vanilla <- lapply(unique(all_cvs[["source_file"]]), function(ith_learner) {
  lapply(1L:5, function(ith_fold) {
    layer_train_data <- filter(layer_dat, 
                               fold != ith_fold,
                               source_file == ith_learner)[, c("target", "fraction_true",
                                                               "pred_mean", "pred_median",
                                                               "n_peptide", "n_pos", "pred_min",
                                                               "pred_max", "n_pred_0.9", "n_pred_0.8")]
    
    layer_test_data <- filter(layer_dat, 
                              fold == ith_fold,
                              source_file == ith_learner)
    
    vanilla_model <- ksvm(x = as.matrix(layer_train_data[,2:10]), 
                          y = as.matrix(layer_train_data[["target"]]), 
                          kernel = "vanilladot",
                          prob.model = TRUE)
    
    pred_df <- mutate(layer_test_data,
                      source_peptide = layer_test_data[["source_peptide"]],
                      pred = predict(vanilla_model, 
                                     as.matrix(layer_test_data[,6:14]), 
                                     type = "probabilities")[, "TRUE"]) %>% 
      mutate(len_group = cut(n_peptide + 9, breaks = c(0, 20, 50, 100, 200, 800),
                             include.lowest = TRUE))
    
    lapply(unique(pred_df[["len_group"]]), function(ith_group) {
      HMeasure(true.class = filter(pred_df, len_group == ith_group)[["target"]],
               scores = filter(pred_df, len_group == ith_group)[["pred"]])[["metrics"]] %>%
        mutate(fold = ith_fold, source_file = ith_learner, len_group = ith_group)
    }) %>% bind_rows
  }) %>% bind_rows
}) %>% bind_rows %>% 
  mutate(len_group = as.character(len_group), 
         len_group = factor(len_group, levels = sort_group(unique(len_group))))


perf_poly <- lapply(unique(all_cvs[["source_file"]]), function(ith_learner) {
  lapply(1L:5, function(ith_fold) {
    layer_train_data <- filter(layer_dat, 
                               fold != ith_fold,
                               source_file == ith_learner)[, c("target", "fraction_true",
                                                               "pred_mean", "pred_median",
                                                               "n_peptide", "n_pos", "pred_min",
                                                               "pred_max", "n_pred_0.9", "n_pred_0.8")]
    
    layer_test_data <- filter(layer_dat, 
                              fold == ith_fold,
                              source_file == ith_learner)
    
    poly_model <- ksvm(x = as.matrix(layer_train_data[,2:10]), 
                       y = as.matrix(layer_train_data[["target"]]), 
                       kernel = "polydot",
                       prob.model = TRUE)
    
    pred_df <- mutate(layer_test_data,
                      source_peptide = layer_test_data[["source_peptide"]],
                      pred = predict(poly_model, 
                                     as.matrix(layer_test_data[,6:14]), 
                                     type = "probabilities")[, "TRUE"]) %>% 
      mutate(len_group = cut(n_peptide + 9, breaks = c(0, 20, 50, 100, 200, 800),
                             include.lowest = TRUE))
    
    lapply(unique(pred_df[["len_group"]]), function(ith_group) {
      HMeasure(true.class = filter(pred_df, len_group == ith_group)[["target"]],
               scores = filter(pred_df, len_group == ith_group)[["pred"]])[["metrics"]] %>%
        mutate(fold = ith_fold, source_file = ith_learner, len_group = ith_group)
    }) %>% bind_rows
  }) %>% bind_rows
}) %>% bind_rows %>% 
  mutate(len_group = as.character(len_group), 
         len_group = factor(len_group, levels = sort_group(unique(len_group))))


perf_rbf <- lapply(unique(all_cvs[["source_file"]]), function(ith_learner) {
  lapply(1L:5, function(ith_fold) {
    layer_train_data <- filter(layer_dat, 
                               fold != ith_fold,
                               source_file == ith_learner)[, c("target", "fraction_true",
                                                               "pred_mean", "pred_median",
                                                               "n_peptide", "n_pos", "pred_min",
                                                               "pred_max", "n_pred_0.9", "n_pred_0.8")]
    
    layer_test_data <- filter(layer_dat, 
                              fold == ith_fold,
                              source_file == ith_learner)
    
    rbf_model <- ksvm(x = as.matrix(layer_train_data[,2:10]), 
                      y = as.matrix(layer_train_data[["target"]]), 
                      kernel = "rbfdot",
                      prob.model = TRUE)
    
    pred_df <- mutate(layer_test_data,
                      source_peptide = layer_test_data[["source_peptide"]],
                      pred = predict(rbf_model, 
                                     as.matrix(layer_test_data[,6:14]), 
                                     type = "probabilities")[, "TRUE"]) %>% 
      mutate(len_group = cut(n_peptide + 9, breaks = c(0, 20, 50, 100, 200, 800),
                             include.lowest = TRUE))
    
    lapply(unique(pred_df[["len_group"]]), function(ith_group) {
      HMeasure(true.class = filter(pred_df, len_group == ith_group)[["target"]],
               scores = filter(pred_df, len_group == ith_group)[["pred"]])[["metrics"]] %>%
        mutate(fold = ith_fold, source_file = ith_learner, len_group = ith_group)
    }) %>% bind_rows
  }) %>% bind_rows
}) %>% bind_rows %>% 
  mutate(len_group = as.character(len_group), 
         len_group = factor(len_group, levels = sort_group(unique(len_group))))

all_perf <- lapply(c("perf_rf", "perf_vanilla", "perf_poly", "perf_rbf"), function(i) {
  get(i) %>% 
    mutate(model = switch(i, "perf_rf" = "rf", "perf_vanilla" = "vanilladot",
                          "perf_poly" = "polydot", "perf_rbf" = "rbfdot")) }) %>% bind_rows()

save.image(file = "report_data.RData")
