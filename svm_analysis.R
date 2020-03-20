library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(ranger)
library(hmeasure)

sort_group <- function(x) {
  splitted_x <- sapply(strsplit(x, split = ","), function(i) i[1])
  x_order <- order(as.numeric(gsub(pattern = "[^0-9]", 
                                   replacement = "", x = splitted_x)))
  x[x_order]
}

get_penultimate <- function(x)
  x[length(x) - 1]

all_cvs <- lapply(list.files("/home/michal/Dropbox/AMP-analysis/AmpGram-analysis/results", 
                             pattern = "csv", full.names = TRUE), function(ith_file) {
                               read.csv(ith_file, stringsAsFactors = FALSE) %>% 
                                 mutate(source_file = get_penultimate(strsplit(ith_file, 
                                                                               split = "[/|]")[[1]]))
                             }) %>% 
  bind_rows() 

# y <- runif(100)
# 
# count_longest <- function(x) {
#   splitted_x <- strsplit(x = paste0(as.numeric(x > 0.5), collapse = ""), 
#                          split = "0")[[1]]
#   len <- unname(sapply(splitted_x, nchar))
#   len[len > 0]
# }
# 
# max(count_longest(y))
# sum(count_longest(y) > 4)


layer_dat <- group_by(all_cvs, source_peptide, target, group, fold, source_file) %>% 
  summarise(fraction_true = mean(pred > 0.5),
            pred_mean = mean(pred),
            pred_median = median(pred),
            n_peptide = length(pred),
            n_pos = sum(pred > 0.5),
            pred_min = min(pred),
            pred_max = max(pred)) %>% 
  ungroup() %>% 
  mutate(target = factor(target))

#filter(all_cvs, source_peptide == "AMP748", source_file == ith_learner)


all_preds <- lapply(unique(all_cvs[["source_file"]]), function(ith_learner) 
  lapply(1L:5, function(ith_fold) {
    ranger_train_data <- filter(layer_dat, 
                                fold != ith_fold,
                                source_file == ith_learner)[, c("target", "fraction_true",
                                                                "pred_mean", "pred_median",
                                                                "n_peptide", "n_pos", "pred_min",
                                                                "pred_max")]
    
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
    
    # lapply(unique(pred_df[["len_group"]]), function(ith_group)
    #   HMeasure(true.class = filter(pred_df, len_group == ith_group)[["target"]],
    #            scores = filter(pred_df, len_group == ith_group)[["pred"]])[["metrics"]] %>%
    #     mutate(fold = ith_fold, source_file = ith_learner, len_group = ith_group)
    # ) %>% bind_rows
    
    pred_df
  }) %>% bind_rows
) %>% bind_rows %>% 
  mutate(len_group = as.character(len_group), 
         len_group = factor(len_group, levels = sort_group(unique(len_group))))


all_perfs <- lapply(unique(all_preds[["source_file"]]), function(ith_learner) 
  lapply(1L:5, function(ith_fold) 
    lapply(unique(all_preds[["len_group"]]), function(ith_group) {
      perf_dat <- filter(all_preds, 
                         len_group == ith_group, 
                         fold == ith_fold,
                         source_file == ith_learner)
      HMeasure(true.class = perf_dat[["target"]],
               scores = perf_dat[["pred"]])[["metrics"]] %>%
        mutate(fold = ith_fold, source_file = ith_learner, len_group = ith_group)
    }) %>% bind_rows
  ) %>% bind_rows
) %>% bind_rows

ggplot(all_perfs, aes(x = source_file, y = AUC)) +
  geom_point() + 
  facet_wrap( ~ len_group, nrow = 1)

ggplot(all_preds, aes(x = len_group, y = pred_median, fill = target)) +
  geom_violin() 

ggplot(all_preds, aes(x = pred_mean, y = fraction_true, color = target)) +
  stat_density2d(aes(alpha = ..level.., fill = target), geom = "polygon",
                 color = "black") + 
  facet_wrap(~ len_group)

filter(all_preds, n_peptide > 100, source_file == "(19,26]") %>% 
  ggplot(aes(x = pred_mean, y = fraction_true, color = target)) +
  stat_density2d(aes(alpha = ..level.., fill = target), geom = "polygon",
                 color = "black") + 
  facet_wrap(~ len_group) 
