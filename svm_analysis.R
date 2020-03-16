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


perf <- lapply(unique(all_cvs[["source_file"]]), function(ith_learner) 
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
    
    lapply(unique(pred_df[["len_group"]]), function(ith_group)
      HMeasure(true.class = filter(pred_df, len_group == ith_group)[["target"]],
               scores = filter(pred_df, len_group == ith_group)[["pred"]])[["metrics"]] %>%
        mutate(fold = ith_fold, source_file = ith_learner, len_group = ith_group)
    ) %>% bind_rows
  }) %>% bind_rows
) %>% bind_rows %>% 
  mutate(len_group = as.character(len_group), 
         len_group = factor(len_group, levels = sort_group(unique(len_group))))

ggplot(perf, aes(x = source_file, y = Spec)) +
  geom_point() + 
  facet_wrap( ~ len_group)
