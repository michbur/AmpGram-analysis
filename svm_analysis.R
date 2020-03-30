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


if(Sys.info()[["nodename"]] %in% c("amyloid", "phobos", "huawei")) {
  data_path <- "/home/michal/Dropbox/AMP-analysis/AmpGram-analysis/"
}
if(Sys.info()[["nodename"]] %in% c("kasia-MACH-WX9", "ryzen")) {
  data_path <- "/home/kasia/Dropbox/AmpGram-analysis/"
}

all_cvs <- lapply(list.files(paste0(data_path, "results"), 
                             pattern = "csv", full.names = TRUE), function(ith_file) {
                               read.csv(ith_file, stringsAsFactors = FALSE) %>% 
                                 mutate(source_file = get_penultimate(strsplit(ith_file, 
                                                                               split = "[/|]")[[1]]))
                             }) %>% 
  bind_rows() 

count_longest <- function(x) {
  splitted_x <- strsplit(x = paste0(as.numeric(x > 0.5), collapse = ""),
                         split = "0")[[1]]
  len <- unname(sapply(splitted_x, nchar))
  if (length(len[len > 0]) == 0) {
    0 } else {
      len[len > 0]
    }
}


layer_dat <- group_by(all_cvs, source_peptide, target, group, fold, source_file) %>% 
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
    lapply(unique(all_preds[["group"]]), function(ith_group) {
      perf_dat <- filter(all_preds, 
                         group == ith_group, 
                         fold == ith_fold,
                         source_file == ith_learner)
      HMeasure(true.class = perf_dat[["target"]],
               scores = perf_dat[["pred"]])[["metrics"]] %>%
        mutate(fold = ith_fold, source_file = ith_learner, group = ith_group)
    }) %>% bind_rows
  ) %>% bind_rows
) %>% 
  bind_rows %>% 
  mutate(group = factor(group, levels = sort_group(unique(group))))

png(filename = "cv_res.png", width = 780, height = 500)
ggplot(all_perfs, aes(x = source_file, y = AUC)) +
  geom_point() + 
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 4) +
  facet_wrap( ~ group, nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))
dev.off()

