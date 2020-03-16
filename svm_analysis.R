library(dplyr)
library(ggplot2)
library(ggbeeswarm)

sort_group <- function(x) {
  splitted_x <- sapply(strsplit(x, split = ","), function(i) i[1])
  x_order <- order(as.numeric(gsub(pattern = "[^0-9]", 
                                   replacement = "", x = splitted_x)))
  x[x_order]
}

all_cvs <- lapply(list.files("/home/michal/Dropbox/AMP-analysis/AmpGram-analysis/results/", 
                             pattern = "csv", full.names = TRUE), function(ith_file) {
                               read.csv(ith_file, stringsAsFactors = FALSE) %>% 
                                 mutate(source_file = strsplit(ith_file, split = "//")[[1]][2])
                             }) %>% 
  bind_rows() 

ith_learner <- unique(all_cvs[["source_file"]])[1]

filter(all_cvs, source_file == ith_learner) %>% 
  group_by(source_peptide, target, group, fold) %>% 
  summarise(fraction_true = mean(pred > 0.5),
            pred_mean = mean(pred),
            pred_median = median(pred),
            n_peptide = length(pred),
            pred_min = min(pred),
            pred_max = max(pred))

