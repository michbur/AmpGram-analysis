library(dplyr)
library(ggplot2)
library(ggbeeswarm)

all_cvs <- lapply(list.files("/home/michal/Dropbox/AMP-analysis/AmpGram-analysis/results/", 
           pattern = "csv", full.names = TRUE), read.csv, stringsAsFactors = FALSE) %>% 
  bind_rows() 

pred_len <- group_by(all_cvs, source_peptide, target) %>% 
  summarise(fraction_true = mean(pred > 0.5),
            len = length(pred) + 9) %>% 
  mutate(cfrac_true = cut(fraction_true, breaks = 0L:5/5, 
                          include.lowest = TRUE)) %>% 
  group_by(target, len, cfrac_true) %>% 
  summarise(n = length(fraction_true)) %>% 
  group_by(target, len) %>% 
  mutate(n_prop = n/sum(n))
   
ggplot(pred_len, aes(x = len, y = n, fill = cfrac_true)) +
  geom_col() +
  facet_wrap(~ target, ncol = 1) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 100))

ggplot(pred_len, aes(x = len, y = n_prop, fill = cfrac_true)) +
  geom_col() +
  facet_wrap(~ target, ncol = 1) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 100))
