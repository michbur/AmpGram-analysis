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
  scale_x_continuous(limits = c(0, 200))


group_by(all_cvs, source_peptide, target) %>% 
  summarise(fraction_true = mean(pred > 0.5))


all_cvs_pos <- group_by(all_cvs, source_peptide) %>% 
  mutate(fraction_true = mean(mean(pred > 0.5))) %>% 
  ungroup() %>% 
  mutate(mer_pos = as.numeric(sapply(strsplit(mer_id, split = "m"), 
                                     last))) %>% 
  mutate(cfrac_true = cut(fraction_true, breaks = 0L:5/5, 
                          include.lowest = TRUE))


set.seed(1)
twenty_peps <- filter(all_cvs_pos, group == "(36,60]") %>% 
  pull(source_peptide) %>% 
  unique %>% 
  sample(20)

filter(all_cvs_pos, source_peptide %in% twenty_peps) %>% 
  ggplot(aes(x = mer_pos, y = pred,
             group = source_peptide)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.5, color = "red") +
  facet_grid(cfrac_true ~ target)

group_by(all_cvs_pos, group, source_peptide, target) %>% 
  summarise(fraction_true = mean(pred > 0.5),
            len = length(pred) + 9) %>% 
  mutate(cfrac_true = cut(fraction_true, breaks = 0L:10/10, 
                          include.lowest = TRUE)) %>% 
  group_by(group, cfrac_true, target) %>% 
  summarise(n = length(target)) %>% 
  ungroup() %>% 
  mutate(group = factor(group, levels = sort_group(unique(group)))) %>% 
  group_by(group, cfrac_true) %>% 
  mutate(n_prop = n/sum(n)) %>% 
  ggplot(aes(x = group, y = n_prop, fill = target)) +
  geom_col() +
  facet_wrap(~ cfrac_true, nrow = 1) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45))
