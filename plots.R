library(protr)
library(dplyr)
library(ggplot2)
library(drake)

neg <- readd(negative_data)
pos <- readd(cdhit_data)
pos_ids <- readd(cdhit_data_ids)
neg_ids <- readd(negative_data_ids)


select_names_by_len_group <- function(dataset, dataset_ids) {
  groups <- names(dataset_ids)
  lapply(groups, function(ith_group) {
    data.frame(names = unlist(unname(dataset_ids[[ith_group]]))) %>% 
                 mutate(len_group = ith_group)
  }) %>% bind_rows()
}

pos_names_by_len_groups <- select_names_by_len_group(pos, pos_ids)
neg_names_by_len_groups <- select_names_by_len_group(neg, neg_ids)

groups <- names(pos_ids)

merge_seqs <- function(groups, dataset, dataset_names_by_len_group) {
  lapply(groups, function(ith_group) {
  names <- filter(dataset_names_by_len_group, len_group == ith_group)[["names"]]
  lapply(names, function(ith_name) {
    unlist(paste(dataset[[ith_name]], collapse = ""))
  }) %>%  
  unlist() %>% 
  paste0(collapse = "")
}) %>% unlist()
}

get_comp_df <- function(groups, dataset, dataset_names_by_len_groups, type) {
  lapply(merge_seqs(groups, dataset, dataset_names_by_len_groups), function(ith_merged) {
  aac <- extractAAC(ith_merged)
  data.frame('aa' = names(aac),
             'freq' = aac,
             'dataset' = type)
}) %>% 
    bind_rows() %>% 
    mutate(len_group = c(sapply(groups, function(i) rep(i, 20))))
}

pos_comp <- get_comp_df(groups, pos, pos_names_by_len_groups, "Positive")
neg_comp <- get_comp_df(groups, neg, neg_names_by_len_groups, "Negative")

composition_df <- bind_rows(pos_comp, neg_comp) %>% 
  mutate(len_group = factor(len_group, levels = sort_group(unique(len_group))))

ggplot(composition_df, aes(x = aa, y = freq)) +
  geom_col() +
  facet_grid(dataset ~ len_group) +
#  facet_wrap(~dataset)
  xlab("Amino acid") +
  ylab("Frequency") +
  theme_bw()




### UniProt length distribution

uniprot_seqs <- read_fasta("/home/kasia/Dropbox/AmpGram-analysis/data/input-seqs.fasta")
length(uniprot_seqs)

lens <- lapply(1L:length(uniprot_seqs), function(i) {
  data.frame(id = names(uniprot_seqs[i]),
             len = length(uniprot_seqs[[i]]))
}) %>% bind_rows()

ggplot(lens, aes(x = len)) +
  geom_density() +
  xlim(c(0, 1000)) +
  geom_vline(xintercept = c(11, 19, 26, 36, 60, 710), color = "red") +
  geom_vline(xintercept = 1000, linetype = "dashed") +
  theme_bw()
