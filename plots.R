library(dplyr)
library(ggplot2)
library(drake)
library(biogram)

select_names_by_len_group <- function(dataset, dataset_ids) {
  groups <- names(dataset_ids)
  lapply(groups, function(ith_group) {
    data.frame(names = unlist(unname(dataset_ids[[ith_group]]))) %>% 
                 mutate(len_group = ith_group)
  }) %>% bind_rows()
}

get_comp_df <- function(groups, dataset, dataset_names_by_len_groups, type) {
  lapply(groups, function(ith_group) {
    names <- filter(dataset_names_by_len_groups, len_group == ith_group)[["names"]]
    aac <- data.frame(table(unlist(dataset[names]))) %>% 
      setNames(c("aa", "freq")) 
    mutate(aac, 
           freq = freq/sum(freq),
           dataset = type,
           len_group = ith_group)
  }) %>% bind_rows()
}

sort_group <- function(x) {
  splitted_x <- sapply(strsplit(x, split = ","), function(i) i[1])
  x_order <- order(as.numeric(gsub(pattern = "[^0-9]", 
                                   replacement = "", x = splitted_x)),
                   na.last = FALSE)
  x[x_order]
}


neg <- readd(negative_data)
pos <- readd(cdhit_data)
pos_ids <- readd(cdhit_data_ids)
neg_ids <- readd(negative_data_ids)

pos_names_by_len_groups <- select_names_by_len_group(pos, pos_ids)
neg_names_by_len_groups <- select_names_by_len_group(neg, neg_ids)

groups <- names(pos_ids)

pos_comp <- get_comp_df(groups, pos, pos_names_by_len_groups, "Positive")
neg_comp <- get_comp_df(groups, neg, neg_names_by_len_groups, "Negative")

composition_df <- bind_rows(pos_comp, neg_comp) %>% 
  mutate(len_group = factor(len_group, levels = sort_group(unique(len_group))),
         dataset = factor(dataset, levels = c("Positive", "Negative")))

ggplot(composition_df, aes(x = aa, y = freq, fill = dataset)) +
  geom_col(position = "dodge") +
  #facet_grid(dataset ~ len_group) +
  facet_wrap(~ len_group, ncol = 2, scales = "free_x") +
  xlab("Amino acid") +
  ylab("Frequency") +
  theme_bw()


### UniProt length distribution
if(Sys.info()[["nodename"]] %in% c("amyloid", "phobos", "huawei")) {
  data_path <- "/home/michal/Dropbox/AMP-analysis/AmpGram-analysis/"
}
if(Sys.info()[["nodename"]] %in% c("kasia-MACH-WX9", "ryzen")) {
  data_path <- "/home/kasia/Dropbox/AmpGram-analysis/"
}

uniprot_seqs <- read_fasta(paste0(data_path, "data/input-seqs.fasta"))
length(uniprot_seqs)

lens <- lapply(1L:length(uniprot_seqs), function(i) {
  data.frame(id = names(uniprot_seqs[i]),
             len = length(uniprot_seqs[[i]]))
}) %>% 
  bind_rows() 

lens_groups <- mutate(lens, len_group = cut(lens[["len"]], c(min(lens[["len"]]), 10, 19, 26, 36, 60, 710, max(lens[["len"]])),
                         include.lowest = TRUE))

lens_groups %>% 
  group_by(len_group) %>% 
  summarise(count = n())

ggplot(lens, aes(x = len)) +
  geom_density() +
  xlim(c(0, 1000)) +
  geom_vline(xintercept = c(11, 19, 26, 36, 60, 710), color = "red") +
  geom_vline(xintercept = 1000, linetype = "dashed") +
  theme_bw()
