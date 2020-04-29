library(dplyr)
library(ggplot2)
library(drake)
library(biogram)
library(tidyr)

if(Sys.info()[["nodename"]] %in% c("amyloid", "phobos", "huawei")) {
  data_path <- "/home/michal/Dropbox/AMP-analysis/AmpGram-analysis/"
}
if(Sys.info()[["nodename"]] %in% c("kasia-MACH-WX9", "ryzen")) {
  data_path <- "/home/kasia/Dropbox/AmpGram-analysis/"
}

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

calculate_len_distribution <- function(lens) {
  lens_groups <- mutate(lens, len_group = cut(lens[["len"]], c(min(lens[["len"]]), 10, 19, 26, 36, 60, 710, max(lens[["len"]])),
                                              include.lowest = TRUE))
  lens_groups %>% 
    group_by(len_group) %>% 
    summarise(count = n())
}


plots_AmpGram <- drake_plan(neg = readd(negative_data),
                            pos = readd(cdhit_data),
                            pos_ids = readd(cdhit_data_ids),
                            neg_ids = readd(negative_data_ids),
                            pos_names_by_len_groups = select_names_by_len_group(pos, pos_ids),
                            neg_names_by_len_groups = select_names_by_len_group(neg, neg_ids),
                            groups = names(pos_ids),
                            pos_comp = get_comp_df(groups, pos, pos_names_by_len_groups, "Positive"),
                            neg_comp = get_comp_df(groups, neg, neg_names_by_len_groups, "Negative"),
                            composition_df = bind_rows(pos_comp, neg_comp) %>% 
                              mutate(len_group = factor(len_group, levels = sort_group(unique(len_group))),
                                     dataset = factor(dataset, levels = c("Positive", "Negative"))),
                            composition_plot = ggplot(composition_df, aes(x = aa, y = freq, fill = dataset)) +
                              geom_col(position = "dodge") +
                              facet_wrap(~ len_group, ncol = 2, scales = "free_x") +
                              xlab("Amino acid") +
                              ylab("Frequency") +
                              theme_bw(),
                            UniProt_seqs = read_fasta(paste0(data_path, "data/input-seqs.fasta")),
                            UniProt_lens = lapply(1L:length(UniProt_seqs), function(i) {
                              data.frame(id = names(UniProt_seqs[i]),
                                         len = length(UniProt_seqs[[i]]))
                            }) %>% 
                              bind_rows(),
                            UniProt_len_distribution = calculate_len_distribution(UniProt_lens),
                            dbAMP_data = read.csv("./data/dbamp_df.csv", stringsAsFactors = FALSE),
                            dbAMP_lens = lapply(dbAMP_data[["Sequence"]], function(ith_seq) {
                              strsplit(ith_seq, "")[[1]] %>% 
                                length() 
                            }) %>% 
                              unlist() %>% 
                              data.frame(len = .),
                            dbAMP_len_distribution = calculate_len_distribution(dbAMP_lens),
                            benchmark_summ_table = pivot_longer(readd(benchmark_summ), c(AUC, MCC, precision, 
                                                                                 sensitivity, specificity)) %>% 
                              mutate(AmpGram = Software == "AmpGram_full"),
                            benchmark_summ_plot = ggplot(benchmark_summ_table, aes(x = Software, y = value, color = AmpGram)) +
                              geom_point() +
                              facet_grid(len_group ~ name, scales = "free_y") +
                              theme(axis.text.x = element_text(angle = 90)))

make(plots_AmpGram, seed = 990)

file.copy(from = ".drake", to = paste0(data_path, "drake-cache"), recursive = TRUE, overwrite = TRUE)