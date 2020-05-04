library(dplyr)
library(ggplot2)
library(drake)
library(biogram)
library(tidyr)
library(ranger)
source("./functions/benchmark_functions.R")
source("./functions/test_alphabet.R")
load("./data/benchmark_data.RData")

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
           Dataset = type,
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

get_prot_preds <- function(protein, regions) {
  prot <- read_fasta(paste0("./data/", protein, ".fasta"))
  mer_df <- unlist(unname(prot)) %>% 
    matrix(nrow = 1) %>% 
    get_single_seq_mers() %>% 
    data.frame(stringsAsFactors = FALSE) %>% 
    mutate(source_peptide = protein,
           mer_id = paste0(source_peptide, "_m", 1L:nrow(.)))
  imp_ngrams <- count_imp_ampgrams(mer_df, imp_bigrams)
  amp_regions <- regions
  pred_mers <- mutate(mer_df, 
                      pred = predict(full_model_mers, as.matrix(imp_ngrams))[["predictions"]][,"TRUE"],
                      pos = 1L:nrow(mer_df),
                      region = ifelse(pos %in% amp_regions, "AMP", "non-AMP"))
}



get_detailed_preds <- function(pred_mers) {
  do.call(rbind, lapply(unique(pred_mers[["source_peptide"]]), function(single_prot) {
    mer_preds <- pred_mers[pred_mers[["source_peptide"]] == single_prot, "pred"]
    pos_matrix <- do.call(cbind, get_ngrams_ind(length(mer_preds) + 9, 10, 0))
    data.frame(
      Protein = single_prot,
      Pos = unique(as.vector(pos_matrix)),
      Probability = unlist(lapply(unique(as.vector(pos_matrix)), function(i)
        max(mer_preds[which(pos_matrix == i, arr.ind = TRUE)[, "row"]])
      )))
  }))
}


plots_AmpGram <- drake_plan(neg = readd(negative_data),
                            pos = readd(cdhit_data),
                            pos_ids = readd(cdhit_data_ids),
                            neg_ids = readd(negative_data_ids),
                            pos_names_by_len_groups = select_names_by_len_group(pos, pos_ids),
                            neg_names_by_len_groups = select_names_by_len_group(neg, neg_ids),
                            groups = names(pos_ids),
                            pos_comp = get_comp_df(groups, pos, pos_names_by_len_groups, "AMP"),
                            neg_comp = get_comp_df(groups, neg, neg_names_by_len_groups, "Non-AMP"),
                            composition_df = bind_rows(pos_comp, neg_comp) %>% 
                              mutate(len_group = factor(len_group, levels = sort_group(unique(len_group)))),
                            composition_plot = ggplot(composition_df, aes(x = aa, y = freq, fill = Dataset)) +
                              geom_col(position = "dodge") +
                              facet_wrap(~ len_group, scales = "free_x", nrow = 2, 
                                         labeller = labeller(len_group = as_labeller(c("[11,19]"="11-19 aa", "(19,26]"="20-26 aa", "(26,36]"="27-36 aa", 
                                                                                       "(36,60]"="37-60 aa", "(60,710]"="61-710 aa")))) +
                              xlab("Amino acid") +
                              ylab("Frequency") +
                              theme_bw() +
                              scale_fill_manual(values = c("#f8766d", "#878787")) +
                              theme(legend.position = c(0.85,0.25)),
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
                            benchmark_summ_table = pivot_longer(readd(benchmark_summ), c(AUC, MCC, Precision, 
                                                                                         Sensitivity, Specificity), values_to = "Value"),
                            benchmark_summ_plot = ggplot(benchmark_summ_table, aes(x = Software, y = Value)) +
                              geom_point() +
                              facet_grid(len_group ~ name, scales = "free_y", labeller = labeller(len_group = as_labeller(
                                c("all"="all lengths", "[11,19]"="11-19 aa", "(19,26]"="20-26 aa", "(26,36]"="27-36 aa", 
                                  "(36,60]"="37-60 aa", "(60,710]"="61-710 aa")))) +
                              theme_bw() +
                              theme(axis.text.x = element_text(angle = 90)),
                            lactoferrin_preds = get_prot_preds("bovine_lactoferrin", c(20:30, 36:60, 287:303)),
                            lactoferrin_detailed_preds = get_detailed_preds(lactoferrin_preds),
                            lactoferrin_profile_plot = ggplot(lactoferrin_detailed_preds, aes(x = Pos, y = Probability, group = Protein)) +
                              geom_ribbon(mapping = aes(xmin = 20, xmax = 30, ymin = 0, ymax = 1), fill = "red") +
                              geom_ribbon(mapping = aes(xmin = 36, xmax = 60, ymin = 0, ymax = 1), fill = "red") +
                              geom_ribbon(mapping = aes(xmin = 287, xmax = 303, ymin = 0, ymax = 1), fill = "red") +
                              geom_point() +
                              geom_line() +
                              geom_hline(yintercept = 0.5, color = "red") +
                              theme_bw(),
                            lactoferrin_mers_plot = ggplot(lactoferrin_preds, aes(x = pos, y = pred)) +
                              geom_hline(yintercept = 0.5, color = "red") +
                              geom_rect(mapping = aes(xmin = 20, xmax = 30, ymin = 0, ymax = 1), fill = "#f8766d") +
                              geom_rect(mapping = aes(xmin = 36, xmax = 60, ymin = 0, ymax = 1), fill = "#f8766d") +
                              geom_rect(mapping = aes(xmin = 287, xmax = 303, ymin = 0, ymax = 1), fill = "#f8766d") +
                              geom_segment(x = 1:length(lactoferrin_preds[["pos"]]), y = lactoferrin_preds[["pred"]],
                                           xend = 10:(length(lactoferrin_preds[["pos"]])+9), yend = lactoferrin_preds[["pred"]],
                                           colour = ifelse(lactoferrin_preds[["pred"]] < 0.5, "#878787", "black")) +
                              xlab("Position") +
                              ylab("Prediction") +
                              ylim(c(0,1)) +
                              scale_x_continuous(breaks = seq(0, length(lactoferrin_preds[["pos"]])+9, by = 20), limits = c(0,708), expand = c(0.01,0.01))+
                              theme_bw(),
                            thrombin_preds = get_prot_preds("thrombin", c(597:622)),
                            thrombin_detailed_preds = get_detailed_preds(thrombin_preds),
                            thrombin_mers_plot = ggplot(thrombin_preds, aes(x = pos, y = pred)) +
                              geom_hline(yintercept = 0.5, color = "red") +
                              geom_rect(mapping = aes(xmin = 597, xmax = 622, ymin = 0, ymax = 1), fill = "#f8766d") +
                              geom_segment(x = 1:length(thrombin_preds[["pos"]]), y = thrombin_preds[["pred"]],
                                           xend = 10:(length(thrombin_preds[["pos"]])+9), yend = thrombin_preds[["pred"]],
                                           colour = ifelse(thrombin_preds[["pred"]] < 0.5, "#878787", "black")) +
                              xlab("Position") +
                              ylab("Prediction") +
                              ylim(c(0,1)) +
                              scale_x_continuous(breaks = seq(0, length(thrombin_preds[["pos"]])+9, by = 20), limits = c(0,622), expand = c(0.01,0.01)) +
                              theme_bw(),
                            Nobles_benchmark_plot = readd(Nobles_datasets_benchmark_res)[,c(1:4, 9:11)] %>% 
                              mutate(AUC = as.double(AUC)) %>% 
                              pivot_longer(c("AUC", "MCC", "Precision", "Sensitivity", "Specificity"), names_to = "Measure", values_to = "Value") %>% 
                              ggplot(aes(x = Software, y = Value)) +
                              geom_point() +
                              facet_grid(Dataset ~ Measure, scales = "free_y") + 
                              theme_bw() +
                              theme(axis.text.x = element_text(angle = 90))
)


make(plots_AmpGram, seed = 990)

file.copy(from = ".drake", to = paste0(data_path, "drake-cache"), recursive = TRUE, overwrite = TRUE)

cairo_ps(filename = "benchmark.eps", width = 10)
readd(benchmark_summ_plot)
dev.off()
  
cairo_ps(filename = "aa_comp.eps", width = 10, height = 4)
readd(composition_plot)
dev.off()

cairo_ps(filename = "Noble_benchmark.eps", width = 10, height = 4)
readd(Nobles_benchmark_plot)
dev.off()

cairo_ps(filename = "lactoferrin_mers.eps", width = 10, height = 2.5)
readd(lactoferrin_mers_plot)
dev.off()

cairo_ps(filename = "thrombin_mers.eps", width = 10, height = 2.5)
readd(thrombin_mers_plot)
dev.off()
