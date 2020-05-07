library(dplyr)
library(ggplot2)
library(drake)
library(biogram)
library(tidyr)
library(ranger)
library(xtable)
library(patchwork)
source("./functions/benchmark_functions.R")
source("./functions/test_alphabet.R")
source("./functions/plot_functions.R")
load("./data/benchmark_data.RData")

if(Sys.info()[["nodename"]] %in% c("amyloid", "phobos", "huawei")) {
  data_path <- "/home/michal/Dropbox/AMP-analysis/AmpGram-analysis/"
}
if(Sys.info()[["nodename"]] %in% c("kasia-MACH-WX9", "ryzen")) {
  data_path <- "/home/kasia/Dropbox/AmpGram-analysis/"
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
                            composition_plot = get_composition_plot(composition_df),
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
                            benchmark_summ_plot = get_benchmark_summ_plot(filter(benchmark_summ_table, 
                                                                                 name != "MCC" & len_group %in% c("all", "(60,710]"))),
                            lactoferrin_preds = get_prot_preds("bovine_lactoferrin", c(20:30, 36:60, 287:303)),
                            lactoferrin_detailed_preds = get_detailed_preds(lactoferrin_preds),
                            lactoferrin_profile_plot = get_lactoferrin_profile_plot(lactoferrin_detailed_preds),
                            lactoferrin_mers_plot = get_lactoferrin_mers_plot(lactoferrin_preds),
                            thrombin_preds = get_prot_preds("thrombin", c(527:622)),
                            thrombin_detailed_preds = get_detailed_preds(thrombin_preds),
                            thrombin_mers_plot = get_thrombin_mers_plot(thrombin_preds),
                            Nobles_benchmark_plot = get_Nobles_benchmark_plot(readd(Nobles_datasets_benchmark_res)),
                            
)


make(plots_AmpGram, seed = 990)

file.copy(from = ".drake", to = paste0(data_path, "drake-cache"), recursive = TRUE, overwrite = TRUE)

cairo_ps(filename = paste0(data_path, "publication-results/benchmark.eps"), width = 10, height = 4.5)
readd(benchmark_summ_plot)
dev.off()

cairo_ps(filename = paste0(data_path, "publication-results/aa_comp.eps"), width = 10, height = 4)
readd(composition_plot)
dev.off()

cairo_ps(filename = paste0(data_path, "publication-results/Noble_benchmark.eps"), width = 10, height = 4)
readd(Nobles_benchmark_plot)
dev.off()

cairo_ps(filename = paste0(data_path, "publication-results/lactoferrin_mers.eps"), width = 10, height = 2.5)
readd(lactoferrin_mers_plot)
dev.off()

cairo_ps(filename = paste0(data_path, "publication-results/thrombin_mers.eps"), width = 10, height = 2.5)
readd(thrombin_mers_plot)
dev.off()

cairo_ps(filename = paste0(data_path, "publication-results/proteins_mers.eps"), width = 10, height = 5)
lactoferrin_mers_plot / thrombin_mers_plot
dev.off()
