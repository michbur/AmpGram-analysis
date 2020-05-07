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

get_composition_plot <- function(composition_df) {
  ggplot(composition_df, aes(x = aa, y = freq, fill = Dataset)) +
    geom_col(position = "dodge") +
    facet_wrap(~ len_group, scales = "free_x", nrow = 2, 
               labeller = labeller(len_group = as_labeller(c("[11,19]"="11-19 aa", "(19,26]"="20-26 aa", "(26,36]"="27-36 aa", 
                                                             "(36,60]"="37-60 aa", "(60,710]"="61-710 aa")))) +
    xlab("Amino acid") +
    ylab("Frequency") +
    theme_bw() +
    scale_fill_manual(values = c("#f8766d", "#878787")) +
    theme(legend.position = c(0.85,0.25))
}


get_benchmark_summ_plot <- function(benchmark_summ_table) {
  ggplot(benchmark_summ_table, aes(x = Software, y = Value)) +
    geom_point() +
    facet_grid(len_group ~ name, scales = "free_y", labeller = labeller(len_group = as_labeller(
      c("all"="all lengths", "[11,19]"="11-19 aa", "(19,26]"="20-26 aa", "(26,36]"="27-36 aa", 
        "(36,60]"="37-60 aa", "(60,710]"="61-710 aa")))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
}

get_lactoferrin_mers_plot <- function(lactoferrin_preds) {
  ggplot(lactoferrin_preds, aes(x = pos, y = pred)) +
    geom_rect(mapping = aes(xmin = 20, xmax = 30, ymin = 0, ymax = 1), fill = "#f8766d") +
    geom_rect(mapping = aes(xmin = 36, xmax = 60, ymin = 0, ymax = 1), fill = "#f8766d") +
    geom_rect(mapping = aes(xmin = 287, xmax = 303, ymin = 0, ymax = 1), fill = "#f8766d") +
    geom_segment(x = 1:length(lactoferrin_preds[["pos"]]), y = lactoferrin_preds[["pred"]],
                 xend = 10:(length(lactoferrin_preds[["pos"]])+9), yend = lactoferrin_preds[["pred"]],
                 colour = ifelse(lactoferrin_preds[["pred"]] < 0.5, "#878787", "black")) +
    geom_hline(yintercept = 0.5, color = "red") +
    xlab("Position") +
    ylab("Prediction") +
    ylim(c(0,1)) +
    labs(tag = "A") +
    ggtitle("Lactoferrin") +
    scale_x_continuous(breaks = seq(0, length(lactoferrin_preds[["pos"]])+9, by = 20), limits = c(0,708), expand = c(0.01,0.01)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}

get_lactoferrin_profile_plot <- function(lactoferrin_detailed_preds) {
  ggplot(lactoferrin_detailed_preds, aes(x = Pos, y = Probability, group = Protein)) +
    geom_ribbon(mapping = aes(xmin = 20, xmax = 30, ymin = 0, ymax = 1), fill = "red") +
    geom_ribbon(mapping = aes(xmin = 36, xmax = 60, ymin = 0, ymax = 1), fill = "red") +
    geom_ribbon(mapping = aes(xmin = 287, xmax = 303, ymin = 0, ymax = 1), fill = "red") +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.5, color = "red") +
    theme_bw()
}


get_thrombin_mers_plot <- function(thrombin_preds) {
  ggplot(thrombin_preds, aes(x = pos, y = pred)) +
    geom_rect(mapping = aes(xmin = 527, xmax = 622, ymin = 0, ymax = 1), fill = "#f8766d") +
    geom_segment(x = 1:length(thrombin_preds[["pos"]]), y = thrombin_preds[["pred"]],
                 xend = 10:(length(thrombin_preds[["pos"]])+9), yend = thrombin_preds[["pred"]],
                 colour = ifelse(thrombin_preds[["pred"]] < 0.5, "#878787", "black")) +
    geom_hline(yintercept = 0.5, color = "red") +
    xlab("Position") +
    ylab("Prediction") +
    ylim(c(0,1)) +
    labs(tag = "B") +
    ggtitle("Thrombin") +
    scale_x_continuous(breaks = seq(0, length(thrombin_preds[["pos"]])+9, by = 20), limits = c(0,622), expand = c(0.01,0.01)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}


get_Nobles_benchmark_plot <- function(Nobles_datasets_benchmark_res) {
  Nobles_datasets_benchmark_res[,c(1:3, 9:11)] %>% 
    mutate(AUC = as.double(AUC)) %>% 
    pivot_longer(c("AUC", "Precision", "Sensitivity", "Specificity"), names_to = "Measure", values_to = "Value") %>% 
    ggplot(aes(x = Software, y = Value)) +
    geom_point() +
    facet_grid(Dataset ~ Measure, scales = "free_y") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
}