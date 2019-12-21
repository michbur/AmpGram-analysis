write_benchmark <- function(pos, pos_id, neg, neg_id) {
  seq_list <- c(pos[unlist(lapply(pos_id, function(ith_len_group) ith_len_group[["benchmark"]]))],
    neg[unlist(lapply(neg_id, function(ith_len_group) ith_len_group[["benchmark"]]))])
  write_fasta(seq_list, file = "results/benchmark.fasta")
}

# write_benchmark(pos = readd(cdhit_data),
#                 pos_id = readd(cdhit_data_ids),
#                 neg = readd(negative_data),
#                 neg_id = readd(negative_data_ids))
