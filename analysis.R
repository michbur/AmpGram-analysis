library(dplyr)
library(magrittr)
library(drake)
library(biogram)
library(pbapply)
library(ranger)
library(cvTools)
library(visNetwork)
library(hmeasure)

if(Sys.info()[["nodename"]] %in% c("amyloid", "phobos", "huawei")) {
  data_path <- "/home/michal/Dropbox/AMP-analysis/AmpGram-analysis/"
}
if(Sys.info()[["nodename"]] == "kasia-MACH-WX9") {
  data_path <- "/home/kasia/Dropbox/AmpGram-analysis/"
}

source("./functions/raw_data.R")
source("./functions/cdhit_data.R")
source("./functions/nonstandard_AMPs.R")
source("./functions/cutting_seqs.R")
source("./functions/holdouts.R")
source("./functions/writing_benchmarks.R")
source("./functions/get_mers.R")
source("./functions/count_ampgrams.R")

analysis_AmpGram <- drake_plan(raw_data = read_raw_data(),
                               nonstandard_AMPs = analyze_nonstandard_AMPs(raw_data),
                               cdhit_data = filter_cdhit(raw_data),
                               negative_data = read_and_cut(data_path, lengths(cdhit_data)),
                               cdhit_data_ids = generate_holdout_groups(cdhit_data),
                               negative_data_ids = generate_holdout_groups(negative_data),
                               benchmark_file = write_benchmark(pos = cdhit_data,
                                                                pos_id = cdhit_data_ids,
                                                                neg = negative_data,
                                                                neg_id = negative_data_ids),
                               mer_df = get_mers(pos = cdhit_data,
                                                 pos_id = cdhit_data_ids,
                                                 neg = negative_data,
                                                 neg_id = negative_data_ids),
                               ngrams12 = count_ampgrams(mer_df, 
                                                         ns = c(1, rep(2, 4)),
                                                         ds = list(0, 0, 1, 2, 3)),
                               ngrams3_1 = count_ampgrams(mer_df, 
                                                          ns = c(3, 3),
                                                          ds = list(c(0, 0), c(0, 1))),
                               ngrams3_2 = count_ampgrams(mer_df, 
                                                          ns = c(3, 3),
                                                          ds = list(c(1, 0), c(1, 1))))

make(analysis_AmpGram, seed = 990)

file.copy(from = ".drake", to = paste0(data_path, "drake-cache"), recursive = TRUE, overwrite = TRUE)

ampgram_cache <- drake_cache(paste0(data_path, "drake-cache/.drake"))
ampgram_cache$import(ampgram_cache)

# drake_cache(paste0(data_path, "drake-cache"))$import(drake_cache(paste0(data_path, "drake-cache")))
# vis_drake_graph(drake_config(analysis_AmpGram))