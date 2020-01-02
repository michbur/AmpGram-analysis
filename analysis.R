library(dplyr)
library(magrittr)
library(drake)
library(biogram)
library(pbapply)
library(ranger)
library(cvTools)
library(visNetwork)
library(hmeasure)

if(Sys.info()[["nodename"]] == "amyloid") {
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

#filter(dbamp_df, Sequence == dbamp_df[["Sequence"]][duplicated(dbamp_df[["Sequence"]])])

analysis_AmpGram <- drake_plan(raw_data = read_raw_data(),
                               nonstandard_AMPs = analyze_nonstandard_AMPs(raw_data),
                               cdhit_data = filter_cdhit(raw_data),
                               negative_data = read_and_cut(data_path, lengths(cdhit_data)),
                               cdhit_data_ids = generate_holdout_groups(cdhit_data),
                               negative_data_ids = generate_holdout_groups(negative_data),
                               benchmark_file = write_benchmark(pos = cdhit_data,
                                                                pos_id = cdhit_data_ids,
                                                                neg = negative_data,
                                                                neg_id = negative_data_ids))

make(analysis_AmpGram, seed = 990)

# vis_drake_graph(drake_config(analysis_AmpGram))
