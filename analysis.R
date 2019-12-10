library(dplyr)
library(magrittr)
library(drake)
library(seqinr)

source("./functions/raw_data.R")
source("./functions/cdhit_data.R")
source("./functions/nonstandard_AMPs.R")

#filter(dbamp_df, Sequence == dbamp_df[["Sequence"]][duplicated(dbamp_df[["Sequence"]])])




analysis_AmpGram <- drake_plan(raw_data = read_raw_data(),
                               nonstandard_AMPs = analyze_nonstandard_AMPs(raw_data),
                               cdhit_data = filter_cdhit(raw_data))

make(analysis_AmpGram)
