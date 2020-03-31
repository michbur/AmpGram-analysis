library(drake)
library(dplyr)
library(biogram)
library(pbapply)
source("./functions/degenerate_ngrams.R")

string2list <- function(x) {
  pasted_group <- strsplit(x, "_", fixed = TRUE)[[1]] %>% 
    toupper()
  res <- strsplit(pasted_group, "")
  names(res) <- 1L:length(res)
  res
}

mer_df <- readd(mer_df)
binary_ngrams <- readd(binary_ngrams)
encodings <- create_encodings()


# there will be error 
test_features(filter(mer_df, group == "[11,19]", fold != 5)[["target"]], 
              degenerate_ngrams(binary_ngrams[mer_df[["group"]] == "[11,19]" & mer_df[["fold"]] != 5, ], 
                                string2list(encodings[1])))

# checking if distr_crit is working for one feature genereting NAs:
targets <- filter(mer_df, group == "[11,19]", fold != 5)[["target"]]
deg_ngrams <- degenerate_ngrams(binary_ngrams[mer_df[["group"]] == "[11,19]" & mer_df[["fold"]] != 5, ], 
                                string2list(encodings[1]))
feature_size  <- slam::col_sums(deg_ngrams)
t <- create_feature_target(feature_size[["1.1.1_0.1"]], abs(sum(targets) -1), 0, abs(length(targets) - sum(targets)))
dists <- distr_crit(t[,1], t[,2], criterion = "ig")

# works fine so checking steps after distr_crit
names(dists) <- feature_size[["1.1.1_0.1"]]
which(colnames(deg_ngrams) == "1.1.1_0.1")
feature <- as.matrix(deg_ngrams[,133])

# here are errors in calc_ig:
estm <- calc_criterion(targets, feature, calc_ig)
calc_ig(feature, targets, length(targets), sum(targets))
