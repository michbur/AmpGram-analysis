if(Sys.info()[["nodename"]] %in% c("amyloid", "phobos", "huawei")) {
  data_path <- "/home/michal/Dropbox/AMP-analysis/AmpGram-analysis/"
}
if(Sys.info()[["nodename"]] == "kasia-MACH-WX9") {
  data_path <- "/home/kasia/Dropbox/AmpGram-analysis/"
}

file.copy(from = paste0(data_path, "drake-cache/.drake"), 
          to = ".drake", recursive = TRUE, overwrite = TRUE)
