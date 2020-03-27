library(dplyr)

paste_to_five <- function(x) {
  paste0(paste0(rep("0", 5 - nchar(x)), collapse = ""), x)
}
  
paste_to_five(12389)s

read.csv("./data/dbamp_df.csv") %>% 
  mutate(id = 1L:nrow(.)) %>% 
  mutate(id = sapply(id, paste_to_five)) %>% 
  mutate(id = paste0("dbAMP_", id)) %>% 
  select(id, Name, Length, Experimental.Evidence) %>% 
  filter(Experimental.Evidence == "YES") %>% 
  arrange(desc(Length)) %>% 
  write.csv(file = "lengths.csv", row.names = FALSE)
