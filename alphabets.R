measures_alphabets <- readd(measures_alphabets)
saveRDS(measures_alphabets, file = "measures_alphabets.rds")

sorted <- measures_alphabets %>% 
  group_by(alphabet, train_group) %>% 
  summarise(mean_AUC = mean(AUC),
            mean_MCC = mean(MCC),
            mean_sens = mean(sensitivity),
            mean_spec = mean(specificity)) %>% 
  arrange(desc(mean_AUC))

best <- filter(measures_alphabets, alphabet %in% sorted[["alphabet"]][1:10]) %>% 
  mutate(alphabet = factor(alphabet, levels = sorted[["alphabet"]][1:10]),
         train_group = factor(train_group, levels = c("[11,19]_(19,26]", "[11,19]", "(19,26]")))

ggplot(best, aes(x = alphabet, y = AUC)) +
  geom_point() + 
  stat_summary(fun = mean, geom = "point", color = "red", size = 4) +
  facet_wrap( ~ train_group, nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

filter(measures_alphabets, alphabet == "c_de_gw_hkr_afilmv_npqsty" & train_group == "[11,19]_(19,26]") %>% 
  pivot_longer(c(AUC, MCC, sensitivity, specificity), names_to = "Measure", values_to = "Value") %>% 
  ggplot(aes(x = Measure, y = Value)) +
  geom_point() + 
  stat_summary(fun = mean, geom = "point", color = "red", size = 4) 
