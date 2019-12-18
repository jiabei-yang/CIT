library(dplyr)
library(ggplot2)

load("../Data/main.RData")

##################################################################################################
########################################### Figure ###############################################
##################################################################################################
res.ct <- res %>% 
  filter(est.mthd %in% c("CT", "CT.default")) 

res.ct <- res.ct %>% 
  mutate(Algorithm = ifelse(est.mthd == "CT", "Best CT", NA)) %>%
  mutate(Algorithm = ifelse(est.mthd == "CT.default", "Original CT", Algorithm)) %>%
  mutate(Method = ifelse(Method == "Misspecified.Honest", "Mis Cov Honest", Method)) %>%
  mutate(Method = ifelse(Method == "Noisy.Honest", "Mis Func Honest", Method)) %>%
  mutate(Method = ifelse(Method == "True.Honest", "True Honest", Method)) %>%
  mutate(Method = ifelse(Method == "Misspecified.NoHonest", "Mis Cov Regular", Method)) %>%
  mutate(Method = ifelse(Method == "Noisy.NoHonest", "Mis Func Regular", Method)) %>%
  mutate(Method = ifelse(Method == "True.NoHonest", "True Regular", Method)) 

res.ct <- res.ct %>%
  mutate(Method    = factor(Method)) %>%
  mutate(Algorithm = factor(Algorithm))
res.ct <- res.ct %>%
  mutate(Method    = factor(Method, levels(res.ct$Method)[c(2, 4, 6, 1, 3, 5)])) %>%
  mutate(setting   = factor(setting, levels(res.ct$setting)[c(2, 1)])) %>%
  mutate(Algorithm = factor(Algorithm, levels(res.ct$Algorithm)[c(2, 1)]))

res.ct <- res.ct %>%
  mutate(color.method = ifelse(grepl("Regular", Method), "Regular", "Honest"))
res.ct <- res.ct %>%
  mutate(Method = sub(" Regular", "", Method)) %>%
  mutate(Method = sub(" Honest", "", Method))


cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

ggplot(res.ct, aes(Method, mse, group = paste(Method, color.method, sep = "."))) +
  geom_boxplot(outlier.size = 0.3, size = 0.3, aes(color = color.method)) + 
  ylab("MSE") + 
  facet_grid(setting ~ Algorithm, scales = "free_x", space = "free_x") +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10(limits = c(2e-9, 2e4)) +
  scale_colour_manual(values = cbbPalette, 
                      name = "Algorithms", 
                      breaks = c("Honest", "Regular"),
                      labels = c("Honest", "Regular")) + 
  theme(legend.position = "bottom")