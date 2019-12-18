library(dplyr)
library(ggplot2)

load("../Data/AppendixC6_BinMixed.RData")

##################################################################################################
################################## MSE Figure (Appendix C.6) #####################################
##################################################################################################
res.select <- res %>%
  filter(!(Method %in% c("True.Honest", "Noisy.Honest", "Misspecified.Honest")))

res.select <- res.select %>%
  mutate(Method = paste(est.mthd, Method, sep = " ")) 

res.select <- res.select %>% 
  mutate(Algorithm = as.character(est.mthd)) %>%
  mutate(Algorithm = ifelse(grepl("CV1", Method), "Main FTS", Algorithm)) %>%
  mutate(Algorithm = ifelse(grepl("CV2", Method), "Alternative FTS", Algorithm)) %>%
  mutate(Method = ifelse(grepl("CT Misspecified.NoHonest", Method), "Misspecified CT", Method)) %>%
  mutate(Method = ifelse(grepl("CT Noisy.NoHonest", Method), "CT", Method)) %>%
  mutate(Method = ifelse(grepl("CT True.NoHonest", Method), "True CT", Method)) %>%
  mutate(Method = ifelse(grepl("IPW InNode.Misspecified", Method), "Misspecified IPW", Method)) %>%
  mutate(Method = ifelse(grepl("IPW InNode.Noisy", Method), "IPW", Method)) %>%
  mutate(Method = ifelse(grepl("IPW InNode.True", Method), "True IPW", Method)) %>%
  mutate(Method = ifelse(grepl("G InNode.Misspecified", Method), "Misspecified G", Method)) %>%
  mutate(Method = ifelse(grepl("G InNode.Noisy", Method), "G", Method)) %>%
  mutate(Method = ifelse(grepl("G InNode.True", Method), "True G", Method)) %>%
  mutate(Method = ifelse(grepl("DR InNode.NoisyIPW.NoisyAdj", Method), "DR", Method)) %>%
  mutate(Method = ifelse(grepl("DR InNode.TrueIPW.NoisyAdj", Method), "True.Propensity DR", Method)) %>%
  mutate(Method = ifelse(grepl("DR InNode.NoisyIPW.TrueAdj", Method), "True.Mean DR", Method)) %>%
  mutate(Method = ifelse(grepl("DR InNode.TrueIPW.TrueAdj", Method), "Both True DR", Method)) %>%
  mutate(Method = ifelse(grepl("DR InNode.MisspecifiedIPW.MisspecifiedAdj", Method), "Both Misspecified DR", Method))             

res.select <- res.select %>%
  mutate(Method    = factor(Method)) %>%
  mutate(Algorithm = factor(Algorithm))
res.select <- res.select %>%
  mutate(Method = factor(Method, levels(res.select$Method)[c(7, 3, 10, 9, 6, 12, 8, 5, 11, 1, 4, 14, 13, 2)])) %>%
  # mutate(Method = factor(Method, levels(res.select$Method)[c(7, 5, 9, 6, 4, 8, 1, 3, 11, 10, 2)])) %>%
  mutate(setting = factor(setting, levels(res.select$setting)[c(2, 1)])) %>%
  mutate(Algorithm = factor(Algorithm, levels(res.select$Algorithm)[c(4, 2, 3, 1)]))
res.select <- res.select %>%
  mutate(color.method = ifelse(grepl("CT", est.mthd), "CT", as.character(est.mthd)))
levels(res.select$Method) <- c("Mis Cov CT", "Mis Func CT", "True CT", 
                               "Mis Cov IPW", "Mis Func IPW", "True IPW",
                               "Mis Cov G", "Mis Func G", "True G",
                               "Both Mis Cov DR", "Both Mis Func DR", "True Prop Mis Func Out DR",
                               "True Out Mis Func Prop DR", "Both True DR")

cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

ggplot(res.select, aes(Method, mse)) +
  geom_boxplot(outlier.size = 0.9, aes(color = color.method)) + 
  ylab("MSE") + 
  facet_grid(setting ~ Algorithm, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_sqrt() +
  scale_colour_manual(values = cbbPalette, 
                      name = "Algorithms", 
                      labels = c("CT", "IPW CAIT", "G CAIT", "DR CAIT")) + 
  theme(legend.position = "bottom")  +
  ggtitle("Binary outcome with both continuous and categorical covariates")

summ.select <- NULL
for (i in c(5, 8, 9, 16, 17)){
  summ.select <- cbind(summ.select,
                       tapply(res.select[, i], 
                              paste(res.select$setting, res.select$Algorithm, res.select$Method, sep = " "), 
                              mean, 
                              na.rm = T))
  
}
colnames(summ.select) <- colnames(res.select)[c(5, 8, 9, 16, 17)]
summ.select <- summ.select[c(1:3 + 11*2 + 3 + 28, 1:3 + 11 + 28, c(5, 7, 9, 4, 6, 8, 1, 2, 11, 10, 3) + 28 + 14, c(5, 7, 9, 4, 6, 8, 1, 2, 11, 10, 3) + 28, 1:3 + 11*2 + 3, 1:3 + 11, c(5, 7, 9, 4, 6, 8, 1, 2, 11, 10, 3) + 14, c(5, 7, 9, 4, 6, 8, 1, 2, 11, 10, 3)), ]

# Generate the table in latex code for paper
summ.latex <- data.frame(Method = rownames(summ.select),
                         round(summ.select, 2))
summ.latex <- summ.latex %>%
  mutate(Method = sub("Homogeneous ", "", Method)) %>%
  mutate(Method = sub("Heterogeneous ", "", Method))
summ.latex <- summ.latex %>%
  mutate(Estimator = ifelse(grepl("IPW", Method), "IPW", NA)) %>%
  mutate(Estimator = ifelse(grepl("G", Method), "G", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("DR", Method), "DR", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("Original CT", Method), "Original CT", Estimator))%>%
  mutate(Estimator = ifelse(grepl("Best CT", Method), "Best CT", Estimator))
summ.latex <- summ.latex %>%
  mutate(Method = sub(" IPW", "", Method)) %>%
  mutate(Method = sub(" G", "", Method)) %>%
  mutate(Method = sub(" DR", "", Method)) %>%
  mutate(Method = sub("Main FTS ", "", Method)) %>%
  mutate(Method = sub("Alternative FTS ", "", Method)) %>%
  mutate(Method = sub("Original CT ", "", Method)) %>%
  mutate(Method = sub("Best CT ", "", Method)) %>%
  mutate(Method = sub(" CT", "", Method)) 
summ.latex <- summ.latex[, c(7, 1:6)]

summ.latex <- data.frame(summ.latex[1:28, ], NA, summ.latex[29:56, 3:7])
summ.latex <- summ.latex %>%
  select(-corr.frst.splt)

print(xtable(summ.latex), include.rownames=FALSE)
