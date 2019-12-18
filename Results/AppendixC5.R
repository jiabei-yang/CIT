library(dplyr)
library(ggplot2)

##################################################################################################
######################################## Out (Appendix C.5) ######################################
##################################################################################################
load("../Data/AppendixC5Out.RData")

res.select <- res %>% 
  filter(Method %in% c("True.NoHonest", "Noisy.NoHonest", "Misspecified.NoHonest", 
                       "Out.True.CV1", "Out.Noisy.CV1", "Out.Misspecified.CV1", 
                       "Out.TrueIPW.TrueAdj.CV1", "Out.NoisyIPW.TrueAdj.CV1",
                       "Out.TrueIPW.NoisyAdj.CV1", "Out.NoisyIPW.NoisyAdj.CV1",
                       "Out.MisspecifiedIPW.MisspecifiedAdj.CV1",
                       "Out.True.CV2", "Out.Noisy.CV2", "Out.Misspecified.CV2", 
                       "Out.TrueIPW.TrueAdj.CV2", "Out.NoisyIPW.TrueAdj.CV2",
                       "Out.TrueIPW.NoisyAdj.CV2", "Out.NoisyIPW.NoisyAdj.CV2",
                       "Out.MisspecifiedIPW.MisspecifiedAdj.CV2")) 

res.select <- res.select %>%
  mutate(Method = paste(est.mthd, Method, sep = " ")) 

res.select <- res.select %>% 
  mutate(Algorithm = ifelse(est.mthd == "CT", "Best CT", NA)) %>%
  mutate(Algorithm = ifelse(grepl("CV1", Method), "Main FTS", Algorithm)) %>%
  mutate(Algorithm = ifelse(grepl("CV2", Method), "Alternative FTS", Algorithm)) %>%
  mutate(Method = ifelse(Method == "CT Misspecified.NoHonest", "Misspecified CT", Method)) %>%
  mutate(Method = ifelse(Method == "CT Noisy.NoHonest", "CT", Method)) %>%
  mutate(Method = ifelse(Method == "CT True.NoHonest", "True CT", Method)) %>%
  mutate(Method = ifelse(grepl("IPW Out.Misspecified", Method), "Misspecified IPW", Method)) %>%
  mutate(Method = ifelse(grepl("IPW Out.Noisy", Method), "IPW", Method)) %>%
  mutate(Method = ifelse(grepl("IPW Out.True", Method), "True IPW", Method)) %>%
  mutate(Method = ifelse(grepl("G Out.Misspecified", Method), "Misspecified G", Method)) %>%
  mutate(Method = ifelse(grepl("G Out.Noisy", Method), "G", Method)) %>%
  mutate(Method = ifelse(grepl("G Out.True", Method), "True G", Method)) %>%
  mutate(Method = ifelse(grepl("DR Out.NoisyIPW.NoisyAdj", Method), "DR", Method)) %>%
  mutate(Method = ifelse(grepl("DR Out.TrueIPW.NoisyAdj", Method), "True.Propensity DR", Method)) %>%
  mutate(Method = ifelse(grepl("DR Out.NoisyIPW.TrueAdj", Method), "True.Mean DR", Method)) %>%
  mutate(Method = ifelse(grepl("DR Out.TrueIPW.TrueAdj", Method), "Both True DR", Method)) %>%
  mutate(Method = ifelse(grepl("DR Out.MisspecifiedIPW.MisspecifiedAdj", Method), "Both Misspecified DR", Method))             

res.select <- res.select %>%
  mutate(Method    = factor(Method)) %>%
  mutate(Algorithm = factor(Algorithm))
res.select <- res.select %>%
  mutate(Method    = factor(Method, levels(res.select$Method)[c(7, 3, 10, 9, 6, 12, 8, 5, 11, 1, 4, 14, 13, 2)])) %>%
  mutate(setting   = factor(setting, levels(res.select$setting)[c(2, 1)])) %>%
  mutate(Algorithm = factor(Algorithm, levels(res.select$Algorithm)[c(2, 3, 1)]))
levels(res.select$Method) <- c("Mis Cov CT", "Mis Func CT", "True CT", 
                               "Mis Cov IPW", "Mis Func IPW", "True IPW",
                               "Mis Cov G", "Mis Func G", "True G",
                               "Both Mis Cov DR", "Both Mis Func DR", "True Prop Mis Func Out DR",
                               "True Out Mis Func Prop DR", "Both True DR")

cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

ggplot(res.select, aes(Method, mse)) +
  geom_boxplot(outlier.size = 0.9, aes(color = est.mthd)) + 
  ylab("MSE") + 
  facet_grid(setting ~ Algorithm, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10(limits = c(2e-9, 2e4)) +
  scale_colour_manual(values = cbbPalette, 
                      name = "Algorithms", 
                      breaks = c("CT", "IPW", "G", "DR"),
                      labels = c("CT", "IPW", "G", "DR")) + 
  theme(legend.position = "bottom") +
  ggtitle("When models are fitted using the whole dataset")

summ.select <- NULL
for (i in c(5, 8, 9, 16, 17)){
  summ.select <- cbind(summ.select,
                       tapply(res.select[, i], 
                              paste(res.select$setting, res.select$Algorithm, res.select$Method, sep = " "), 
                              mean, 
                              na.rm = T))
  
}
colnames(summ.select) <- colnames(res.select)[c(5, 8, 9, 16, 17)]
summ.select <- summ.select[c(c(5, 7, 9, 4, 6, 8, 1, 2, 11, 10, 3) + 14 + 25, c(5, 7, 9, 4, 6, 8, 1, 2, 11, 10, 3) + 25, c(5, 7, 9, 4, 6, 8, 1, 2, 11, 10, 3) + 14, c(5, 7, 9, 4, 6, 8, 1, 2, 11, 10, 3)), ]

# Generate the table in latex code for paper
summ.latex <- data.frame(Method = rownames(summ.select),
                         round(summ.select, 2))
summ.latex <- summ.latex %>%
  mutate(Method = sub("Homogeneous ", "", Method)) %>%
  mutate(Method = sub("Heterogeneous ", "", Method))
summ.latex <- summ.latex %>%
  mutate(Estimator = ifelse(grepl("IPW", Method), "IPW", NA)) %>%
  mutate(Estimator = ifelse(grepl("G", Method), "G", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("DR", Method), "DR", Estimator))
summ.latex <- summ.latex %>%
  mutate(Method = sub(" IPW", "", Method)) %>%
  mutate(Method = sub(" G", "", Method)) %>%
  mutate(Method = sub(" DR", "", Method)) %>%
  mutate(Method = sub("Main FTS ", "", Method)) %>%
  mutate(Method = sub("Alternative FTS", "", Method))
summ.latex <- summ.latex[, c(7, 1:6)]

summ.latex <- data.frame(summ.latex[1:22, ], NA, summ.latex[23:44, 3:7])
summ.latex <- summ.latex %>%
  select(-corr.frst.splt)

print(xtable(summ.latex), include.rownames=FALSE)

##################################################################################################
####################################### Insplit (Appendix C.5) ###################################
##################################################################################################
load("../Data/AppendixC5Insplit.RData")

res.select <- res %>% 
  filter(Method %in% c("True.NoHonest", "Noisy.NoHonest", "Misspecified.NoHonest", 
                       "InSplit.True.CV1", "InSplit.Noisy.CV1", "InSplit.Misspecified.CV1", 
                       "InSplit.TrueIPW.TrueAdj.CV1", "InSplit.NoisyIPW.TrueAdj.CV1",
                       "InSplit.TrueIPW.NoisyAdj.CV1", "InSplit.NoisyIPW.NoisyAdj.CV1",
                       "InSplit.MisspecifiedIPW.MisspecifiedAdj.CV1",
                       "InSplit.True.CV2", "InSplit.Noisy.CV2", "InSplit.Misspecified.CV2", 
                       "InSplit.TrueIPW.TrueAdj.CV2", "InSplit.NoisyIPW.TrueAdj.CV2",
                       "InSplit.TrueIPW.NoisyAdj.CV2", "InSplit.NoisyIPW.NoisyAdj.CV2",
                       "InSplit.MisspecifiedIPW.MisspecifiedAdj.CV2")) 

res.select <- res.select %>%
  mutate(Method = paste(est.mthd, Method, sep = " ")) 

res.select <- res.select %>% 
  mutate(Algorithm = ifelse(est.mthd == "CT", "Best CT", NA)) %>%
  mutate(Algorithm = ifelse(grepl("CV1", Method), "Main FTS", Algorithm)) %>%
  mutate(Algorithm = ifelse(grepl("CV2", Method), "Alternative FTS", Algorithm)) %>%
  mutate(Method = ifelse(Method == "CT Misspecified.NoHonest", "Misspecified CT", Method)) %>%
  mutate(Method = ifelse(Method == "CT Noisy.NoHonest", "CT", Method)) %>%
  mutate(Method = ifelse(Method == "CT True.NoHonest", "True CT", Method)) %>%
  mutate(Method = ifelse(grepl("IPW InSplit.Misspecified", Method), "Misspecified IPW", Method)) %>%
  mutate(Method = ifelse(grepl("IPW InSplit.Noisy", Method), "IPW", Method)) %>%
  mutate(Method = ifelse(grepl("IPW InSplit.True", Method), "True IPW", Method)) %>%
  mutate(Method = ifelse(grepl("G InSplit.Misspecified", Method), "Misspecified G", Method)) %>%
  mutate(Method = ifelse(grepl("G InSplit.Noisy", Method), "G", Method)) %>%
  mutate(Method = ifelse(grepl("G InSplit.True", Method), "True G", Method)) %>%
  mutate(Method = ifelse(grepl("DR InSplit.NoisyIPW.NoisyAdj", Method), "DR", Method)) %>%
  mutate(Method = ifelse(grepl("DR InSplit.TrueIPW.NoisyAdj", Method), "True.Propensity DR", Method)) %>%
  mutate(Method = ifelse(grepl("DR InSplit.NoisyIPW.TrueAdj", Method), "True.Mean DR", Method)) %>%
  mutate(Method = ifelse(grepl("DR InSplit.TrueIPW.TrueAdj", Method), "Both True DR", Method)) %>%
  mutate(Method = ifelse(grepl("DR InSplit.MisspecifiedIPW.MisspecifiedAdj", Method), "Both Misspecified DR", Method))             

res.select <- res.select %>%
  mutate(Method    = factor(Method)) %>%
  mutate(Algorithm = factor(Algorithm))
res.select <- res.select %>%
  mutate(Method    = factor(Method, levels(res.select$Method)[c(7, 3, 10, 9, 6, 12, 8, 5, 11, 1, 4, 14, 13, 2)])) %>%
  mutate(setting   = factor(setting, levels(res.select$setting)[c(2, 1)])) %>%
  mutate(Algorithm = factor(Algorithm, levels(res.select$Algorithm)[c(2, 3, 1)]))
levels(res.select$Method) <- c("Mis Cov CT", "Mis Func CT", "True CT", 
                               "Mis Cov IPW", "Mis Func IPW", "True IPW",
                               "Mis Cov G", "Mis Func G", "True G",
                               "Both Mis Cov DR", "Both Mis Func DR", "True Prop Mis Func Out DR",
                               "True Out Mis Func Prop DR", "Both True DR")

cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

ggplot(res.select, aes(Method, mse)) +
  geom_boxplot(outlier.size = 0.9, aes(color = est.mthd)) + 
  ylab("MSE") + 
  facet_grid(setting ~ Algorithm, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10(limits = c(2e-9, 2e4)) +
  scale_colour_manual(values = cbbPalette, 
                      name = "Algorithms", 
                      breaks = c("CT", "IPW", "G", "DR"),
                      labels = c("CT", "IPW", "G", "DR")) + 
  theme(legend.position = "bottom") +
  ggtitle("When models are fitted separately in the left and right subgroups")

summ.select <- NULL
for (i in c(5, 8, 9, 16, 17)){
  summ.select <- cbind(summ.select,
                       tapply(res.select[, i], 
                              paste(res.select$setting, res.select$Algorithm, res.select$Method, sep = " "), 
                              mean, 
                              na.rm = T))
  
}
colnames(summ.select) <- colnames(res.select)[c(5, 8, 9, 16, 17)]
summ.select <- summ.select[c(c(5, 7, 9, 4, 6, 8, 1, 2, 11, 10, 3) + 14 + 25, c(5, 7, 9, 4, 6, 8, 1, 2, 11, 10, 3) + 25, c(5, 7, 9, 4, 6, 8, 1, 2, 11, 10, 3) + 14, c(5, 7, 9, 4, 6, 8, 1, 2, 11, 10, 3)), ]

# Generate the table in latex code for paper
summ.latex <- data.frame(Method = rownames(summ.select),
                         round(summ.select, 2))
summ.latex <- summ.latex %>%
  mutate(Method = sub("Homogeneous ", "", Method)) %>%
  mutate(Method = sub("Heterogeneous ", "", Method))
summ.latex <- summ.latex %>%
  mutate(Estimator = ifelse(grepl("IPW", Method), "IPW", NA)) %>%
  mutate(Estimator = ifelse(grepl("G", Method), "G", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("DR", Method), "DR", Estimator))
summ.latex <- summ.latex %>%
  mutate(Method = sub(" IPW", "", Method)) %>%
  mutate(Method = sub(" G", "", Method)) %>%
  mutate(Method = sub(" DR", "", Method)) %>%
  mutate(Method = sub("Main FTS ", "", Method)) %>%
  mutate(Method = sub("Alternative FTS", "", Method))
summ.latex <- summ.latex[, c(7, 1:6)]

summ.latex <- data.frame(summ.latex[1:22, ], NA, summ.latex[23:44, 3:7])
summ.latex <- summ.latex %>%
  select(-corr.frst.splt)

print(xtable(summ.latex), include.rownames=FALSE)



