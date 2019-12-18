library(dplyr)
library(ggplot2)

load("../Data/main.RData")

##################################################################################################
############################################ CV2 #################################################
##################################################################################################
res.cv2 <- res %>% 
  filter(Method %in% c("True.NoHonest", "Noisy.NoHonest", "Misspecified.NoHonest", 
                       "InNode.True.CV2", "InNode.Noisy.CV2", "InNode.Misspecified.CV2", 
                       "InNode.TrueIPW.TrueAdj.CV2", "InNode.NoisyIPW.TrueAdj.CV2",
                       "InNode.TrueIPW.NoisyAdj.CV2", "InNode.NoisyIPW.NoisyAdj.CV2",
                       "InNode.MisspecifiedIPW.MisspecifiedAdj.CV2")) %>%
  filter(!is.na(mse))
res.cv2 <- res.cv2 %>%
  filter(est.mthd != "CT.default")

res.cv2 <- res.cv2 %>%
  mutate(Method = paste(est.mthd, Method, sep = " ")) 

res.cv2 <- res.cv2 %>% 
  mutate(Algorithm = ifelse(est.mthd == "CT", "Best CT", NA)) %>%
  mutate(Algorithm = ifelse(grepl("CV2", Method), "Alternative FTS", Algorithm)) %>%
  mutate(Method = ifelse(Method == "CT Misspecified.NoHonest", "Misspecified CT", Method)) %>%
  mutate(Method = ifelse(Method == "CT Noisy.NoHonest", "CT", Method)) %>%
  mutate(Method = ifelse(Method == "CT True.NoHonest", "True CT", Method)) %>%
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

res.cv2 <- res.cv2 %>%
  mutate(Method    = factor(Method)) %>%
  mutate(Algorithm = factor(Algorithm))
res.cv2 <- res.cv2 %>%
  mutate(Method    = factor(Method, levels(res.cv2$Method)[c(7, 3, 10, 9, 6, 12, 8, 5, 11, 1, 4, 14, 13, 2)])) %>%
  mutate(setting   = factor(setting, levels(res.cv2$setting)[c(2, 1)])) %>%
  mutate(Algorithm = factor(Algorithm, levels(res.cv2$Algorithm)[c(2, 1)]))
levels(res.cv2$Method) <- c("Mis Cov CT", "Mis Func CT", "True CT", 
                            "Mis Cov IPW", "Mis Func IPW", "True IPW",
                            "Mis Cov G", "Mis Func G", "True G",
                            "Both Mis Cov DR", "Both Mis Func DR", "True Prop Mis Func Out DR",
                            "True Out Mis Func Prop DR", "Both True DR")

cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

ggplot(res.cv2, aes(Method, mse)) +
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
  theme(legend.position = "bottom")

summ.select <- NULL
for (i in c(5, 8, 9, 16, 17)){
  summ.select <- cbind(summ.select,
                       tapply(res.cv2[, i], 
                              paste(res.cv2$setting, res.cv2$Algorithm, res.cv2$Method, sep = " "), 
                              mean, 
                              na.rm = T))
  
}
colnames(summ.select) <- colnames(res.cv2)[c(5, 8, 9, 16, 17)]
summ.select <- summ.select[c(c(8, 10, 12, 7, 9, 11, 4, 5, 14, 13, 6) - 3 + 14, c(8, 10, 12, 7, 9, 11, 4, 5, 14, 13, 6) - 3), ]
round(summ.select, 2)

# Generate the Latex code
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
  mutate(Method = sub("Alternative FTS", "", Method))
summ.latex <- summ.latex[, c(7, 1:6)]

summ.latex <- data.frame(summ.latex[1:11, ], NA, summ.latex[12:22, 3:7])
summ.latex <- summ.latex %>%
  select(-corr.frst.splt)

print(xtable(summ.latex), include.rownames=FALSE)
