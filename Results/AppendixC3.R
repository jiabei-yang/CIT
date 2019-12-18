library(dplyr)
library(ggplot2)

load("../Data/main.RData")

##################################################################################################
############################################ Table ###############################################
##################################################################################################
res.select <- res %>% 
  filter(Method %in% c("True.NoHonest", "Noisy.NoHonest", "Misspecified.NoHonest", 
                       "InNode.True.CV1", "InNode.Noisy.CV1", "InNode.Misspecified.CV1", 
                       "InNode.TrueIPW.TrueAdj.CV1", "InNode.NoisyIPW.TrueAdj.CV1",
                       "InNode.TrueIPW.NoisyAdj.CV1", "InNode.NoisyIPW.NoisyAdj.CV1",
                       "InNode.MisspecifiedIPW.MisspecifiedAdj.CV1")) %>%
  filter(!is.na(mse))

res.select <- res.select %>%
  mutate(Method = paste(est.mthd, Method, sep = " ")) 

res.select <- res.select %>% 
  mutate(Algorithm = ifelse(est.mthd == "CT", "Best CT", NA)) %>%
  mutate(Algorithm = ifelse(est.mthd == "CT.default", "Original CT", Algorithm)) %>%
  mutate(Algorithm = ifelse(grepl("CV1", Method), "Observational Data IT", Algorithm)) %>%
  mutate(Method = ifelse(Method == "CT Misspecified.NoHonest", "Misspecified CT", Method)) %>%
  mutate(Method = ifelse(Method == "CT Noisy.NoHonest", "CT", Method)) %>%
  mutate(Method = ifelse(Method == "CT True.NoHonest", "True CT", Method)) %>%
  mutate(Method = ifelse(Method == "CT.default Misspecified.NoHonest", "Misspecified CT", Method)) %>%
  mutate(Method = ifelse(Method == "CT.default Noisy.NoHonest", "CT", Method)) %>%
  mutate(Method = ifelse(Method == "CT.default True.NoHonest", "True CT", Method)) %>%
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
  mutate(Method    = factor(Method, levels(res.select$Method)[c(7, 3, 10, 9, 6, 12, 8, 5, 11, 1, 4, 14, 13, 2)])) %>%
  mutate(setting   = factor(setting, levels(res.select$setting)[c(2, 1)])) %>%
  mutate(Algorithm = factor(Algorithm, levels(res.select$Algorithm)[c(3, 1, 2)]))
levels(res.select$Method) <- c("Mis Cov CT", "Mis Func CT", "True CT", 
                               "Mis Cov IPW", "Mis Func IPW", "True IPW",
                               "Mis Cov G", "Mis Func G", "True G",
                               "Both Mis Cov DR", "Both Mis Func DR", "True Prop Mis Func Out DR",
                               "True Out Mis Func Prop DR", "Both True DR")
res.select <- res.select %>%
  mutate(color.method = ifelse(est.mthd == "CT.default", "CT", as.character(est.mthd)))

# Average running time
summ.select <- NULL
summ.select <- cbind(summ.select,
                     tapply(res.select[, 16], 
                            paste(res.select$setting, res.select$Algorithm, res.select$Method, sep = " "), 
                            mean, 
                            na.rm = T))

colnames(summ.select) <- colnames(res.select)[16]
summ.select <- data.frame(summ.select[c(c(15:17, 1:3, c(8, 10, 12, 7, 9, 11, 4, 5, 14, 13, 6)) + 17, c(15:17, 1:3, c(8, 10, 12, 7, 9, 11, 4, 5, 14, 13, 6))), ])
round(summ.select, 2)

# Generate the Latex code
summ.latex <- data.frame(Method = rownames(summ.select),
                         Time   = round(summ.select, 2))
colnames(summ.latex)[2] <- "Time"
summ.latex <- summ.latex %>%
  mutate(Method = sub("Homogeneous ", "", Method)) %>%
  mutate(Method = sub("Heterogeneous ", "", Method))
summ.latex <- summ.latex %>%
  mutate(Estimator = ifelse(grepl("IPW", Method), "IPW", NA)) %>%
  mutate(Estimator = ifelse(grepl("G", Method), "G", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("DR", Method), "DR", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("CT", Method), "CT", Estimator))
summ.latex <- summ.latex %>%
  mutate(Method = sub("Observational Data IT", "", Method)) %>%
  mutate(Method = sub("Original CT", "", Method)) %>%
  mutate(Method = sub("Best CT", "", Method))
summ.latex <- summ.latex %>%
  mutate(Method = sub(" IPW", "", Method)) %>%
  mutate(Method = sub(" G", "", Method)) %>%
  mutate(Method = sub(" DR", "", Method)) %>%
  mutate(Method = sub(" CT", "", Method))
summ.latex <- summ.latex[, c(3, 1:2)]

summ.latex <- data.frame(summ.latex[1:17, ], NA, summ.latex[18:34, 3])
colnames(summ.latex)[5] <- "Time"

print(xtable(summ.latex), include.rownames=FALSE)

