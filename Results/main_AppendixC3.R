library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)

# original ct
load("../Data/main/CausalTreeOriginal.RData")
performance.orig.ct <- performance.ct
dim(performance.orig.ct)

# best ct
load("../Data/main/CausalTreeBest.RData")
performance.best.ct <- performance.ct
dim(performance.best.ct)

# ipw
load("../Data/main/MainIpw.RData")
performance.ipw <- performance.ipwInnd
dim(performance.ipw)

# g
load("../Data/main/MainG.RData")
performance.g <- performance.gInnd
dim(performance.g)

# dr 
load("../Data/main/MainDr.RData")
performance.dr <- performance.drInnd
dim(performance.dr)

performance.orig.ct <- rbind(performance.orig.ct[, 8 * 0 + 1:8],
                             performance.orig.ct[, 8 * 2 + 1:8],
                             performance.orig.ct[, 8 * 4 + 1:8],
                             cbind(performance.orig.ct[, 8 * 6 + 7 * 0 + 1:7], NA),
                             cbind(performance.orig.ct[, 8 * 6 + 7 * 2 + 1:7], NA),
                             cbind(performance.orig.ct[, 8 * 6 + 7 * 4 + 1:7], NA))
performance.best.ct <- rbind(performance.best.ct[, 8 * 0 + 1:8],
                             performance.best.ct[, 8 * 2 + 1:8],
                             performance.best.ct[, 8 * 4 + 1:8],
                             cbind(performance.best.ct[, 8 * 6 + 7 * 0 + 1:7], NA),
                             cbind(performance.best.ct[, 8 * 6 + 7 * 2 + 1:7], NA),
                             cbind(performance.best.ct[, 8 * 6 + 7 * 4 + 1:7], NA))

performance.ipw <- rbind(performance.ipw[, 8 * 0 + 1:8],
                         performance.ipw[, 8 * 1 + 1:8],
                         performance.ipw[, 8 * 2 + 1:8],
                         cbind(performance.ipw[, 8 * 3 + 7 * 0 + 1:7], NA),
                         cbind(performance.ipw[, 8 * 3 + 7 * 1 + 1:7], NA),
                         cbind(performance.ipw[, 8 * 3 + 7 * 2 + 1:7], NA))
performance.g <- rbind(performance.g[, 8 * 0 + 1:8],
                       performance.g[, 8 * 1 + 1:8],
                       performance.g[, 8 * 2 + 1:8],
                       cbind(performance.g[, 8 * 3 + 7 * 0 + 1:7], NA),
                       cbind(performance.g[, 8 * 3 + 7 * 1 + 1:7], NA),
                       cbind(performance.g[, 8 * 3 + 7 * 2 + 1:7], NA))
performance.dr <- rbind(performance.dr[, 8 * 0 + 1:8],
                        performance.dr[, 8 * 1 + 1:8],
                        performance.dr[, 8 * 2 + 1:8],
                        performance.dr[, 8 * 3 + 1:8],
                        performance.dr[, 8 * 4 + 1:8],
                        cbind(performance.dr[, 8 * 5 + 7 * 0 + 1:7], NA),
                        cbind(performance.dr[, 8 * 5 + 7 * 1 + 1:7], NA),
                        cbind(performance.dr[, 8 * 5 + 7 * 2 + 1:7], NA),
                        cbind(performance.dr[, 8 * 5 + 7 * 3 + 1:7], NA),
                        cbind(performance.dr[, 8 * 5 + 7 * 4 + 1:7], NA))

performance.all <- rbind(performance.orig.ct, 
                         performance.best.ct, 
                         performance.ipw,
                         performance.g,
                         performance.dr)
rm(performance.orig.ct, performance.best.ct, 
   performance.ipw, performance.g, performance.dr)
dim(performance.all)
colnames(performance.all)
colnames(performance.all) <- gsub("hetero.propsc.true.nohonest.", "", colnames(performance.all))

# scenarios
algorithm <- c("Original CT", "Best CT")
setting   <- c("Heterogeneous", "Homogeneous")
Method    <- c("True", "Mis Func", "Unmeasured Cov")

scnrs.ct <- expand.grid(Method, setting, algorithm)
colnames(scnrs.ct) <- c("Method", "setting", "algorithm")

algorithm <- c("Causal Interaction Tree (CIT)")
setting   <- c("Heterogeneous", "Homogeneous")
Method   <- c("True IPW-CIT", "Mis Func IPW-CIT", "Unmeasured Cov IPW-CIT", 
              "True g-CIT", "Mis Func g-CIT", "Unmeasured Cov g-CIT", 
              "Both True DR-CIT", "True Out Mis Func Treat DR-CIT", "True Treat Mis Func Out DR-CIT", "Both Mis Func DR-CIT", "Both Unmeasured Cov DR-CIT")

scnrs.cit <- expand.grid(Method, setting, algorithm)
scnrs.cit <- scnrs.cit[c(1:3, 12:14, 4:6, 15:17, 7:11, 18:22), ]
colnames(scnrs.cit) <- c("Method", "setting", "algorithm")

scnrs <- rbind(scnrs.ct, scnrs.cit)
performance.all <- cbind(scnrs[rep(1:nrow(scnrs), each = 10^4), ],
                         performance.all)

# relevel factors
performance.all <- performance.all %>%
  mutate(Method    = factor(Method, c("Unmeasured Cov", "Mis Func", "True", 
                                      "Unmeasured Cov IPW-CIT", "Mis Func IPW-CIT", "True IPW-CIT", 
                                      "Unmeasured Cov g-CIT", "Mis Func g-CIT", "True g-CIT", 
                                      "Both Unmeasured Cov DR-CIT", "Both Mis Func DR-CIT", "True Treat Mis Func Out DR-CIT", "True Out Mis Func Treat DR-CIT", "Both True DR-CIT")),
         setting   = factor(setting, c("Homogeneous", "Heterogeneous")),
         algorithm = factor(algorithm, c("Original CT", "Best CT", "Causal Interaction Tree (CIT)")))

# plot
range(performance.all$mse)

# Figure 1
ggplot(performance.all, aes(Method, mse)) +
  geom_boxplot(outlier.size = 0.4, fatten = 2) +
  ylab("MSE") +
  facet_grid(setting ~ algorithm, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10(limits = c(6.22e-13, 1.86e6)) +
  theme(legend.position = "none")

# Table 1
summ.select <- NULL
for (i in c(5, 8, 9, 11)) {
  summ.select <- cbind(summ.select,
                       tapply(performance.all[, i],
                              paste(performance.all$setting, performance.all$algorithm, performance.all$Method, sep = " "),
                              mean))
  
}
colnames(summ.select) <- colnames(performance.all)[c(5, 8, 9, 11)]
summ.select <- summ.select[c(c(c(17, 15, 16), c(3, 1, 2), c(14, 8, 10), c(14, 8, 10) - 1, c(6, 4, 12, 11, 5)) + 17, 
                             c(c(17, 15, 16), c(3, 1, 2), c(14, 8, 10), c(14, 8, 10) - 1, c(6, 4, 12, 11, 5))), ]
round(summ.select, 2)

# Generate the Latex code
summ.latex <- data.frame(Method = rownames(summ.select),
                         round(summ.select, 2))
summ.latex <- summ.latex %>%
  mutate(Method = sub("Homogeneous ", "", Method)) %>%
  mutate(Method = sub("Heterogeneous ", "", Method))
summ.latex <- summ.latex %>%
  mutate(Estimator = ifelse(grepl("IPW", Method), "IPW-CIT", NA)) %>%
  mutate(Estimator = ifelse(grepl("g", Method), "g-CIT", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("DR", Method), "DR-CIT", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("Original CT", Method), "Original CT", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("Best CT", Method), "Best CT", Estimator))
summ.latex <- summ.latex %>%
  mutate(Method = sub("Causal Interaction Tree \\(CIT\\)", "", Method)) %>%
  mutate(Method = sub("Original CT", "", Method)) %>%
  mutate(Method = sub("Best CT", "", Method))
summ.latex <- summ.latex %>%
  mutate(Method = sub(" IPW-CIT", "", Method)) %>%
  mutate(Method = sub(" g-CIT", "", Method)) %>%
  mutate(Method = sub(" DR-CIT", "", Method)) 
  # mutate(Method = sub(" CT", "", Method))
summ.latex <- summ.latex[, c(6, 1:5)]

summ.latex <- data.frame(summ.latex[1:17, ], NA, summ.latex[18:34, 3:6])
summ.latex <- summ.latex %>%
  dplyr::select(-corr.frst.splt)

print(xtable(summ.latex), include.rownames=FALSE)

# time
# Table S2
summ.select <- NULL
for (i in c(10)) {
  summ.select <- cbind(summ.select,
                       tapply(performance.all[, i],
                              paste(performance.all$setting, performance.all$algorithm, performance.all$Method, sep = " "),
                              mean))
  
}
colnames(summ.select) <- colnames(performance.all)[10]
summ.select <- summ.select[c(c(c(17, 15, 16), c(3, 1, 2), c(14, 8, 10), c(14, 8, 10) - 1, c(6, 4, 12, 11, 5)) + 17, 
                             c(c(17, 15, 16), c(3, 1, 2), c(14, 8, 10), c(14, 8, 10) - 1, c(6, 4, 12, 11, 5))), ]
round(summ.select, 2)

# Generate the Latex code
summ.latex <- data.frame(Method = names(summ.select),
                         round(summ.select, 2))
summ.latex <- summ.latex %>%
  mutate(Method = sub("Homogeneous ", "", Method)) %>%
  mutate(Method = sub("Heterogeneous ", "", Method))
summ.latex <- summ.latex %>%
  mutate(Estimator = ifelse(grepl("IPW", Method), "IPW-CIT", NA)) %>%
  mutate(Estimator = ifelse(grepl("g", Method), "g-CIT", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("DR", Method), "DR-CIT", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("Original CT", Method), "Original CT", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("Best CT", Method), "Best CT", Estimator))
summ.latex <- summ.latex %>%
  mutate(Method = sub("Causal Interaction Tree \\(CIT\\)", "", Method)) %>%
  mutate(Method = sub("Original CT", "", Method)) %>%
  mutate(Method = sub("Best CT", "", Method))
summ.latex <- summ.latex %>%
  mutate(Method = sub(" IPW-CIT", "", Method)) %>%
  mutate(Method = sub(" g-CIT", "", Method)) %>%
  mutate(Method = sub(" DR-CIT", "", Method)) 
# mutate(Method = sub(" CT", "", Method))
summ.latex <- summ.latex[, c(3, 1:2)]

summ.latex <- data.frame(summ.latex[1:17, ], NA, summ.latex[18:34, 3])
print(xtable(summ.latex), include.rownames=FALSE)

# Rhc Root Bootstrap Interval
load("../Data/main/RhcRootBi.RData")
trt.eff.ci <- data.frame(apply(trt.eff, 2, quantile, probs = c(0.025, 0.5, 0.975)))
colnames(trt.eff.ci) <- c("IPW", "g", "DR")
round(trt.eff.ci, 3)

# Rhc First Split Bootstrap Interval
load("../Data/main/RhcFrstSpltBi.RData")
apply(dr.ci.frst.splt, 2, quantile, probs = c(0.025, 0.5, 0.975))
