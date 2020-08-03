library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)

# best ct
load("../Data/main/CausalTreeBest.RData")
performance.best.ct <- performance.ct
dim(performance.best.ct)

# cv2
load("../Data/AppendixC4/Cv2Ipw.RData")
load("../Data/AppendixC4/Cv2G.RData")
load("../Data/AppendixC4/Cv2Dr.RData")

performance.ipw <- performance.ipwInnd
performance.g   <- performance.gInnd
performance.dr  <- performance.drInnd

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

performance.all <- rbind(performance.best.ct, 
                         performance.ipw,
                         performance.g,
                         performance.dr)
rm(performance.best.ct, performance.ct,
   performance.ipw, performance.g, performance.dr,
   performance.ipwInnd, performance.gInnd, performance.drInnd)
dim(performance.all)
colnames(performance.all)
colnames(performance.all) <- gsub("hetero.propsc.true.nohonest.", "", colnames(performance.all))

# scenarios
algorithm <- c("Best CT")
setting   <- c("Heterogeneous", "Homogeneous")
Method    <- c("True", "Mis Func", "Unmeasured Cov")

scnrs.ct <- expand.grid(Method, setting, algorithm)
colnames(scnrs.ct) <- c("Method", "setting", "algorithm")

algorithm <- c("Alternative FTS")
setting   <- c("Heterogeneous", "Homogeneous")
Method   <- c("True IPW-CIT", "Mis Func IPW-CIT", "Unmeasured Cov IPW-CIT", 
              "True g-CIT", "Mis Func g-CIT", "Unmeasured Cov g-CIT", 
              "Both True DR-CIT", "True Out Mis Func Treat DR-CIT", "True Treat Mis Func Out DR-CIT", "Both Mis Func DR-CIT", "Both Unmeasured Cov DR-CIT")

scnrs.cit <- expand.grid(Method, setting, algorithm)
scnrs.cit <- scnrs.cit[c(1:3, 12:14, 4:6, 15:17, 7:11, 18:22), ]
colnames(scnrs.cit) <- c("Method", "setting", "algorithm")

scnrs <- rbind(scnrs.ct, scnrs.cit)
performance.all <- cbind(scnrs[c(rep(1:nrow(scnrs.ct), each = 10^4), 
                                 rep((nrow(scnrs.ct)+1):nrow(scnrs), each = 10^3)), ],
                         performance.all)

# relevel factors
performance.all <- performance.all %>%
  mutate(Method    = factor(Method, c("Unmeasured Cov", "Mis Func", "True", 
                                      "Unmeasured Cov IPW-CIT", "Mis Func IPW-CIT", "True IPW-CIT", 
                                      "Unmeasured Cov g-CIT", "Mis Func g-CIT", "True g-CIT", 
                                      "Both Unmeasured Cov DR-CIT", "Both Mis Func DR-CIT", "True Treat Mis Func Out DR-CIT", "True Out Mis Func Treat DR-CIT", "Both True DR-CIT")),
         setting   = factor(setting, c("Homogeneous", "Heterogeneous")),
         algorithm = factor(algorithm, c("Best CT", "Alternative FTS")))

performance.all <- performance.all %>%
  mutate(est.mthd = ifelse(algorithm == "Best CT", "CT", NA)) %>%
  mutate(est.mthd = ifelse(grepl("IPW-CIT", Method), "IPW-CIT", est.mthd)) %>%
  mutate(est.mthd = ifelse(grepl("g-CIT", Method), "g-CIT", est.mthd)) %>%
  mutate(est.mthd = ifelse(grepl("DR-CIT", Method), "DR-CIT", est.mthd))
performance.all <- performance.all %>%
  mutate(est.mthd = factor(est.mthd, c("CT", "IPW-CIT", "g-CIT", "DR-CIT")))

# plot
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

# Figure S3
ggplot(performance.all, aes(Method, mse)) +
  geom_boxplot(outlier.size = 0.6, aes(color = est.mthd)) + 
  ylab("MSE") + 
  facet_grid(setting ~ algorithm, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10(limits = c(6.22e-13, 1.86e6)) +
  # scale_y_log10() +
  scale_colour_manual(values = cbbPalette, 
                      name = "Algorithms") + 
  theme(legend.position = "bottom")

# Table S3
summ.select <- NULL
for (i in c(5, 8:11)) {
  summ.select <- cbind(summ.select,
                       tapply(performance.all[, i],
                              paste(performance.all$setting, performance.all$algorithm, performance.all$Method, sep = " "),
                              mean))
  
}
colnames(summ.select) <- colnames(performance.all)[c(5, 8:11)]
summ.select <- summ.select[c(c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 14, 
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2))), ]
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
  mutate(Estimator = ifelse(grepl("DR", Method), "DR-CIT", Estimator)) 
summ.latex <- summ.latex %>%
  mutate(Method = sub("Alternative FTS", "", Method)) 
summ.latex <- summ.latex %>%
  mutate(Method = sub(" IPW-CIT", "", Method)) %>%
  mutate(Method = sub(" g-CIT", "", Method)) %>%
  mutate(Method = sub(" DR-CIT", "", Method)) 
# mutate(Method = sub(" CT", "", Method))
summ.latex <- summ.latex[, c(7, 1:6)]

summ.latex <- data.frame(summ.latex[1:11, ], NA, summ.latex[12:22, 3:ncol(summ.latex)])
summ.latex <- summ.latex %>%
  dplyr::select(-corr.frst.splt)

print(xtable(summ.latex), include.rownames=FALSE)

