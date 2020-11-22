library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)

# best ct
load("../Data/main/CausalTreeBest.RData")
performance.best.ct <- performance.ct
dim(performance.best.ct)

# Rf Out
load("../Data/AppendixC6/Rf.RData")
performance.dr  <- performance.drOut

performance.best.ct <- rbind(performance.best.ct[, 8 * 0 + 1:8],
                             performance.best.ct[, 8 * 2 + 1:8],
                             performance.best.ct[, 8 * 4 + 1:8],
                             cbind(performance.best.ct[, 8 * 6 + 7 * 0 + 1:7], NA),
                             cbind(performance.best.ct[, 8 * 6 + 7 * 2 + 1:7], NA),
                             cbind(performance.best.ct[, 8 * 6 + 7 * 4 + 1:7], NA))
performance.dr <- rbind(performance.dr[, 8 * 0 + 1:8],
                        performance.dr[, 8 * 1 + 1:8],
                        performance.dr[, 8 * 2 + 1:8],
                        performance.dr[, 8 * 3 + 1:8],
                        cbind(performance.dr[, 8 * 4 + 7 * 0 + 1:7], NA),
                        cbind(performance.dr[, 8 * 4 + 7 * 1 + 1:7], NA),
                        cbind(performance.dr[, 8 * 4 + 7 * 2 + 1:7], NA),
                        cbind(performance.dr[, 8 * 4 + 7 * 3 + 1:7], NA))

performance.all <- rbind(performance.best.ct, 
                         performance.dr)
rm(performance.best.ct, performance.ct,
   performance.dr, performance.drOut)
dim(performance.all)
colnames(performance.all)
colnames(performance.all) <- gsub("hetero.propsc.true.nohonest.", "", colnames(performance.all))

# scenarios
algorithm <- c("Optimized CT")
setting   <- c("Heterogeneous", "Homogeneous")
Method    <- c("True", "Mis Func", "Unmeasured Cov")

scnrs.ct <- expand.grid(Method, setting, algorithm)
colnames(scnrs.ct) <- c("Method", "setting", "algorithm")

algorithm <- c("Main FTS", "Alternative FTS")
setting   <- c("Heterogeneous", "Homogeneous")
Method   <- c("DR-CIT", "Unmeasured Cov DR-CIT")
scnrs.cit <- expand.grid(Method, setting, algorithm)
scnrs.cit <- scnrs.cit[c(1:2, 1:2 + 4, 3:4, 3:4 + 4), ]
colnames(scnrs.cit) <- c("Method", "setting", "algorithm")

scnrs <- rbind(scnrs.ct, scnrs.cit)
performance.all <- cbind(scnrs[c(rep(1:nrow(scnrs.ct), each = 10^4), 
                                 rep((nrow(scnrs.ct)+1):nrow(scnrs), each = 10^3)), ],
                         performance.all)

# relevel factors
performance.all <- performance.all %>%
  mutate(Method    = factor(Method, c("Unmeasured Cov", "Mis Func", "True", 
                                      "Unmeasured Cov DR-CIT", "DR-CIT")),
         setting   = factor(setting, c("Homogeneous", "Heterogeneous")),
         algorithm = factor(algorithm, c("Optimized CT", "Main FTS", "Alternative FTS")))

performance.all <- performance.all %>%
  mutate(est.mthd = ifelse(algorithm == "Optimized CT", "CT", NA)) %>%
  mutate(est.mthd = ifelse(grepl("DR-CIT", Method), "DR-CIT", est.mthd))
performance.all <- performance.all %>%
  mutate(est.mthd = factor(est.mthd, c("CT", "DR-CIT")))

# plot
# cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")
cbbPalette <- c("#999999", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

# Figure S6
ggplot(performance.all, aes(Method, mse)) +
  geom_boxplot(outlier.size = 0.6, aes(color = est.mthd)) + 
  ylab("MSE") + 
  ggtitle("When models are estimated using random forests") +
  facet_grid(setting ~ algorithm, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  # scale_y_log10(limits = c(6.22e-13, 1.86e6)) +
  scale_y_log10() +
  scale_colour_manual(values = cbbPalette, 
                      name = "Algorithms") + 
  theme(legend.position = "bottom")

# Table S6
summ.select <- NULL
for (i in c(5, 8, 9:11)) {
  summ.select <- cbind(summ.select,
                       tapply(performance.all[, i],
                              paste(performance.all$setting, performance.all$algorithm, performance.all$Method, sep = " "),
                              mean))
  
}
colnames(summ.select) <- colnames(performance.all)[c(5, 8, 9:11)]
summ.select <- summ.select[c(c(2, 1) + 7 + 2,
                             c(2, 1) + 7,
                             c(2, 1) + 2,
                             c(2, 1)), ]
round(summ.select, 2)

# Generate the Latex code
summ.latex <- data.frame(Method = rownames(summ.select),
                         round(summ.select, 2))
summ.latex <- summ.latex %>%
  mutate(Method = sub("Homogeneous ", "", Method)) %>%
  mutate(Method = sub("Heterogeneous ", "", Method))
# summ.latex <- summ.latex %>%
#   # mutate(Estimator = ifelse(grepl("IPW", Method), "IPW-CIT", NA)) %>%
#   # mutate(Estimator = ifelse(grepl("G", Method), "G-CIT", Estimator)) %>%
#   mutate(Estimator = ifelse(grepl("DR", Method), "DR-CIT", Estimator)) 
summ.latex <- summ.latex %>%
  mutate(Method = sub("Main FTS", "", Method)) %>%
  mutate(Method = sub("Alternative FTS", "", Method))
# summ.latex <- summ.latex %>%
#   # mutate(Method = sub(" IPW-CIT", "", Method)) %>%
#   # mutate(Method = sub(" G-CIT", "", Method)) %>%
#   mutate(Method = sub(" DR-CIT", "", Method)) 
# mutate(Method = sub(" CT", "", Method))
# summ.latex <- summ.latex[, c(7, 1:6)]

summ.latex <- data.frame(summ.latex[1:4, ], NA, summ.latex[5:8, 2:ncol(summ.latex)])
summ.latex <- summ.latex %>%
  dplyr::select(-corr.frst.splt)

print(xtable(summ.latex), include.rownames=FALSE)





