library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)

load("../Data/AppendixC9/LnrSpltCtSettings.RData")

# scnrs in causal tree
comb.scnr <- expand.grid(cv_honest    = c(T, F),
                         cv_option    = c("TOT", "matching", "CT", "fit"),
                         split_honest = c(T, F),
                         split_rule   = c("TOT", "CT", "fit", "tstats"))

comb.scnr <- comb.scnr %>%
  mutate(split_honest = ifelse(split_rule %in% c("TOT"), NA, split_honest)) %>%
  mutate(cv_honest    = ifelse(cv_option %in% c("TOT", "matching"), NA, cv_honest))
comb.scnr <- comb.scnr[!duplicated(comb.scnr),]
rownames(comb.scnr) <- 1:nrow(comb.scnr)
comb.scnr <- comb.scnr[, 4:1]

# expect ncol same as the calculated number below
dim(performance.ct)
(7+4) * 42 + 7 * 42

# generate the dataset
performance.ct.used <- NULL
for (comb.scnr.i in 1:nrow(comb.scnr)) {
  performance.ct.used <- rbind(performance.ct.used,
                               performance.ct[, 11 * (comb.scnr.i - 1) + 1:11])
}
for (comb.scnr.i in 1:nrow(comb.scnr)) {
  performance.ct.used <- rbind(performance.ct.used,
                               cbind(performance.ct[, 11 * 42 + 7 * (comb.scnr.i - 1) + 1:7], 
                                     matrix(NA, nrow = 10^3, ncol = 4)))
}

# combine scnrs with performance
comb.scnr <- comb.scnr %>%
  unite("Method", 1:4, sep = ".")
comb.scnr <- comb.scnr[c(rep(1:nrow(comb.scnr), each = 1000), 
                         rep(1:nrow(comb.scnr), each = 1000)), ]
comb.scnr <- cbind(data.frame(setting = c(rep("Heterogeneous", 42 * 1000),
                                          rep("Homogeneous", 42 * 1000))),
                   comb.scnr)

performance.ct.used <- cbind(comb.scnr, performance.ct.used)
colnames(performance.ct.used)[2] <- "Method"
colnames(performance.ct.used)[3:ncol(performance.ct.used)] <- gsub("hetero.tot.tot.", "", colnames(performance.ct.used)[3:ncol(performance.ct.used)])

# relevel factors
performance.ct.used <- performance.ct.used %>%
  mutate(setting = factor(setting, c("Homogeneous", "Heterogeneous")))

# Appendix causal tree settings
ggplot(performance.ct.used, aes(x = Method, y = mse)) +
  geom_boxplot(outlier.size = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(setting ~ .) + 
  # scale_y_log10(limits = c(2e-9, 2e4)) +
  scale_y_log10() +
  ylab("MSE") +
  xlab("Method")

summ <- NULL
for (i in 3:ncol(performance.ct.used)){
  summ <- cbind(summ,
                tapply(performance.ct.used[, i], 
                       paste(performance.ct.used$setting, performance.ct.used$Method, sep = " "), 
                       mean))
  
}
colnames(summ) <- colnames(performance.ct.used)[3:ncol(performance.ct.used)]

# find the best CT
# sort((rank(summ[1:42, "mse"]) + rank(summ[1:42 + 42, "mse"])) / 2)[1:10]
best.combs <- sort((rank(summ[1:42, "mse"]) + rank(summ[1:42 + 42, "mse"]) + 
                      rank(summ[1:42, "numb.noise"]) + rank(summ[1:42 + 42, "numb.noise"]) + 
                      (43 - rank(summ[1:42, "p.corr.3lvl.splts.X4"])) + (43 - rank(summ[1:42, "p.corr.3lvl.splts.X6"]))) / 6)[1:10]
best.combs <- sub("Heterogeneous ", "", names(best.combs))

summ.homo <- summ[1:42 + 42, ]
rownames(summ.homo) <- sub("Homogeneous ", "", rownames(summ.homo))
summ.hetero <- summ[1:42, ]
rownames(summ.hetero) <- sub("Heterogeneous ", "", rownames(summ.hetero))

summ.select <- rbind(summ.homo[rownames(summ.homo) %in% best.combs, ], 
                     summ.hetero[rownames(summ.hetero) %in% best.combs, ])

# generate the Latex code
summ.latex <- data.frame(Method = rownames(summ.select),
                         round(summ.select, 2))
# summ.latex <- summ.latex %>%
#   mutate(Method = sub("Homogeneous ", "", Method)) %>%
#   mutate(Method = sub("Heterogeneous ", "", Method))

summ.latex <- data.frame(summ.latex[1:10, 1:8], NA, summ.latex[11:20, 2:ncol(summ.latex)])
# summ.latex <- summ.latex %>%
#   dplyr::select(-corr.frst.splt)

summ.latex <- summ.latex[match(best.combs, summ.latex$Method), ]
print(xtable(summ.latex), include.rownames=FALSE)

# Regular or Honest
R <- 10^3
load("../Data/AppendixC9/LnrSpltCtOrig.RData")
performance.orig.ct <- performance.ct

load("../Data/AppendixC9/LnrSpltCtBest.RData")
performance.best.ct <- performance.ct

# ncol should be the same as calculated number below
dim(performance.orig.ct)
dim(performance.best.ct)

performance.orig.ct <- rbind(performance.orig.ct[, 11 * 0 + 1:11],
                             performance.orig.ct[, 11 * 1 + 1:11],
                             performance.orig.ct[, 11 * 2 + 1:11],
                             performance.orig.ct[, 11 * 3 + 1:11],
                             performance.orig.ct[, 11 * 4 + 1:11],
                             performance.orig.ct[, 11 * 5 + 1:11],
                             cbind(performance.orig.ct[, 11 * 6 + 7 * 0 + 1:7], 
                                   matrix(NA, nrow = R, ncol = 4)),
                             cbind(performance.orig.ct[, 11 * 6 + 7 * 1 + 1:7],                                     
                                   matrix(NA, nrow = R, ncol = 4)),
                             cbind(performance.orig.ct[, 11 * 6 + 7 * 2 + 1:7],                                     
                                   matrix(NA, nrow = R, ncol = 4)),
                             cbind(performance.orig.ct[, 11 * 6 + 7 * 3 + 1:7],                                     
                                   matrix(NA, nrow = R, ncol = 4)),
                             cbind(performance.orig.ct[, 11 * 6 + 7 * 4 + 1:7],                                     
                                   matrix(NA, nrow = R, ncol = 4)),
                             cbind(performance.orig.ct[, 11 * 6 + 7 * 5 + 1:7],                                     
                                   matrix(NA, nrow = R, ncol = 4)))
performance.best.ct <- rbind(performance.best.ct[, 11 * 0 + 1:11],
                             performance.best.ct[, 11 * 1 + 1:11],
                             performance.best.ct[, 11 * 2 + 1:11],
                             performance.best.ct[, 11 * 3 + 1:11],
                             performance.best.ct[, 11 * 4 + 1:11],
                             performance.best.ct[, 11 * 5 + 1:11],
                             cbind(performance.best.ct[, 11 * 6 + 7 * 0 + 1:7],                                     
                                   matrix(NA, nrow = R, ncol = 4)),
                             cbind(performance.best.ct[, 11 * 6 + 7 * 1 + 1:7],                                     
                                   matrix(NA, nrow = R, ncol = 4)),
                             cbind(performance.best.ct[, 11 * 6 + 7 * 2 + 1:7],                                     
                                   matrix(NA, nrow = R, ncol = 4)),
                             cbind(performance.best.ct[, 11 * 6 + 7 * 3 + 1:7],                                     
                                   matrix(NA, nrow = R, ncol = 4)),
                             cbind(performance.best.ct[, 11 * 6 + 7 * 4 + 1:7],                                     
                                   matrix(NA, nrow = R, ncol = 4)),
                             cbind(performance.best.ct[, 11 * 6 + 7 * 5 + 1:7],                                     
                                   matrix(NA, nrow = R, ncol = 4)))
performance.ct <- rbind(performance.orig.ct, performance.best.ct)
rm(performance.orig.ct, performance.best.ct)
dim(performance.ct)
colnames(performance.ct)
colnames(performance.ct) <- gsub("hetero.propsc.true.nohonest.", "", colnames(performance.ct))

# scenarios
algorithm <- c("Original CT", "Best CT")
setting   <- c("Heterogeneous", "Homogeneous")
Method    <- c("True", "Mis Func", "Unmeasured Cov")
i_honest  <- c("Regular", "Honest")

scnrs <- expand.grid(i_honest, Method, setting, algorithm)
colnames(scnrs) <- c("i_honest", "Method", "setting", "algorithm")

performance.ct <- cbind(scnrs[rep(1:nrow(scnrs), each = R), ],
                        performance.ct)
colnames(performance.ct)

# relevel factors
performance.ct <- performance.ct %>%
  mutate(i_honest  = factor(i_honest, c("Regular", "Honest")),
         Method    = factor(Method, c("Unmeasured Cov", "Mis Func", "True")),
         setting   = factor(setting, c("Homogeneous", "Heterogeneous")),
         algorithm = factor(algorithm, c("Original CT", "Best CT")))

# plot
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

ggplot(performance.ct, aes(Method, mse)) +
  geom_boxplot(outlier.size = 0.3, size = 0.3, aes(color = i_honest)) + 
  ylab("MSE") + 
  facet_grid(setting ~ algorithm, scales = "free_x", space = "free_x") +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  # scale_y_log10(limits = c(2e-9, 2e4)) +
  scale_y_log10() +
  scale_colour_manual(values = cbbPalette, 
                      name = "Algorithms") + 
  theme(legend.position = "bottom")

range(performance.ct$mse)

# original honest better, best regular better

# Comparisons
R <- 10^3
# original ct
load("../Data/AppendixC9/LnrSpltCtOrig.RData")
performance.orig.ct <- performance.ct

# best ct
load("../Data/AppendixC9/LnrSpltCtBest.RData")
performance.best.ct <- performance.ct

# cit
load("../Data/AppendixC9/LnrSpltIpw.RData")
load("../Data/AppendixC9/LnrSpltG.RData")
load("../Data/AppendixC9/LnrSpltDr.RData")

performance.ipw <- performance.ipwInnd
performance.g   <- performance.gInnd
performance.dr  <- performance.drInnd

performance.orig.ct <- rbind(performance.orig.ct[, 11 * 1 + 1:11],
                             performance.orig.ct[, 11 * 3 + 1:11],
                             performance.orig.ct[, 11 * 5 + 1:11],
                             cbind(performance.orig.ct[, 11 * 6 + 7 * 1 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                             cbind(performance.orig.ct[, 11 * 6 + 7 * 3 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                             cbind(performance.orig.ct[, 11 * 6 + 7 * 5 + 1:7], matrix(NA, nrow = R, ncol = 4)))
performance.best.ct <- rbind(performance.best.ct[, 11 * 1 + 1:11],
                             performance.best.ct[, 11 * 3 + 1:11],
                             performance.best.ct[, 11 * 5 + 1:11],
                             cbind(performance.best.ct[, 11 * 6 + 7 * 1 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                             cbind(performance.best.ct[, 11 * 6 + 7 * 3 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                             cbind(performance.best.ct[, 11 * 6 + 7 * 5 + 1:7], matrix(NA, nrow = R, ncol = 4)))

performance.ipw <- rbind(performance.ipw[, 11 * 0 + 1:11],
                         performance.ipw[, 11 * 1 + 1:11],
                         performance.ipw[, 11 * 2 + 1:11],
                         performance.ipw[, 11 * 3 + 1:11],
                         performance.ipw[, 11 * 4 + 1:11],
                         performance.ipw[, 11 * 5 + 1:11],
                         cbind(performance.ipw[, 11 * 6 + 7 * 0 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                         cbind(performance.ipw[, 11 * 6 + 7 * 1 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                         cbind(performance.ipw[, 11 * 6 + 7 * 2 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                         cbind(performance.ipw[, 11 * 6 + 7 * 3 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                         cbind(performance.ipw[, 11 * 6 + 7 * 4 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                         cbind(performance.ipw[, 11 * 6 + 7 * 5 + 1:7], matrix(NA, nrow = R, ncol = 4)))
performance.g <- rbind(performance.g[, 11 * 0 + 1:11],
                       performance.g[, 11 * 1 + 1:11],
                       performance.g[, 11 * 2 + 1:11],
                       performance.g[, 11 * 3 + 1:11],
                       performance.g[, 11 * 4 + 1:11],
                       performance.g[, 11 * 5 + 1:11],
                       cbind(performance.g[, 11 * 6 + 7 * 0 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                       cbind(performance.g[, 11 * 6 + 7 * 1 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                       cbind(performance.g[, 11 * 6 + 7 * 2 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                       cbind(performance.g[, 11 * 6 + 7 * 3 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                       cbind(performance.g[, 11 * 6 + 7 * 4 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                       cbind(performance.g[, 11 * 6 + 7 * 5 + 1:7], matrix(NA, nrow = R, ncol = 4)))
performance.dr <- rbind(performance.dr[, 11 * 0 + 1:11],
                        performance.dr[, 11 * 1 + 1:11],
                        performance.dr[, 11 * 2 + 1:11],
                        performance.dr[, 11 * 3 + 1:11],
                        performance.dr[, 11 * 4 + 1:11],
                        performance.dr[, 11 * 5 + 1:11],
                        performance.dr[, 11 * 6 + 1:11],
                        performance.dr[, 11 * 7 + 1:11],
                        performance.dr[, 11 * 8 + 1:11],
                        performance.dr[, 11 * 9 + 1:11],
                        cbind(performance.dr[, 11 * 10 + 7 * 0 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                        cbind(performance.dr[, 11 * 10 + 7 * 1 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                        cbind(performance.dr[, 11 * 10 + 7 * 2 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                        cbind(performance.dr[, 11 * 10 + 7 * 3 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                        cbind(performance.dr[, 11 * 10 + 7 * 4 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                        cbind(performance.dr[, 11 * 10 + 7 * 5 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                        cbind(performance.dr[, 11 * 10 + 7 * 6 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                        cbind(performance.dr[, 11 * 10 + 7 * 7 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                        cbind(performance.dr[, 11 * 10 + 7 * 8 + 1:7], matrix(NA, nrow = R, ncol = 4)),
                        cbind(performance.dr[, 11 * 10 + 7 * 9 + 1:7], matrix(NA, nrow = R, ncol = 4)))

performance.all <- rbind(performance.orig.ct,
                         performance.best.ct, 
                         performance.ipw,
                         performance.g,
                         performance.dr)
rm(performance.orig.ct, performance.best.ct, performance.ct,
   performance.ipw, performance.g, performance.dr,
   performance.ipwInnd, performance.gInnd, performance.drInnd)
dim(performance.all)
colnames(performance.all)
colnames(performance.all) <- gsub("hetero.propsc.true.honest.", "", colnames(performance.all))

# scenarios
algorithm <- c("Original CT", "Best CT")
setting   <- c("Heterogeneous", "Homogeneous")
Method    <- c("True", "Mis Func", "Unmeasured Cov")

scnrs.ct <- expand.grid(Method, setting, algorithm)
colnames(scnrs.ct) <- c("Method", "setting", "algorithm")

algorithm <- c("Main FTS", "Alternative FTS")
setting   <- c("Heterogeneous", "Homogeneous")
Method   <- c("True IPW-CIT", "Mis Func IPW-CIT", "Unmeasured Cov IPW-CIT", 
              "True g-CIT", "Mis Func g-CIT", "Unmeasured Cov g-CIT", 
              "Both True DR-CIT", "True Out Mis Func Treat DR-CIT", "True Treat Mis Func Out DR-CIT", "Both Mis Func DR-CIT", "Both Unmeasured Cov DR-CIT")
scnrs.cit <- expand.grid(Method, setting, algorithm)
scnrs.cit <- scnrs.cit[c(1:3, 1:3 + 22, 12:14, 12:14 + 22, 
                         4:6, 4:6 + 22, 15:17, 15:17 + 22, 
                         7:11, 7:11 + 22, 18:22, 18:22 + 22), ]
colnames(scnrs.cit) <- c("Method", "setting", "algorithm")

scnrs <- rbind(scnrs.ct, scnrs.cit)
performance.all <- cbind(scnrs[rep(1:nrow(scnrs), each = R), ],
                         performance.all)

# relevel factors
performance.all <- performance.all %>%
  mutate(Method    = factor(Method, c("Unmeasured Cov", "Mis Func", "True", 
                                      "Unmeasured Cov IPW-CIT", "Mis Func IPW-CIT", "True IPW-CIT", 
                                      "Unmeasured Cov g-CIT", "Mis Func g-CIT", "True g-CIT", 
                                      "Both Unmeasured Cov DR-CIT", "Both Mis Func DR-CIT", "True Treat Mis Func Out DR-CIT", "True Out Mis Func Treat DR-CIT", "Both True DR-CIT")),
         setting   = factor(setting, c("Homogeneous", "Heterogeneous")),
         algorithm = factor(algorithm, c("Original CT", "Best CT", "Main FTS", "Alternative FTS")))

performance.all <- performance.all %>%
  mutate(est.mthd = ifelse(algorithm %in% c("Original CT", "Best CT"), "CT", NA)) %>%
  mutate(est.mthd = ifelse(grepl("IPW-CIT", Method), "IPW-CIT", est.mthd)) %>%
  mutate(est.mthd = ifelse(grepl("g-CIT", Method), "g-CIT", est.mthd)) %>%
  mutate(est.mthd = ifelse(grepl("DR-CIT", Method), "DR-CIT", est.mthd))
performance.all <- performance.all %>%
  mutate(est.mthd = factor(est.mthd, c("CT", "IPW-CIT", "g-CIT", "DR-CIT")))

# plot
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

# Figure S9
ggplot(performance.all, aes(Method, mse)) +
  geom_boxplot(outlier.size = 0.6, aes(color = est.mthd)) + 
  ylab("MSE") + 
  ggtitle("When true underlying subgroups do not follow a tree structure") +
  facet_grid(setting ~ algorithm, scales = "free", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  # scale_y_log10(limits = c(2e-9, 2e4)) +
  scale_y_log10() +
  scale_colour_manual(values = cbbPalette, 
                      name = "Algorithms", 
                      breaks = c("CT", "IPW-CIT", "g-CIT", "DR-CIT"),
                      labels = c("CT", "IPW-CIT", "g-CIT", "DR-CIT")) + 
  theme(legend.position = "bottom")

# exact.corr, size = 1 for homogeneous
identical(performance.all$exact.corr[performance.all$setting == "Homogeneous"], as.numeric(performance.all$size.tree[performance.all$setting == "Homogeneous"] == 1))

# Table S10
summ.select <- NULL
for (i in c(5:6, 8, 10:14)) {
  summ.select <- cbind(summ.select,
                       tapply(performance.all[, i],
                              paste(performance.all$setting, performance.all$algorithm, performance.all$Method, sep = " "),
                              mean))
  
}
colnames(summ.select) <- colnames(performance.all)[c(5:6, 8, 10:14)]
summ.select <- summ.select[c(c(3, 1, 2) + 28 + 25,
                             c(3, 1, 2) + 28 + 11,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 28 + 11 + 3,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 28,
                             c(3, 1, 2) + 22 + 3,
                             c(3, 1, 2) + 11,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 14,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2))), ]
round(summ.select, 2)

# generate the Latex code
summ.latex <- data.frame(Method = rownames(summ.select),
                         round(summ.select, 2))
summ.latex <- summ.latex %>%
  mutate(Method = sub("Homogeneous ", "", Method)) %>%
  mutate(Method = sub("Heterogeneous ", "", Method))
summ.latex <- summ.latex %>%
  mutate(Estimator = ifelse(grepl("Original CT", Method), "Original CT", NA)) %>%
  mutate(Estimator = ifelse(grepl("Best CT", Method), "Best CT", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("IPW", Method), "IPW-CIT", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("g", Method), "g-CIT", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("DR", Method), "DR-CIT", Estimator)) 
summ.latex <- summ.latex %>%
  mutate(Method = sub("Original CT", "", Method)) %>%
  mutate(Method = sub("Best CT", "", Method)) %>%
  mutate(Method = sub("Main FTS", "", Method)) %>%
  mutate(Method = sub("Alternative FTS", "", Method))
summ.latex <- summ.latex %>%
  mutate(Method = sub(" IPW-CIT", "", Method)) %>%
  mutate(Method = sub(" g-CIT", "", Method)) %>%
  mutate(Method = sub(" DR-CIT", "", Method)) 
# mutate(Method = sub(" CT", "", Method))
summ.latex <- summ.latex %>%
  mutate(p.corr.3lvl = p.corr.3lvl.splts.X4 + p.corr.3lvl.splts.X6, 
         p.corr.4lvl = p.corr.4lvl.splts.X4 + p.corr.4lvl.splts.X6)
summ.latex <- summ.latex[, c("Estimator", "Method", "exact.corr", "size.tree", "numb.noise", "t", "p.corr.3lvl", "p.corr.4lvl")]


summ.latex <- data.frame(summ.latex[1:28, c("Estimator", "Method", "exact.corr", "size.tree", "numb.noise", "t")], NA, 
                         summ.latex[29:56, 4:ncol(summ.latex)])
# summ.latex <- summ.latex %>%
#   select(-corr.frst.splt)

print(xtable(summ.latex), include.rownames=FALSE)

