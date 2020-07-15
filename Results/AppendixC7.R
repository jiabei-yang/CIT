library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)

load("../Data/AppendixC7/BinCausalTreeSettings.RData")

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

# generate the dataset
performance.ct.used <- NULL
for (comb.scnr.i in 1:nrow(comb.scnr)) {
  performance.ct.used <- rbind(performance.ct.used,
                               performance.ct[, 8 * (comb.scnr.i - 1) + 1:8])
}
for (comb.scnr.i in 1:nrow(comb.scnr)) {
  performance.ct.used <- rbind(performance.ct.used,
                               cbind(performance.ct[, 8 * 42 + 7 * (comb.scnr.i - 1) + 1:7], NA))
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
colnames(performance.ct.used)[3:10] <- gsub("hetero.tot.tot.", "", colnames(performance.ct.used)[3:10])

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
for (i in 3:10){
  summ <- cbind(summ,
                tapply(performance.ct.used[, i], 
                       paste(performance.ct.used$setting, performance.ct.used$Method, sep = " "), 
                       mean))
  
}
colnames(summ) <- colnames(performance.ct.used)[3:10]
best.combs <- sort((rank(summ[1:42, 1]) + rank(summ[1:42 + 42, 1])) / 2)[1:10]
best.combs <- sub("Heterogeneous ", "", names(best.combs))

summ.homo <- summ[1:42 + 42, ]
rownames(summ.homo) <- sub("Homogeneous ", "", rownames(summ.homo))
summ.hetero <- summ[1:42, ]
rownames(summ.hetero) <- sub("Heterogeneous ", "", rownames(summ.hetero))

summ.select <- rbind(summ.homo[rownames(summ.homo) %in% best.combs, ], 
                     summ.hetero[rownames(summ.hetero) %in% best.combs, ])

# Generate the Latex code
summ.latex <- data.frame(Method = rownames(summ.select),
                         round(summ.select, 2))
# summ.latex <- summ.latex %>%
#   mutate(Method = sub("Homogeneous ", "", Method)) %>%
#   mutate(Method = sub("Heterogeneous ", "", Method))

summ.latex <- data.frame(summ.latex[1:10, ], NA, summ.latex[11:20, 2:7])
summ.latex <- summ.latex %>%
  dplyr::select(-corr.frst.splt)

summ.latex <- summ.latex[match(best.combs, summ.latex$Method), ]
print(xtable(summ.latex), include.rownames=FALSE)



# Regular or Honest
load("../Data/AppendixC7/BinCausalTreeOriginal.RData")
performance.orig.ct <- performance.ct

load("../Data/AppendixC7/BinCausalTreeBest.RData")
performance.best.ct <- performance.ct

performance.orig.ct <- rbind(performance.orig.ct[, 8 * 0 + 1:8],
                             performance.orig.ct[, 8 * 1 + 1:8],
                             performance.orig.ct[, 8 * 2 + 1:8],
                             performance.orig.ct[, 8 * 3 + 1:8],
                             performance.orig.ct[, 8 * 4 + 1:8],
                             performance.orig.ct[, 8 * 5 + 1:8],
                             cbind(performance.orig.ct[, 8 * 6 + 7 * 0 + 1:7], NA),
                             cbind(performance.orig.ct[, 8 * 6 + 7 * 1 + 1:7], NA),
                             cbind(performance.orig.ct[, 8 * 6 + 7 * 2 + 1:7], NA),
                             cbind(performance.orig.ct[, 8 * 6 + 7 * 3 + 1:7], NA),
                             cbind(performance.orig.ct[, 8 * 6 + 7 * 4 + 1:7], NA),
                             cbind(performance.orig.ct[, 8 * 6 + 7 * 5 + 1:7], NA))
performance.best.ct <- rbind(performance.best.ct[, 8 * 0 + 1:8],
                             performance.best.ct[, 8 * 1 + 1:8],
                             performance.best.ct[, 8 * 2 + 1:8],
                             performance.best.ct[, 8 * 3 + 1:8],
                             performance.best.ct[, 8 * 4 + 1:8],
                             performance.best.ct[, 8 * 5 + 1:8],
                             cbind(performance.best.ct[, 8 * 6 + 7 * 0 + 1:7], NA),
                             cbind(performance.best.ct[, 8 * 6 + 7 * 1 + 1:7], NA),
                             cbind(performance.best.ct[, 8 * 6 + 7 * 2 + 1:7], NA),
                             cbind(performance.best.ct[, 8 * 6 + 7 * 3 + 1:7], NA),
                             cbind(performance.best.ct[, 8 * 6 + 7 * 4 + 1:7], NA),
                             cbind(performance.best.ct[, 8 * 6 + 7 * 5 + 1:7], NA))
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

performance.ct <- cbind(scnrs[rep(1:nrow(scnrs), each = 10^3), ],
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

# Comparisons
load("../Data/AppendixC7/BinCausalTreeOriginal.RData")
performance.orig.ct <- performance.ct

load("../Data/AppendixC7/BinCausalTreeBest.RData")
performance.best.ct <- performance.ct

# cit
# IPW
load("../Data/AppendixC7/BinMixedIpw.RData")

# G
load("../Data/AppendixC7/BinMixedG.RData")

# DR
load("../Data/AppendixC7/BinMixedDr.RData")
dim(performance.drInnd)

performance.ipw <- performance.ipwInnd
performance.g   <- performance.gInnd
performance.dr <- performance.drInnd

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
                         performance.ipw[, 8 * 3 + 1:8],
                         performance.ipw[, 8 * 4 + 1:8],
                         performance.ipw[, 8 * 5 + 1:8],
                         cbind(performance.ipw[, 8 * 6 + 7 * 0 + 1:7], NA),
                         cbind(performance.ipw[, 8 * 6 + 7 * 1 + 1:7], NA),
                         cbind(performance.ipw[, 8 * 6 + 7 * 2 + 1:7], NA),
                         cbind(performance.ipw[, 8 * 6 + 7 * 3 + 1:7], NA),
                         cbind(performance.ipw[, 8 * 6 + 7 * 4 + 1:7], NA),
                         cbind(performance.ipw[, 8 * 6 + 7 * 5 + 1:7], NA))
performance.g <- rbind(performance.g[, 8 * 0 + 1:8],
                       performance.g[, 8 * 1 + 1:8],
                       performance.g[, 8 * 2 + 1:8],
                       performance.g[, 8 * 3 + 1:8],
                       performance.g[, 8 * 4 + 1:8],
                       performance.g[, 8 * 5 + 1:8],
                       cbind(performance.g[, 8 * 6 + 7 * 0 + 1:7], NA),
                       cbind(performance.g[, 8 * 6 + 7 * 1 + 1:7], NA),
                       cbind(performance.g[, 8 * 6 + 7 * 2 + 1:7], NA),
                       cbind(performance.g[, 8 * 6 + 7 * 3 + 1:7], NA),
                       cbind(performance.g[, 8 * 6 + 7 * 4 + 1:7], NA),
                       cbind(performance.g[, 8 * 6 + 7 * 5 + 1:7], NA))
performance.dr <- rbind(performance.dr[, 8 * 0 + 1:8],
                        performance.dr[, 8 * 1 + 1:8],
                        performance.dr[, 8 * 2 + 1:8],
                        performance.dr[, 8 * 3 + 1:8],
                        performance.dr[, 8 * 4 + 1:8],
                        performance.dr[, 8 * 5 + 1:8],
                        performance.dr[, 8 * 6 + 1:8],
                        performance.dr[, 8 * 7 + 1:8],
                        performance.dr[, 8 * 8 + 1:8],
                        performance.dr[, 8 * 9 + 1:8],
                        cbind(performance.dr[, 8 * 10 + 7 * 0 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 1 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 2 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 3 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 4 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 5 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 6 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 7 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 8 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 9 + 1:7], NA))

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
colnames(performance.all) <- gsub("hetero.propsc.true.nohonest.", "", colnames(performance.all))

# scenarios
algorithm <- c("Original CT", "Best CT")
setting   <- c("Heterogeneous", "Homogeneous")
Method    <- c("True", "Mis Func", "Unmeasured Cov")

scnrs.ct <- expand.grid(Method, setting, algorithm)
colnames(scnrs.ct) <- c("Method", "setting", "algorithm")

algorithm <- c("Main FTS", "Alternative FTS")
setting   <- c("Heterogeneous", "Homogeneous")
Method   <- c("True IPW-CIT", "Mis Func IPW-CIT", "Unmeasured Cov IPW-CIT", 
              "True G-CIT", "Mis Func G-CIT", "Unmeasured Cov G-CIT", 
              "Both True DR-CIT", "True Out Mis Func Prop DR-CIT", "True Prop Mis Func Out DR-CIT", "Both Mis Func DR-CIT", "Both Unmeasured Cov DR-CIT")
scnrs.cit <- expand.grid(Method, setting, algorithm)
scnrs.cit <- scnrs.cit[c(1:3, 1:3 + 22, 12:14, 12:14 + 22, 
                         4:6, 4:6 + 22, 15:17, 15:17 + 22, 
                         7:11, 7:11 + 22, 18:22, 18:22 + 22), ]
colnames(scnrs.cit) <- c("Method", "setting", "algorithm")

scnrs <- rbind(scnrs.ct, scnrs.cit)
performance.all <- cbind(scnrs[rep(1:nrow(scnrs), each = 10^3), ],
                         performance.all)

# relevel factors
performance.all <- performance.all %>%
  mutate(Method    = factor(Method, c("Unmeasured Cov", "Mis Func", "True", 
                                      "Unmeasured Cov IPW-CIT", "Mis Func IPW-CIT", "True IPW-CIT", 
                                      "Unmeasured Cov G-CIT", "Mis Func G-CIT", "True G-CIT", 
                                      "Both Unmeasured Cov DR-CIT", "Both Mis Func DR-CIT", "True Prop Mis Func Out DR-CIT", "True Out Mis Func Prop DR-CIT", "Both True DR-CIT")),
         setting   = factor(setting, c("Homogeneous", "Heterogeneous")),
         algorithm = factor(algorithm, c("Original CT", "Best CT", "Main FTS", "Alternative FTS")))

performance.all <- performance.all %>%
  mutate(est.mthd = ifelse(algorithm %in% c("Original CT", "Best CT"), "CT", NA)) %>%
  mutate(est.mthd = ifelse(grepl("IPW-CIT", Method), "IPW-CIT", est.mthd)) %>%
  mutate(est.mthd = ifelse(grepl("G-CIT", Method), "G-CIT", est.mthd)) %>%
  mutate(est.mthd = ifelse(grepl("DR-CIT", Method), "DR-CIT", est.mthd))
performance.all <- performance.all %>%
  mutate(est.mthd = factor(est.mthd, c("CT", "IPW-CIT", "G-CIT", "DR-CIT")))

# plot
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

# Figure S7
ggplot(performance.all, aes(Method, mse)) +
  geom_boxplot(outlier.size = 0.6, aes(color = est.mthd)) + 
  ylab("MSE") + 
  ggtitle("Binary outcome with both continuous and categorical covariates") +
  facet_grid(setting ~ algorithm, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  # scale_y_log10(limits = c(2e-9, 2e4)) +
  scale_y_sqrt() +
  scale_colour_manual(values = cbbPalette, 
                      name = "Algorithms", 
                      breaks = c("CT", "IPW-CIT", "G-CIT", "DR-CIT"),
                      labels = c("CT", "IPW-CIT", "G-CIT", "DR-CIT")) + 
  theme(legend.position = "bottom")

# Table S7
summ.select <- NULL
for (i in c(5, 8, 9:11)) {
  summ.select <- cbind(summ.select,
                       tapply(performance.all[, i],
                              paste(performance.all$setting, performance.all$algorithm, performance.all$Method, sep = " "),
                              mean))
  
}
colnames(summ.select) <- colnames(performance.all)[c(5, 8, 9:11)]
summ.select <- summ.select[c(c(3, 1, 2) + 28 + 25,
                             c(3, 1, 2) + 28 + 11,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 28 + 11 + 3,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 28,
                             c(3, 1, 2) + 22 + 3,
                             c(3, 1, 2) + 11,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 14,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2))), ]
round(summ.select, 2)

# Generate the Latex code
summ.latex <- data.frame(Method = rownames(summ.select),
                         round(summ.select, 2))
summ.latex <- summ.latex %>%
  mutate(Method = sub("Homogeneous ", "", Method)) %>%
  mutate(Method = sub("Heterogeneous ", "", Method))
summ.latex <- summ.latex %>%
  mutate(Estimator = ifelse(grepl("Original CT", Method), "Original CT", NA)) %>%
  mutate(Estimator = ifelse(grepl("Best CT", Method), "Best CT", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("IPW", Method), "IPW-CIT", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("G", Method), "G-CIT", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("DR", Method), "DR-CIT", Estimator)) 
summ.latex <- summ.latex %>%
  mutate(Method = sub("Original CT", "", Method)) %>%
  mutate(Method = sub("Best CT", "", Method)) %>%
  mutate(Method = sub("Main FTS", "", Method)) %>%
  mutate(Method = sub("Alternative FTS", "", Method))
summ.latex <- summ.latex %>%
  mutate(Method = sub(" IPW-CIT", "", Method)) %>%
  mutate(Method = sub(" G-CIT", "", Method)) %>%
  mutate(Method = sub(" DR-CIT", "", Method)) 
# mutate(Method = sub(" CT", "", Method))
summ.latex <- summ.latex[, c(7, 1:6)]

summ.latex <- data.frame(summ.latex[1:28, ], NA, summ.latex[29:56, 3:ncol(summ.latex)])
summ.latex <- summ.latex %>%
  dplyr::select(-corr.frst.splt)

print(xtable(summ.latex), include.rownames=FALSE)

# Coverage
load("../Data/AppendixC7/BinMixedIpw.RData")
performance.ipw <- performance.ipwInnd

load("../Data/AppendixC7/BinMixedG.RData")
performance.g <- performance.gInnd

load("../Data/AppendixC7/BinMixedDr.RData")
dim(performance.drInnd)
performance.dr <- performance.drInnd

performance.ipw <- rbind(performance.ipw[, 8 * 0 + 1:8],
                         performance.ipw[, 8 * 1 + 1:8],
                         performance.ipw[, 8 * 2 + 1:8],
                         performance.ipw[, 8 * 3 + 1:8],
                         performance.ipw[, 8 * 4 + 1:8],
                         performance.ipw[, 8 * 5 + 1:8],
                         cbind(performance.ipw[, 8 * 6 + 7 * 0 + 1:7], NA),
                         cbind(performance.ipw[, 8 * 6 + 7 * 1 + 1:7], NA),
                         cbind(performance.ipw[, 8 * 6 + 7 * 2 + 1:7], NA),
                         cbind(performance.ipw[, 8 * 6 + 7 * 3 + 1:7], NA),
                         cbind(performance.ipw[, 8 * 6 + 7 * 4 + 1:7], NA),
                         cbind(performance.ipw[, 8 * 6 + 7 * 5 + 1:7], NA))
performance.g <- rbind(performance.g[, 8 * 0 + 1:8],
                       performance.g[, 8 * 1 + 1:8],
                       performance.g[, 8 * 2 + 1:8],
                       performance.g[, 8 * 3 + 1:8],
                       performance.g[, 8 * 4 + 1:8],
                       performance.g[, 8 * 5 + 1:8],
                       cbind(performance.g[, 8 * 6 + 7 * 0 + 1:7], NA),
                       cbind(performance.g[, 8 * 6 + 7 * 1 + 1:7], NA),
                       cbind(performance.g[, 8 * 6 + 7 * 2 + 1:7], NA),
                       cbind(performance.g[, 8 * 6 + 7 * 3 + 1:7], NA),
                       cbind(performance.g[, 8 * 6 + 7 * 4 + 1:7], NA),
                       cbind(performance.g[, 8 * 6 + 7 * 5 + 1:7], NA))
performance.dr <- rbind(performance.dr[, 8 * 0 + 1:8],
                        performance.dr[, 8 * 1 + 1:8],
                        performance.dr[, 8 * 2 + 1:8],
                        performance.dr[, 8 * 3 + 1:8],
                        performance.dr[, 8 * 4 + 1:8],
                        performance.dr[, 8 * 5 + 1:8],
                        performance.dr[, 8 * 6 + 1:8],
                        performance.dr[, 8 * 7 + 1:8],
                        performance.dr[, 8 * 8 + 1:8],
                        performance.dr[, 8 * 9 + 1:8],
                        cbind(performance.dr[, 8 * 10 + 7 * 0 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 1 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 2 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 3 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 4 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 5 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 6 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 7 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 8 + 1:7], NA),
                        cbind(performance.dr[, 8 * 10 + 7 * 9 + 1:7], NA))

performance.all <- rbind(performance.ipw,
                         performance.g,
                         performance.dr)
rm(performance.ipw, performance.g, performance.dr,
   performance.ipwInnd, performance.gInnd, performance.drInnd)

colnames(performance.all)
colnames(performance.all) <- gsub("hetero.ipw.glm.propscinnd.true.cv1.", "", colnames(performance.all))

algorithm <- c("Main FTS", "Alternative FTS")
setting   <- c("Heterogeneous", "Homogeneous")
Method   <- c("True IPW-CIT", "Mis Func IPW-CIT", "Unmeasured Cov IPW-CIT", 
              "True G-CIT", "Mis Func G-CIT", "Unmeasured Cov G-CIT", 
              "Both True DR-CIT", "True Out Mis Func Prop DR-CIT", "True Prop Mis Func Out DR-CIT", "Both Mis Func DR-CIT", "Both Unmeasured Cov DR-CIT")
scnrs.cit <- expand.grid(Method, setting, algorithm)
scnrs.cit <- scnrs.cit[c(1:3, 1:3 + 22, 12:14, 12:14 + 22, 
                         4:6, 4:6 + 22, 15:17, 15:17 + 22, 
                         7:11, 7:11 + 22, 18:22, 18:22 + 22), ]
colnames(scnrs.cit) <- c("Method", "setting", "algorithm")

performance.all <- cbind(scnrs.cit[rep(1:nrow(scnrs.cit), each = 10^3), ],
                         performance.all)

# confidence intervals
load("../Data/AppendixC7/BinMixedIpwCover.RData")
load("../Data/AppendixC7/BinMixedGCover.RData")
load("../Data/AppendixC7/BinMixedDrCover.RData")

ci.ipw <- rbind(ci.ipw[, 6 * 0 + 1:6],
                ci.ipw[, 6 * 1 + 1:6],
                ci.ipw[, 6 * 2 + 1:6],
                ci.ipw[, 6 * 0 + 1:6],
                ci.ipw[, 6 * 1 + 1:6],
                ci.ipw[, 6 * 2 + 1:6],
                cbind(ci.ipw[, 6 * 3 + 3 * 0 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
                cbind(ci.ipw[, 6 * 3 + 3 * 1 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
                cbind(ci.ipw[, 6 * 3 + 3 * 2 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
                cbind(ci.ipw[, 6 * 3 + 3 * 0 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
                cbind(ci.ipw[, 6 * 3 + 3 * 1 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
                cbind(ci.ipw[, 6 * 3 + 3 * 2 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)))

ci.g <- rbind(ci.g[, 6 * 0 + 1:6],
              ci.g[, 6 * 1 + 1:6],
              ci.g[, 6 * 2 + 1:6],
              ci.g[, 6 * 0 + 1:6],
              ci.g[, 6 * 1 + 1:6],
              ci.g[, 6 * 2 + 1:6],
              cbind(ci.g[, 6 * 3 + 3 * 0 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
              cbind(ci.g[, 6 * 3 + 3 * 1 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
              cbind(ci.g[, 6 * 3 + 3 * 2 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
              cbind(ci.g[, 6 * 3 + 3 * 0 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
              cbind(ci.g[, 6 * 3 + 3 * 1 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
              cbind(ci.g[, 6 * 3 + 3 * 2 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)))

ci.dr <- rbind(ci.dr[, 6 * 0 + 1:6],
               ci.dr[, 6 * 1 + 1:6],
               ci.dr[, 6 * 2 + 1:6],
               ci.dr[, 6 * 3 + 1:6],
               ci.dr[, 6 * 4 + 1:6],
               ci.dr[, 6 * 0 + 1:6],
               ci.dr[, 6 * 1 + 1:6],
               ci.dr[, 6 * 2 + 1:6],
               ci.dr[, 6 * 3 + 1:6],
               ci.dr[, 6 * 4 + 1:6],
               cbind(ci.dr[, 6 * 5 + 3 * 0 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
               cbind(ci.dr[, 6 * 5 + 3 * 1 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
               cbind(ci.dr[, 6 * 5 + 3 * 2 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
               cbind(ci.dr[, 6 * 5 + 3 * 3 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
               cbind(ci.dr[, 6 * 5 + 3 * 4 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
               cbind(ci.dr[, 6 * 5 + 3 * 0 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
               cbind(ci.dr[, 6 * 5 + 3 * 1 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
               cbind(ci.dr[, 6 * 5 + 3 * 2 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
               cbind(ci.dr[, 6 * 5 + 3 * 3 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)),
               cbind(ci.dr[, 6 * 5 + 3 * 4 + 1:3], matrix(NA, nrow = 10^3, ncol = 3)))
ci.all <- rbind(ci.ipw, ci.g, ci.dr)
dim(ci.all)

performance.all <- cbind(performance.all, ci.all)
colnames(performance.all)[(ncol(performance.all) - 6 + 1):ncol(performance.all)] <- c(paste0("l", c("025", "50", "975")),
                                                                                      paste0("r", c("025", "50", "975")))

performance.all <- performance.all %>%
  mutate(l025 = ifelse(exact.corr == 0, NA, l025),
         l50  = ifelse(exact.corr == 0, NA, l50),
         l975  = ifelse(exact.corr == 0, NA, l975),
         r025 = ifelse(exact.corr == 0, NA, r025),
         r50  = ifelse(exact.corr == 0, NA, r50),
         r975  = ifelse(exact.corr == 0, NA, r975))

performance.all <- performance.all %>%
  mutate(cover_l = (l025 <= 0.1) & (l975 >= 0.1),
         cover_r = (r025 <= -0.3) & (r975 >= -0.3))

# Table S8
summ.select <- NULL
summ.select <- cbind(summ.select, tapply(performance.all$cover_l, 
                                         paste(performance.all$setting, performance.all$algorithm, performance.all$Method, sep = " "),
                                         function(x) mean(x, na.rm = T)))
summ.select <- cbind(summ.select, tapply(performance.all$cover_r, 
                                         paste(performance.all$setting, performance.all$algorithm, performance.all$Method, sep = " "),
                                         function(x) mean(x, na.rm = T)))

summ.select <- summ.select[c(c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 22 + 11,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 22,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 11,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2))), ]

# Generate the Latex code
summ.latex <- data.frame(Method = rownames(summ.select),
                         round(summ.select, 2))
summ.latex <- summ.latex %>%
  mutate(Method = sub("Homogeneous ", "", Method)) %>%
  mutate(Method = sub("Heterogeneous ", "", Method))
summ.latex <- summ.latex %>%
  mutate(Estimator = ifelse(grepl("IPW", Method), "IPW-CIT", NA)) %>%
  mutate(Estimator = ifelse(grepl("G", Method), "G-CIT", Estimator)) %>%
  mutate(Estimator = ifelse(grepl("DR", Method), "DR-CIT", Estimator)) 
summ.latex <- summ.latex %>%
  mutate(Method = sub("Main FTS", "", Method)) %>%
  mutate(Method = sub("Alternative FTS", "", Method))
summ.latex <- summ.latex %>%
  mutate(Method = sub(" IPW-CIT", "", Method)) %>%
  mutate(Method = sub(" G-CIT", "", Method)) %>%
  mutate(Method = sub(" DR-CIT", "", Method)) 
# mutate(Method = sub(" CT", "", Method))
summ.latex <- summ.latex[, c(4, 1:3)]

summ.latex <- data.frame(summ.latex[1:22, ], summ.latex[23:44, 3:ncol(summ.latex)])
# summ.latex <- summ.latex %>%
#   dplyr::select(-corr.frst.splt)

print(xtable(summ.latex), include.rownames=FALSE)



