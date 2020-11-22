library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)
library(grid)

# best ct
load("../Data/main/CausalTreeBest.RData")
performance.best.ct <- performance.ct
dim(performance.best.ct)

# out
load("../Data/AppendixC5/OutIpw.RData")
load("../Data/AppendixC5/OutG.RData")
load("../Data/AppendixC5/OutDr.RData")

performance.ipw <- performance.ipwOut
performance.g   <- performance.gOut
performance.dr  <- performance.drOut

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

performance.all <- rbind(performance.best.ct, 
                         performance.ipw,
                         performance.g,
                         performance.dr)
rm(performance.best.ct, performance.ct,
   performance.ipw, performance.g, performance.dr,
   performance.ipwOut, performance.gOut, performance.drOut)
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
Method   <- c("True IPW-CIT", "Mis Func IPW-CIT", "Unmeasured Cov IPW-CIT", 
              "True g-CIT", "Mis Func g-CIT", "Unmeasured Cov g-CIT", 
              "Both True DR-CIT", "True Out Mis Func Treat DR-CIT", "True Treat Mis Func Out DR-CIT", "Both Mis Func DR-CIT", "Both Unmeasured Cov DR-CIT")
scnrs.cit <- expand.grid(Method, setting, algorithm)
scnrs.cit <- scnrs.cit[c(1:3, 1:3 + 22, 12:14, 12:14 + 22, 
                         4:6, 4:6 + 22, 15:17, 15:17 + 22, 
                         7:11, 7:11 + 22, 18:22, 18:22 + 22), ]
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
         algorithm = factor(algorithm, c("Optimized CT", "Main FTS", "Alternative FTS")))

performance.all <- performance.all %>%
  mutate(est.mthd = ifelse(algorithm == "Optimized CT", "CT", NA)) %>%
  mutate(est.mthd = ifelse(grepl("IPW-CIT", Method), "IPW-CIT", est.mthd)) %>%
  mutate(est.mthd = ifelse(grepl("g-CIT", Method), "g-CIT", est.mthd)) %>%
  mutate(est.mthd = ifelse(grepl("DR-CIT", Method), "DR-CIT", est.mthd))
performance.all <- performance.all %>%
  mutate(est.mthd = factor(est.mthd, c("CT", "IPW-CIT", "g-CIT", "DR-CIT")))

# plot
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

# Figure S4
p <- ggplot(performance.all, aes(Method, mse)) +
  geom_boxplot(outlier.size = 0.6, aes(color = est.mthd)) + 
  ylab("MSE") + 
  ggtitle("When models are fitted using the whole dataset") +
  facet_grid(setting ~ algorithm, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10(limits = c(6.22e-13, 1.86e6)) +
  # scale_y_log10() +
  scale_colour_manual(values = cbbPalette, 
                      name = "Algorithms") + 
  theme(legend.position = "bottom")
g         <- ggplot_gtable(ggplot_build(p))
strip_x   <- which(grepl('strip-t', g$layout$name))

l <- which(grepl('text', g$grobs[[strip_x[1]]]$grobs[[1]]$childrenOrder))
g$grobs[[strip_x[1]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$fontsize <- 8

grid.draw(g)

# Table S4
summ.select <- NULL
for (i in c(5, 8:11)) {
  summ.select <- cbind(summ.select,
                       tapply(performance.all[, i],
                              paste(performance.all$setting, performance.all$algorithm, performance.all$Method, sep = " "),
                              mean))
  
}
colnames(summ.select) <- colnames(performance.all)[c(5, 8:11)]
summ.select <- summ.select[c(c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 14 + 11 + 11,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 14 + 11,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 11,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2))), ]
round(summ.select, 2)

# generate the Latex code
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
  mutate(Method = sub("Main FTS", "", Method)) %>%
  mutate(Method = sub("Alternative FTS", "", Method))
summ.latex <- summ.latex %>%
  mutate(Method = sub(" IPW-CIT", "", Method)) %>%
  mutate(Method = sub(" g-CIT", "", Method)) %>%
  mutate(Method = sub(" DR-CIT", "", Method)) 
# mutate(Method = sub(" CT", "", Method))
summ.latex <- summ.latex[, c(7, 1:6)]

summ.latex <- data.frame(summ.latex[1:22, ], NA, summ.latex[23:44, 3:ncol(summ.latex)])
summ.latex <- summ.latex %>%
  dplyr::select(-corr.frst.splt)

print(xtable(summ.latex), include.rownames=FALSE)



# In Split
# best ct
load("../Data/main/CausalTreeBest.RData")
performance.best.ct <- performance.ct
dim(performance.best.ct)

# split
load("../Data/AppendixC5/SplitIpw.RData")
load("../Data/AppendixC5/SplitG.RData")
load("../Data/AppendixC5/SplitDr.RData")

performance.ipw <- performance.ipwInsplt
performance.g   <- performance.gInsplt
performance.dr  <- performance.drInsplt

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

performance.all <- rbind(performance.best.ct, 
                         performance.ipw,
                         performance.g,
                         performance.dr)
rm(performance.best.ct, performance.ct,
   performance.ipw, performance.g, performance.dr,
   performance.ipwInsplt, performance.gInsplt, performance.drInsplt)
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
Method   <- c("True IPW-CIT", "Mis Func IPW-CIT", "Unmeasured Cov IPW-CIT", 
              "True g-CIT", "Mis Func g-CIT", "Unmeasured Cov g-CIT", 
              "Both True DR-CIT", "True Out Mis Func Treat DR-CIT", "True Treat Mis Func Out DR-CIT", "Both Mis Func DR-CIT", "Both Unmeasured Cov DR-CIT")
scnrs.cit <- expand.grid(Method, setting, algorithm)
scnrs.cit <- scnrs.cit[c(1:3, 1:3 + 22, 12:14, 12:14 + 22, 
                         4:6, 4:6 + 22, 15:17, 15:17 + 22, 
                         7:11, 7:11 + 22, 18:22, 18:22 + 22), ]
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
         algorithm = factor(algorithm, c("Optimized CT", "Main FTS", "Alternative FTS")))

performance.all <- performance.all %>%
  mutate(est.mthd = ifelse(algorithm == "Optimized CT", "CT", NA)) %>%
  mutate(est.mthd = ifelse(grepl("IPW-CIT", Method), "IPW-CIT", est.mthd)) %>%
  mutate(est.mthd = ifelse(grepl("g-CIT", Method), "g-CIT", est.mthd)) %>%
  mutate(est.mthd = ifelse(grepl("DR-CIT", Method), "DR-CIT", est.mthd))
performance.all <- performance.all %>%
  mutate(est.mthd = factor(est.mthd, c("CT", "IPW-CIT", "g-CIT", "DR-CIT")))

# plot
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

# Figure S5
p <- ggplot(performance.all, aes(Method, mse)) +
  geom_boxplot(outlier.size = 0.6, aes(color = est.mthd)) + 
  ylab("MSE") + 
  ggtitle("When models are fitted separately within each child node") +
  facet_grid(setting ~ algorithm, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10(limits = c(6.22e-13, 1.86e6)) +
  # scale_y_log10() +
  scale_colour_manual(values = cbbPalette, 
                      name = "Algorithms") + 
  theme(legend.position = "bottom")
g         <- ggplot_gtable(ggplot_build(p))
strip_x   <- which(grepl('strip-t', g$layout$name))

l <- which(grepl('text', g$grobs[[strip_x[1]]]$grobs[[1]]$childrenOrder))
g$grobs[[strip_x[1]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$fontsize <- 8

grid.draw(g)

# Table S5
summ.select <- NULL
for (i in c(5, 8:11)) {
  summ.select <- cbind(summ.select,
                       tapply(performance.all[, i],
                              paste(performance.all$setting, performance.all$algorithm, performance.all$Method, sep = " "),
                              mean))
  
}
colnames(summ.select) <- colnames(performance.all)[c(5, 8:11)]
summ.select <- summ.select[c(c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 14 + 11 + 11,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 14 + 11,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2)) + 11,
                             c(c(11, 5, 7), c(11, 5, 7) - 1, c(3, 1, 9, 8, 2))), ]
round(summ.select, 2)

# generate the Latex code
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
  mutate(Method = sub("Main FTS", "", Method)) %>%
  mutate(Method = sub("Alternative FTS", "", Method))
summ.latex <- summ.latex %>%
  mutate(Method = sub(" IPW-CIT", "", Method)) %>%
  mutate(Method = sub(" g-CIT", "", Method)) %>%
  mutate(Method = sub(" DR-CIT", "", Method)) 
# mutate(Method = sub(" CT", "", Method))
summ.latex <- summ.latex[, c(7, 1:6)]

summ.latex <- data.frame(summ.latex[1:22, ], NA, summ.latex[23:44, 3:ncol(summ.latex)])
summ.latex <- summ.latex %>%
  dplyr::select(-corr.frst.splt)

print(xtable(summ.latex), include.rownames=FALSE)




