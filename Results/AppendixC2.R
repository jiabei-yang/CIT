library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)

load("../Data/main/CausalTreeOriginal.RData")
performance.orig.ct <- performance.ct
dim(performance.orig.ct)

load("../Data/main/CausalTreeBest.RData")
performance.best.ct <- performance.ct
dim(performance.best.ct)

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

performance.ct <- cbind(scnrs[rep(1:nrow(scnrs), each = 10^4), ],
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

# Figure S2
ggplot(performance.ct, aes(Method, mse)) +
  geom_boxplot(outlier.size = 0.3, size = 0.3, aes(color = i_honest)) + 
  ylab("MSE") + 
  facet_grid(setting ~ algorithm, scales = "free_x", space = "free_x") +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10(limits = c(6.22e-13, 1.86e6)) +
  # scale_y_log10() +
  scale_colour_manual(values = cbbPalette, 
                      name = "Algorithms") + 
  theme(legend.position = "bottom")


