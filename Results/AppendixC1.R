library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)

load("../Data/AppendixC1/CausalTreeSettings.RData")

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

performance.ct.used <- NULL
for (comb.scnr.i in 1:nrow(comb.scnr)) {
  performance.ct.used <- rbind(performance.ct.used,
                               performance.ct[, 8 * (comb.scnr.i - 1) + 1:8])
}
for (comb.scnr.i in 1:nrow(comb.scnr)) {
  performance.ct.used <- rbind(performance.ct.used,
                               cbind(performance.ct[, 8 * 42 + 7 * (comb.scnr.i - 1) + 1:7], NA))
}

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

performance.ct.used <- performance.ct.used %>%
  mutate(setting = factor(setting, c("Homogeneous", "Heterogeneous")))

# Appendix causal tree settings
# Figure S1
ggplot(performance.ct.used, aes(x = Method, y = mse)) +
  geom_boxplot(outlier.size = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(setting ~ .) + 
  scale_y_log10(limits = c(6.22e-13, 1.86e6)) +
  # scale_y_log10() +
  ylab("MSE") +
  xlab("Method")

# Table S1
summ <- NULL
for (i in c(3:4, 7:10)){
  summ <- cbind(summ,
                tapply(performance.ct.used[, i], 
                       paste(performance.ct.used$setting, performance.ct.used$Method, sep = " "), 
                       mean))
  
}
colnames(summ) <- colnames(performance.ct.used)[c(3:4, 7:10)]
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
# table same as before so do not change

