library(dplyr)
library(ggplot2)

load("../Data/AppendixC1.RData")

##################################################################################################
############################################# Figure #############################################
##################################################################################################
res <- res %>%
  mutate(method = paste(split_rule, split_honest, cv_option, cv_honest, sep = "."))
res <- res %>%
  mutate(setting = factor(setting, levels(res$setting)[c(2, 1)]))

ggplot(res, aes(x = method, y = mse)) +
  geom_boxplot(outlier.size = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(setting ~ .) + 
  scale_y_log10(limits = c(2e-9, 2e4)) +
  ylab("MSE") +
  xlab("Method")

##################################################################################################
############################################# Table ##############################################
##################################################################################################
summ <- NULL
for (i in c(6, 7, 10, 11, 18, 19)){
  summ <- cbind(summ,
                tapply(res[, i], 
                       paste(res$setting, res$split_rule, res$split_honest, res$cv_option, res$cv_honest, sep = " "), 
                       mean))
  
}
colnames(summ) <- colnames(res)[c(6, 7, 10, 11, 18, 19)]
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
  select(-corr.frst.splt)

summ.latex <- summ.latex[match(best.combs, summ.latex$Method), ]
print(xtable(summ.latex), include.rownames=FALSE)



