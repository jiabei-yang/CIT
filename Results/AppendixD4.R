library(dplyr)
library(gridExtra)
library(ggplot2)

load("../Data/AppendixD4/RhcSimSubFtrDrAFrst5YLnrFrst5.RData")
load("../Data/AppendixD4/RhcSimSubFtrIpwAFrst5YLnrFrst5.RData")
load("../Data/AppendixD4/RhcSimSubFtrGAFrst5YLnrFrst5.RData")
dim(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5)
dim(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5)
dim(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5)

res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5 <- data.frame(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5)
res.rhcSimSubFtr.g.AFrst5.YLnrFrst5   <- data.frame(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5)
res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5  <- data.frame(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5)

# Subset features, A first 5, Y first 5
colnames(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5)[1:9] <- paste0("hetero_", c("t1", "t2", "size_tree", "size_tree_2ndLast", "mse", 
                                                                            "exact_corr", "numb_noise", "pps", "corr_frst_splt"))
colnames(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5)[1:8 + 9 + 12 * 7] <- paste0("homo_", c("t1", "t2", "size_tree", "size_tree_2ndLast", "mse", 
                                                                                       "exact_corr", "numb_noise", "pps"))

colnames(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5)[1:9] <- paste0("hetero_", c("t1", "t2", "size_tree", "size_tree_2ndLast", "mse", 
                                                                          "exact_corr", "numb_noise", "pps", "corr_frst_splt"))
colnames(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5)[1:8 + 9 + 12 * 7] <- paste0("homo_", c("t1", "t2", "size_tree", "size_tree_2ndLast", "mse", 
                                                                                     "exact_corr", "numb_noise", "pps"))

colnames(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5)[1:12] <- paste0("hetero_", c("t1", "t2", paste0("size_tree_", 1:4), 
                                                                            "size_tree_2ndLast", "mse", "exact_corr", "numb_noise", "pps", "corr_frst_splt"))
colnames(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5)[1:11 + 12 + 12 * 7] <- paste0("homo_", c("t1", "t2", paste0("size_tree_", 1:4), 
                                                                                        "size_tree_2ndLast", "mse", "exact_corr", "numb_noise", "pps"))

# Summary Paper
table(as.numeric(as.character(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5$hetero_exact_corr))) / 1000
table(as.numeric(as.character(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5$homo_exact_corr))) / 1000
table(as.numeric(as.character(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5$hetero_exact_corr))) / 1000
table(as.numeric(as.character(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5$homo_exact_corr))) / 1000
table(as.numeric(as.character(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5$hetero_exact_corr))) / 1000
table(as.numeric(as.character(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5$homo_exact_corr))) / 1000

mean(as.numeric(as.character(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5$homo_numb_noise)))
mean(as.numeric(as.character(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5$homo_numb_noise)))
mean(as.numeric(as.character(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5$homo_numb_noise)))
mean(as.numeric(as.character(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5$hetero_numb_noise)))
mean(as.numeric(as.character(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5$hetero_numb_noise)))
mean(as.numeric(as.character(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5$hetero_numb_noise)))

mean(as.numeric(as.character(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5$homo_pps)))
mean(as.numeric(as.character(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5$homo_pps)))
mean(as.numeric(as.character(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5$homo_pps)))
mean(as.numeric(as.character(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5$hetero_pps)))
mean(as.numeric(as.character(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5$hetero_pps)))
mean(as.numeric(as.character(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5$hetero_pps)))

mean(as.numeric(as.character(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5$hetero_corr_frst_splt)))
mean(as.numeric(as.character(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5$hetero_corr_frst_splt)))
mean(as.numeric(as.character(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5$hetero_corr_frst_splt)))

table(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5[, 10])
table(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5[, 10])
table(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5[, 13])

table(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5[, 9 + 9 + 12 * 7])
table(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5[, 9 + 9 + 12 * 7])
table(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5[, 12 + 12 + 12 * 7])


mse.rhcSimSubFtr.AFrst5.YLnrFrst5 <- data.frame(setting   = rep(c("Heterogeneous", "Homogeneous"), each = 1000 * 3),
                                                Algorithm = rep(rep(c("IPW-CIT", "g-CIT", "DR-CIT"), each = 1000), 2),
                                                mse       = c(as.numeric(as.character(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5$hetero_mse)), 
                                                              as.numeric(as.character(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5$hetero_mse)), 
                                                              as.numeric(as.character(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5$hetero_mse)), 
                                                              as.numeric(as.character(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5$homo_mse)), 
                                                              as.numeric(as.character(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5$homo_mse)), 
                                                              as.numeric(as.character(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5$homo_mse))))
sizeTree.rhcSimSubFtr.AFrst5.YLnrFrst5 <- data.frame(setting = rep(c("Heterogeneous", "Homogeneous"), each = 1000 * 3),
                                                     Algorithm = rep(rep(c("IPW-CIT", "g-CIT", "DR-CIT"), each = 1000), 2),
                                                     size    = c(as.numeric(as.character(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5$hetero_size_tree)),
                                                                 as.numeric(as.character(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5$hetero_size_tree)),
                                                                 as.numeric(as.character(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5$hetero_size_tree_4)), 
                                                                 as.numeric(as.character(res.rhcSimSubFtr.ipw.AFrst5.YLnrFrst5$homo_size_tree)), 
                                                                 as.numeric(as.character(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5$homo_size_tree)), 
                                                                 as.numeric(as.character(res.rhcSimSubFtr.dr.AFrst5.YLnrFrst5$homo_size_tree_4))))

tapply(sizeTree.rhcSimSubFtr.AFrst5.YLnrFrst5$size, paste0(sizeTree.rhcSimSubFtr.AFrst5.YLnrFrst5$setting, ".", sizeTree.rhcSimSubFtr.AFrst5.YLnrFrst5$Algorithm), table)
sizeTree.rhcSimSubFtr.AFrst5.YLnrFrst5 <- sizeTree.rhcSimSubFtr.AFrst5.YLnrFrst5 %>%
  mutate(size_fig = ifelse(size %in% 1:5, size, NA)) %>%
  mutate(size_fig = ifelse(size %in% 6:10, "6-10", size_fig)) %>%
  mutate(size_fig = ifelse(size %in% 11:15, "11-15", size_fig)) %>%
  mutate(size_fig = ifelse(size %in% 16:25, "16-25", size_fig)) %>%
  mutate(size_fig = ifelse(size >= 26, ">25", size_fig))
sizeTree.rhcSimSubFtr.AFrst5.YLnrFrst5 <- sizeTree.rhcSimSubFtr.AFrst5.YLnrFrst5 %>%
  mutate(size_fig = factor(size_fig, c(1:5, "6-10", "11-15", "16-25", ">25")))
tapply(sizeTree.rhcSimSubFtr.AFrst5.YLnrFrst5$size_fig, paste0(sizeTree.rhcSimSubFtr.AFrst5.YLnrFrst5$setting, ".", sizeTree.rhcSimSubFtr.AFrst5.YLnrFrst5$Algorithm), table)

mse.rhcSimSubFtr.AFrst5.YLnrFrst5 <- mse.rhcSimSubFtr.AFrst5.YLnrFrst5 %>%
  mutate(setting   = factor(setting, c("Homogeneous", "Heterogeneous")),
         Algorithm = factor(Algorithm, c("IPW-CIT", "g-CIT", "DR-CIT")))

sizeTree.rhcSimSubFtr.AFrst5.YLnrFrst5 <- sizeTree.rhcSimSubFtr.AFrst5.YLnrFrst5 %>%
  mutate(setting   = factor(setting, c("Homogeneous", "Heterogeneous")),
         Algorithm = factor(Algorithm, c("IPW-CIT", "g-CIT", "DR-CIT")))

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#999999")

p1 <- ggplot(mse.rhcSimSubFtr.AFrst5.YLnrFrst5, aes(x = Algorithm)) +
  geom_boxplot(aes(y = mse)) + 
  scale_y_log10() +
  facet_grid(rows = vars(setting)) + 
  # xlab("Setting") +
  ylab("MSE") +
  ggtitle("MSE") +
  theme_bw()

p2 <- ggplot(sizeTree.rhcSimSubFtr.AFrst5.YLnrFrst5, aes(x = size_fig, fill = setting)) +
  geom_bar(position = position_dodge2(preserve = "single", padding = 0)) + 
  facet_grid(Algorithm ~ .) +
  theme_bw() +
  # scale_x_log10() +
  # scale_y_sqrt(name   = "Count",
  #              breaks = c(0, 10, seq(250, 1000, by = 250)), 
  #              labels = c(0, 10, seq(250, 1000, by = 250)),
  #              limits = c(0, 1000)) +
  ylim(0, 1000) +
  ylab("Count") +
  xlab("Size of trees") +
  ggtitle("Size of trees") +
  # scale_colour_manual(values = cbbPalette, 
  #                     name = "Algorithms", 
  #                     breaks = c("IPW-CIT", "G-CIT", "DR-CIT"),
  #                     labels = c("IPW-CIT", "G-CIT", "DR-CIT")) + 
  scale_fill_manual(values = cbbPalette[4:5], 
                    name = "Setting", 
                    breaks = c("Homogeneous", "Heterogeneous"),
                    labels = c("Homogeneous", "Heterogeneous")) +
  theme(legend.position = "bottom") +
  geom_text(stat='count', 
            aes(label=..count.., group = setting), 
            vjust    = -0.5, 
            position = position_dodge2(preserve = "single", padding = 0, width = 0.9), 
            size = 2.5) 

# Figure S13
grid.arrange(p1, p2, ncol = 2)
