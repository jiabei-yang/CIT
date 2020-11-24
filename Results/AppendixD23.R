library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

load("../Data/AppendixD2/RhcBootIpwSplit.RData")
rhc.boot.ipw <- rhc.boot.ipw.result
load("../Data/AppendixD2/RhcBootGSplit.RData")
rhc.boot.g <- rhc.boot.g.result
load("../Data/AppendixD2/RhcBootDrSplit.RData")
rhc.boot.dr <- rhc.boot.dr.result

rm(rhc.boot.ipw.result, rhc.boot.dr.result, rhc.boot.g.result)

# combine results from CITs
rhc.boot.cit <- rbind(cbind("IPW-CIT", rhc.boot.ipw),
                      cbind("g-CIT", rhc.boot.g),
                      cbind("DR-CIT", rhc.boot.dr))
rhc.boot.cit <- data.frame(rhc.boot.cit)
colnames(rhc.boot.cit)[1:4] <- c("estimator", "t1", "t2", "size_tree")

rhc.boot.cit <- rhc.boot.cit %>%
  mutate(size_tree = as.numeric(levels(size_tree))[size_tree])
rhc.boot.cit <- rhc.boot.cit %>%
  mutate(estimator = factor(estimator, c("IPW-CIT", "g-CIT", "DR-CIT")))

rhc.boot.cit <- rhc.boot.cit %>%
  mutate(size_tree_fig = ifelse(size_tree == 1, "1", NA)) %>%
  mutate(size_tree_fig = ifelse((size_tree >= 2) & (size_tree <= 10), "2 - 10", size_tree_fig)) %>%
  mutate(size_tree_fig = ifelse((size_tree >= 11) & (size_tree <= 30), "11 - 30", size_tree_fig)) %>%
  mutate(size_tree_fig = ifelse((size_tree >= 31) & (size_tree <= 60), "31 - 60", size_tree_fig)) %>%
  mutate(size_tree_fig = ifelse(size_tree >= 61, "> 60", size_tree_fig))
rhc.boot.cit <- rhc.boot.cit %>%
  mutate(size_tree_fig = factor(size_tree_fig, c("1", "2 - 10", "11 - 30", "31 - 60", "> 60")))

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#999999")

# Figure S10
ggplot(rhc.boot.cit, aes(x = size_tree_fig, fill = estimator)) +
  geom_bar(position = position_dodge2(preserve = "single", padding = 0)) + 
  # facet_grid(estimator ~ .) +
  theme_bw() +
  # scale_x_log10(breaks = c(1:3, 5, 10, 50, 100)) +
  # scale_y_sqrt(name   = "Count",
  #              breaks = c(0, 10, seq(250, 1000, by = 250)), 
  #              labels = c(0, 10, seq(250, 1000, by = 250)),
  #              limits = c(0, 1000)) +
  # ylim(0, 1000) +
  scale_y_continuous(breaks = seq(0, 1000, 250), 
                     limits = c(0, 1030)) +
  ylab("Count") +
  xlab("Size of trees") +
  # scale_colour_manual(values = cbbPalette, 
  #                     name = "Algorithms", 
  #                     breaks = c("IPW-CIT", "g-CIT", "DR-CIT"),
  #                     labels = c("IPW-CIT", "g-CIT", "DR-CIT")) + 
  scale_fill_manual(values = cbbPalette, 
                    name = "Algorithms", 
                    breaks = c("IPW-CIT", "g-CIT", "DR-CIT"),
                    labels = c("IPW-CIT", "g-CIT", "DR-CIT")) +
  theme(legend.position = "bottom") +
  annotate(geom = "text", x = c(0.7, 1, 1.3), y = c(828, 997, 879) + 25, label = c(828, 997, 879)) 

tapply(rhc.boot.cit$size_tree_fig, rhc.boot.cit$estimator, table)

# bootstrap confidence interval
load("../Data/AppendixD3/RhcSpltRootCi.RData")
ci.ipw <- ci.cit[c(rhc.boot.cit$size_tree[rhc.boot.cit$estimator == "IPW-CIT"] == 1), 1:3]
ci.g   <- ci.cit[c(rhc.boot.cit$size_tree[rhc.boot.cit$estimator == "g-CIT"] == 1), 4:6]
ci.dr  <- ci.cit[c(rhc.boot.cit$size_tree[rhc.boot.cit$estimator == "DR-CIT"] == 1), 7:9]

apply(ci.ipw, 2, median)
apply(ci.g, 2, median)
apply(ci.dr, 2, median)
apply(ci.dr, 2, mean)

ci.cit <- rbind(cbind("IPW-CIT", ci.ipw),
                cbind("g-CIT", ci.g),
                cbind("DR-CIT", ci.dr))
ci.cit <- data.frame(ci.cit)
colnames(ci.cit) <- c("estimator", "ci025", "ci50", "ci975")

ci.cit <- ci.cit %>%
  mutate(ci025 = as.numeric(levels(ci025))[ci025],
         ci50  = as.numeric(levels(ci50))[ci50],
         ci975 = as.numeric(levels(ci975))[ci975])
ci.cit <- ci.cit %>%
  mutate(estimator = factor(estimator, c("IPW-CIT", "g-CIT", "DR-CIT")))

ci.cit <- ci.cit %>%
  gather(limit, value, ci025:ci975)
ci.cit <- ci.cit %>%
  mutate(limit = factor(limit, 
                        levels = c("ci025", "ci50", "ci975"),
                        labels = c("2.5%", "Median", "97.5%")))

# Figure S11
ggplot(ci.cit, aes(x = value, fill = estimator)) +
  geom_histogram(position = position_dodge2(preserve = "single", padding = 0)) + 
  facet_grid(limit ~ .) +
  theme_bw() +
  ylab("Count") +
  xlab("Treatment effect") +
  scale_fill_manual(values = cbbPalette, 
                    name = "Algorithms", 
                    breaks = c("IPW-CIT", "g-CIT", "DR-CIT"),
                    labels = c("IPW-CIT", "g-CIT", "DR-CIT")) +
  theme(legend.position = "bottom") 

# Lower, upper bounds of the treatment effect estimates in two terminal nodes 
load("../Data/AppendixD3/RhcSpltFrstSpltCi.RData")
ci.frst.splt.cit <- data.frame(ci.frst.splt.cit)
colnames(ci.frst.splt.cit) <- c("include",
                                paste0(rep(c("ipw_", "g_", "dr_"), each = 3 * 3), 
                                       rep(rep(c("l_", "r_", "diff_"), each = 3), 3), 
                                       rep(rep(c("025", "50", "975"), 3), 3)))
ci.frst.splt.cit <- ci.frst.splt.cit %>%
  filter(include == 1) %>%
  select(starts_with("dr"))

apply(ci.frst.splt.cit, 2, median)
apply(ci.frst.splt.cit, 2, mean)

ci.frst.splt.cit <- ci.frst.splt.cit %>%
  dplyr::select(-dr_diff_025, -dr_diff_50, -dr_diff_975)
ci.frst.splt.cit <- ci.frst.splt.cit %>%
  gather(limit, value)
ci.frst.splt.cit <- ci.frst.splt.cit %>%
  mutate(node = ifelse(grepl("_l_", limit), "P(survival) < 0.85", "P(survival) >= 0.85")) %>%
  mutate(bound = ifelse(grepl("025", limit), "2.5%", NA)) %>%
  mutate(bound = ifelse(grepl("50", limit), "Median", bound)) %>%
  mutate(bound = ifelse(grepl("975", limit), "97.5%", bound)) 

ci.frst.splt.cit <- ci.frst.splt.cit %>%
  mutate(node = factor(node, c("P(survival) < 0.85", "P(survival) >= 0.85"))) %>%
  mutate(bound = factor(bound, c("2.5%", "Median", "97.5%")))

# Figure S12
ggplot(ci.frst.splt.cit, aes(x = value, fill = node)) +
  geom_histogram(position = position_dodge2(preserve = "single", padding = 0)) + 
  facet_grid(bound ~ .) +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 5.5, 11, 16.5),
                     labels = paste0(round(c(0, 5.5, 11, 16.5) / 44, 4) * 100, "%")) +
  ylab("Count") +
  xlab("Treatment effect") +
  scale_fill_manual(values = cbbPalette[6:7], 
                    name = "Node", 
                    breaks = c("P(survival) < 0.85", "P(survival) >= 0.85")) +
  theme(legend.position = "bottom") 

