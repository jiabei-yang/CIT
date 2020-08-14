#!/usr/bin/env Rscript
library(causalTree)

load("../seed1000.rda")

setwd("../")
folder <- paste(getwd(), "/Functions/", sep="")
functions <- list.files(folder)
functions <- grep("main.R", functions, invert = T, value = T)
functions <- paste(folder, functions, sep = "")
for (i in functions){
  source(i)
}

setwd("Appendix/")

# Rscript running part
registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))

# Read in the arguments from command line
option_list = list(
  make_option(c("-s", "--start"), type="integer"),
  make_option(c("-e", "--end"), type = "integer"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
start <- opt$start
end   <- opt$end

# Data Preparation
rhc.orig <- read.csv("../Data/RHC-cleanup/rhc.csv")

# Remove the first column since it contains unnecessary row numbers
rhc.used <- rhc.orig[, -1]

rhc.used <- rhc.used %>%
  mutate(SadmdDte = as.Date(sadmdte,  origin = "1960-01-01"),
         DschDte  = as.Date(dschdte,  origin = "1960-01-01"),
         DthDte   = as.Date(dthdte,   origin = "1960-01-01"),
         LstctDte = as.Date(lstctdte, origin = "1960-01-01"))

# Factor
for (i in 9:20) {
  rhc.used[, i] <- as.factor(rhc.used[, i])
}

rhc.used <- rhc.used %>%
  mutate(adld3p = as.factor(adld3p))

# Replace NA in cat2 since they have practical meaning
rhc.used <- rhc.used %>%
  mutate(imputed.cat2 = as.factor(ifelse(is.na(cat2), "NA", as.character(cat2))))

# Change factor treatment and outcome to numeric treatment
rhc.used <- rhc.used %>%
  mutate(num.swang1 = ifelse(swang1 == "No RHC", 0, 1))
rhc.used <- rhc.used %>%
  mutate(num.death = ifelse(death == "No", 0, 1))
rhc.used <- rhc.used %>%
  mutate(num.dth30 = ifelse(dth30 == "No", 0, 1))

rhc.used <- rhc.used %>%
  mutate(sqrt.wblc1 = sqrt(wblc1),
         log.alb1   = log(alb1),
         log.bili1  = log(bili1),
         log.crea1  = log(crea1),
         sqrt.urin1 = sqrt(urin1))

# Create wt0 variable
num.col <- ncol(rhc.used)
rhc.used <- rhc.used %>%
  mutate(wt0 = ifelse(wtkilo1 == 0, 1, 0))

rhc.used <- data.frame(A = rhc.used$num.swang1,
                       Y = rhc.used$num.dth30,
                       rhc.used[, c("cat1", "imputed.cat2", "ca", 
                                    "cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx", 
                                    "renalhx", "liverhx", "gibledhx", "malighx", "immunhx",
                                    "transhx", "amihx", "age", "sex", "edu", 
                                    "surv2md1", "das2d3pc", "aps1", "scoma1", "meanbp1", 
                                    "wblc1", "hrt1", "resp1", "temp1", "pafi1", 
                                    "alb1", "hema1", "bili1", "crea1", "sod1", 
                                    "pot1", "paco21", "ph1", "wtkilo1", "dnr1",
                                    "ninsclas", "resp", "card", "neuro", "gastr",
                                    "renal", "meta", "hema", "seps", "trauma",
                                    "ortho", "race", "income", "wt0")])

# fit.rhc.y.full <- lm(Y ~ ., data = rhc.used)
# fit.rhc.y.null <- lm(Y ~ A, data = rhc.used)
# 
# fit.rhc.y.select <- step(fit.rhc.y.null, scope = list(lower = fit.rhc.y.null, upper = fit.rhc.y.full), direction = "both")
fit.rhc.y.select <- lm(Y ~ surv2md1 + dnr1 + bili1 + cat1 + imputed.cat2 + A, data = rhc.used)
round(summary(fit.rhc.y.select)$coef, 2)

# fit.rhc.trt.full <- glm(A ~ ., family = binomial(link = "logit"), data = rhc.used[, colnames(rhc.used) != "Y"])
# fit.rhc.trt.null <- glm(A ~ 1, family = binomial(link = "logit"), data = rhc.used[, colnames(rhc.used) != "Y"])
# 
# fit.rhc.trt.select <- step(fit.rhc.trt.null, scope = list(lower = fit.rhc.trt.null, upper = fit.rhc.trt.full), direction = "both")
fit.rhc.trt.select <- glm(A ~ cat1 + pafi1 + meanbp1 + imputed.cat2 + resp1, family = binomial(link = "logit"), data = rhc.used)
round(summary(fit.rhc.trt.select)$coef, 2)

res.rhcSimSubFtr.g.AFrst5.YLnrFrst5 <- 
  foreach(i = start:end,
          .combine = "rbind") %dopar% {
            
            # set seed
            set.seed(a[i]) 
            
            # treatment, simulated first 5
            rhc.used.sim <- rhc.used %>%
              mutate(prob.sim.A = predict(fit.rhc.trt.select, type = "response"))
            rhc.used.sim <- rhc.used.sim %>%
              mutate(sim.A = rbinom(nrow(rhc.used.sim), size = 1, prob = prob.sim.A))
            
            # outcome, linear first 5
            coef.rhc.y.select <- coefficients(fit.rhc.y.select)
            lnr.y <- model.matrix(fit.rhc.y.select) %*% coef.rhc.y.select
            
            # rescale probability y to range [0.05, 0.95]
            # no effect
            lnr.noeff.y <- (lnr.y - min(lnr.y)) / (max(lnr.y) - min(lnr.y)) * (0.95 - 0.05) + 0.05
            # effect, signal 0.4, same as simulation
            lnr.eff.y   <- (lnr.y - min(lnr.y)) / (max(lnr.y) - min(lnr.y)) * (0.55 - 0.05) + 0.05 + 
              0.4 * rhc.used.sim$sim.A * (grepl("MOSF", rhc.used.sim$cat1) | (rhc.used.sim$cat1 == "ARF"))
            
            rhc.used.sim <- rhc.used.sim %>%
              mutate(prob.noeff.y = lnr.noeff.y,
                     prob.eff.y   = lnr.eff.y)
            
            rhc.used.sim <- rhc.used.sim %>%
              mutate(sim.Y.noeff = rbinom(nrow(rhc.used.sim), size = 1, prob = prob.noeff.y),
                     sim.Y.eff   = rbinom(nrow(rhc.used.sim), size = 1, prob = prob.eff.y))
            
            rhc.model.noeff <- data.frame(A = rhc.used.sim$sim.A,
                                          Y = rhc.used.sim$sim.Y.noeff,
                                          rhc.used.sim[, c("cat1", "imputed.cat2", "age", "surv2md1", "meanbp1", 
                                                           "resp1", "pafi1", "bili1", "dnr1")])
            
            rhc.model.eff <- data.frame(A = rhc.used.sim$sim.A,
                                        Y = rhc.used.sim$sim.Y.eff,
                                        rhc.used.sim[, c("cat1", "imputed.cat2", "age", "surv2md1", "meanbp1", 
                                                         "resp1", "pafi1", "bili1", "dnr1")])
            
            # Subsampling
            set.seed(a[i])
            
            ind.tree  <- sample(1:nrow(rhc.model.eff), ceiling(nrow(rhc.model.eff) * 0.5), replace = F)
            ind.train <- sample(ind.tree, round(length(ind.tree) * 0.8), replace = F)
            ind.val   <- ind.tree[!(ind.tree %in% ind.train)]
            ind.est   <- (1:nrow(rhc.model.eff))[!((1:nrow(rhc.model.eff)) %in% ind.tree)]
            
            rhc.eff.tree  <- rhc.model.eff[ind.tree, ]
            rhc.eff.train <- rhc.model.eff[ind.train, ]
            rhc.eff.val   <- rhc.model.eff[ind.val, ]
            rhc.eff.est   <- rhc.model.eff[ind.est, ]
            
            rhc.noeff.tree  <- rhc.model.noeff[ind.tree, ]
            rhc.noeff.train <- rhc.model.noeff[ind.train, ]
            rhc.noeff.val   <- rhc.model.noeff[ind.val, ]
            rhc.noeff.est   <- rhc.model.noeff[ind.est, ]
            
            res <- NULL
            
            ############################################################ 
            ################### Heterogeneous effect ###################
            ############################################################
            t0 <- Sys.time()
            seq.created.eff <- create.sequence(data.used         = rhc.eff.train,
                                               est.used          = "G",
                                               type.var          = "bin",
                                               adj.mod.loc       = "node",
                                               adj.mthd          = "GLM",
                                               adj.form.true     = NULL,
                                               num.truc.obs      = 30,
                                               min.node          = 20)
            t1 <- Sys.time()
            
            final.tree.eff <- EstG.CvMethod1(data.used         = rhc.eff.train,
                                             tree.list         = seq.created.eff$tree.list,
                                             lambda.used       = qchisq(0.95, 1),
                                             val.sample        = rhc.eff.val,              
                                             type.var          = "bin",
                                             adj.mod.loc       = "node",
                                             adj.mthd          = "GLM",
                                             adj.form.true     = NULL)
            
            t2 <- Sys.time()
            
            # time
            t.step.12 <- as.numeric(difftime(t1, t0, units = "secs"))
            t.step.3  <- as.numeric(difftime(t2, t1, units = "secs")) 
            
            # size of final tree
            size.tree <- sum(final.tree.eff[[1]]$frame$var == "<leaf>")
            
            size.tree.2ndLast <- sum(seq.created.eff$tree.list[[length(seq.created.eff$tree.list) - 1]]$frame$var == "<leaf>")
            
            res <- c(res, t.step.12, t.step.3, size.tree, size.tree.2ndLast)
            
            # treatment effect prediction error
            # correct tree
            # correct first split
            # number of noise variables
            trt.eff <- coef.rhc.y.select[names(coef.rhc.y.select) == "A"]
            eval.final.eff <- eval.measures.eff(final.tree   = final.tree.eff[[1]],
                                                test.data    = rhc.eff.est,
                                                true.trt.eff = trt.eff / (max(lnr.y) - min(lnr.y)) * (0.55 - 0.05) + 
                                                  0.4 * (grepl("MOSF", rhc.eff.est$cat1) | (rhc.eff.est$cat1 == "ARF")),
                                                noise.var    = colnames(rhc.eff.est)[!(colnames(rhc.eff.est) %in% c("A", "Y", "cat1"))],
                                                corr.split   = "cat1",
                                                where.split  = list(c(1)),
                                                dir.split    = list(c(NULL)),
                                                split.cate   = list(c(1, 8, 9)))
            
            res <- c(res, eval.final.eff$mse, eval.final.eff$exact.corr, eval.final.eff$numb.noise, eval.final.eff$pps)
            
            corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.eff$tree.list[[1]],
                                                       corr.split = "cat1",
                                                       split.cate = list(c(1, 8, 9)))
            res <- c(res, corr.frst.splt)
            
            # first split
            all.nodes <- rownames(seq.created.eff$tree.list[[1]]$frame)
            
            ind.node.1  <- which(all.nodes == 1)
            frst.splits <- as.character(seq.created.eff$tree.list[[1]]$frame$var[ind.node.1])
            ind.splits  <- which(seq.created.eff$tree.list[[1]]$splits[, 3] == seq.created.eff$tree.list[[1]]$frame$split.stat[ind.node.1])
            dir.split   <- seq.created.eff$tree.list[[1]]$splits[ind.splits, 2]
            ind.csplit  <- seq.created.eff$tree.list[[1]]$splits[ind.splits, 4]
            # sometimes csplit do not exist
            ncol.csplit <- ncol(seq.created.eff$tree.list[[1]]$csplit)
            # ncol.csplit <- 3
            
            # record the corresponding splitting point
            # need both 1) index of csplit be an integer and 2) ncat greater than or equal to 2
            if (((ind.csplit %% 1) == 0) & (dir.split >= 2)) {
              # placeholder for continuous split
              frst.splits <- c(frst.splits, rep(NA, 2))
              frst.splits <- c(frst.splits, seq.created.eff$tree.list[[1]]$csplit[ind.csplit, ])
            } else {
              frst.splits <- c(frst.splits, seq.created.eff$tree.list[[1]]$splits[ind.splits, c(2, 4)])
              frst.splits <- c(frst.splits, rep(NA, ncol.csplit))
            }
            
            for (split.i in 2:7) {
              
              if (split.i %in% all.nodes) {
                
                ind.node  <- which(all.nodes == split.i)
                split.var <- as.character(seq.created.eff$tree.list[[1]]$frame$var[ind.node])
                
                # if terminal node, put NA as placeholders
                if (split.var == "<leaf>") {
                  frst.splits <- c(frst.splits, rep(NA, 1 + 2 + ncol.csplit))
                  
                } else { # assign split information if not terminal node
                  frst.splits <- c(frst.splits, split.var)
                  ind.splits  <- which(seq.created.eff$tree.list[[1]]$splits[, 3] == seq.created.eff$tree.list[[1]]$frame$split.stat[ind.node])
                  dir.split   <- seq.created.eff$tree.list[[1]]$splits[ind.splits, 2]
                  ind.csplit  <- seq.created.eff$tree.list[[1]]$splits[ind.splits, 4]
                  
                  # record the corresponding splitting point
                  # need both 1) index of csplit be an integer and 2) ncat greater than or equal to 2
                  if (((ind.csplit %% 1) == 0) & (dir.split >= 2)) {
                    # placeholder for continuous split
                    frst.splits <- c(frst.splits, rep(NA, 2))
                    frst.splits <- c(frst.splits, seq.created.eff$tree.list[[1]]$csplit[ind.csplit, ])
                  } else {
                    frst.splits <- c(frst.splits, seq.created.eff$tree.list[[1]]$splits[ind.splits, c(2, 4)])
                    frst.splits <- c(frst.splits, rep(NA, ncol.csplit))
                  }
                  
                }
                
              } else {
                frst.splits <- c(frst.splits, rep(NA, 1 + 2 + ncol.csplit))
              }
              
            }
            
            res <- c(res, frst.splits)
            
            ############################################################ 
            #################### Homogeneous effect ####################
            ############################################################
            t0 <- Sys.time()
            seq.created.noeff <- create.sequence(data.used         = rhc.noeff.train,
                                                 est.used          = "G",
                                                 type.var          = "bin",
                                                 adj.mod.loc       = "node",
                                                 adj.mthd          = "GLM",
                                                 adj.form.true     = NULL,
                                                 num.truc.obs      = 30,
                                                 min.node          = 20)
            t1 <- Sys.time()
            
            final.tree.noeff <- EstG.CvMethod1(data.used         = rhc.noeff.train,
                                               tree.list         = seq.created.noeff$tree.list,
                                               lambda.used       = qchisq(0.95, 1),
                                               val.sample        = rhc.noeff.val,              
                                               type.var          = "bin",
                                               adj.mod.loc       = "node",
                                               adj.mthd          = "GLM",
                                               adj.form.true     = NULL)
            
            t2 <- Sys.time()
            
            # time
            t.step.12 <- as.numeric(difftime(t1, t0, units = "secs"))
            t.step.3  <- as.numeric(difftime(t2, t1, units = "secs")) 
            
            # size of final tree
            size.tree <- sum(final.tree.noeff[[1]]$frame$var == "<leaf>")
            
            size.tree.2ndLast <- sum(seq.created.noeff$tree.list[[length(seq.created.noeff$tree.list) - 1]]$frame$var == "<leaf>")
            
            res <- c(res, t.step.12, t.step.3, size.tree, size.tree.2ndLast)
            
            # treatment effect prediction error
            # correct tree
            # correct first split
            # number of noise variables
            trt.noeff <- coef.rhc.y.select[names(coef.rhc.y.select) == "A"]
            eval.final.noeff <- eval.measures.eff(final.tree   = final.tree.noeff[[1]],
                                                  test.data    = rhc.noeff.est,
                                                  true.trt.eff = rep(trt.noeff / (max(lnr.y) - min(lnr.y)) * (0.95 - 0.05), nrow(rhc.noeff.est)),
                                                  noise.var    = colnames(rhc.noeff.est)[!(colnames(rhc.noeff.est) %in% c("A", "Y"))],
                                                  corr.split   = NULL,
                                                  where.split  = NULL,
                                                  dir.split    = NULL,
                                                  split.cate   = NULL)
            
            res <- c(res, eval.final.noeff$mse, eval.final.noeff$exact.corr, eval.final.noeff$numb.noise, eval.final.noeff$pps)
            
            # first split
            all.nodes <- rownames(seq.created.noeff$tree.list[[1]]$frame)
            
            ind.node.1  <- which(all.nodes == 1)
            frst.splits <- as.character(seq.created.noeff$tree.list[[1]]$frame$var[ind.node.1])
            ind.splits  <- which(seq.created.noeff$tree.list[[1]]$splits[, 3] == seq.created.noeff$tree.list[[1]]$frame$split.stat[ind.node.1])
            dir.split   <- seq.created.noeff$tree.list[[1]]$splits[ind.splits, 2]
            ind.csplit  <- seq.created.noeff$tree.list[[1]]$splits[ind.splits, 4]
            # sometimes csplit do not exist
            ncol.csplit <- ncol(seq.created.noeff$tree.list[[1]]$csplit)
            # ncol.csplit <- 3
            
            # record the corresponding splitting point
            # need both 1) index of csplit be an integer and 2) ncat greater than or equal to 2
            if (((ind.csplit %% 1) == 0) & (dir.split >= 2)) {
              # placeholder for continuous split
              frst.splits <- c(frst.splits, rep(NA, 2))
              frst.splits <- c(frst.splits, seq.created.noeff$tree.list[[1]]$csplit[ind.csplit, ])
            } else {
              frst.splits <- c(frst.splits, seq.created.noeff$tree.list[[1]]$splits[ind.splits, c(2, 4)])
              frst.splits <- c(frst.splits, rep(NA, ncol.csplit))
            }
            
            for (split.i in 2:7) {
              
              if (split.i %in% all.nodes) {
                
                ind.node  <- which(all.nodes == split.i)
                split.var <- as.character(seq.created.noeff$tree.list[[1]]$frame$var[ind.node])
                
                # if terminal node, put NA as placeholders
                if (split.var == "<leaf>") {
                  frst.splits <- c(frst.splits, rep(NA, 1 + 2 + ncol.csplit))
                  
                } else { # assign split information if not terminal node
                  frst.splits <- c(frst.splits, split.var)
                  ind.splits  <- which(seq.created.noeff$tree.list[[1]]$splits[, 3] == seq.created.noeff$tree.list[[1]]$frame$split.stat[ind.node])
                  dir.split   <- seq.created.noeff$tree.list[[1]]$splits[ind.splits, 2]
                  ind.csplit  <- seq.created.noeff$tree.list[[1]]$splits[ind.splits, 4]
                  
                  # record the corresponding splitting point
                  # need both 1) index of csplit be an integer and 2) ncat greater than or equal to 2
                  if (((ind.csplit %% 1) == 0) & (dir.split >= 2)) {
                    # placeholder for continuous split
                    frst.splits <- c(frst.splits, rep(NA, 2))
                    frst.splits <- c(frst.splits, seq.created.noeff$tree.list[[1]]$csplit[ind.csplit, ])
                  } else {
                    frst.splits <- c(frst.splits, seq.created.noeff$tree.list[[1]]$splits[ind.splits, c(2, 4)])
                    frst.splits <- c(frst.splits, rep(NA, ncol.csplit))
                  }
                  
                }
                
              } else {
                frst.splits <- c(frst.splits, rep(NA, 1 + 2 + ncol.csplit))
              }
              
            }
            
            res <- c(res, frst.splits)
            
            # print must be put before output, otherwise the output will be print output
            if (i%%20 == 0) {print(i)}
            
            res
            
          }

save(res.rhcSimSubFtr.g.AFrst5.YLnrFrst5, file = paste0("../Data/AppendixD4/RhcSimSubFtrGAFrst5YLnrFrst5.RData"))
