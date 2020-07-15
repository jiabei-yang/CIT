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

rhc.model <- data.frame(A = rhc.used$num.swang1,
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
nrow.rhc.model <- nrow(rhc.model)

load("../Data/main/RhcOptDr.RData")

load("../Data/AppendixD2/RhcBootDrSplit.RData")
res.boot.dr <- rhc.boot.dr.result

res.boot.dr <- data.frame(res.boot.dr[, 1:15])
colnames(res.boot.dr)[1:6] <- c("t1", "t2", "size_tree", "frst_splt", "drct",
                                "splt_pt")
res.boot.dr <- res.boot.dr %>%
  mutate(drct    = as.numeric(levels(drct))[drct],
         splt_pt = as.numeric(levels(splt_pt))[splt_pt])

# sort(table(res.boot.dr$frst_splt))
ind.surv2md1         <- which(res.boot.dr$frst_splt == "surv2md1")
res.boot.dr.surv2md1 <- res.boot.dr[ind.surv2md1, ]
# hist(res.boot.dr.surv2md1$splt_pt, breaks = 30, xlim = c(0, 1))
# table(res.boot.dr.surv2md1$splt_pt)

ind.surv2md1.75 <- which((res.boot.dr$frst_splt == "surv2md1") & (res.boot.dr$splt_pt > 0.75))
splt.pt <- as.numeric(names(table(res.boot.dr.surv2md1$splt_pt))[which.max(table(res.boot.dr.surv2md1$splt_pt))])
# save(ind.surv2md1, ind.surv2md1.75, splt.pt, file = "../Data/Revision/RhcSplitBoot/RhcFrstSpltIndSurv2md120200708.RData")

boot.R <- 10^3
propsc.mthd <- "GLM"
adj.mthd    <- "GLM"
type.var    <- "bin"

ci.frst.splt.cit <- 
  foreach(i = 1:length(ind.surv2md1),
          .combine  = "rbind") %dopar% {
            
            set.seed(a[ind.surv2md1[i]])
            
            ind.tree  <- sample(1:nrow.rhc.model, ceiling(nrow.rhc.model * 0.5), replace = F)
            ind.est   <- (1:nrow.rhc.model)[!((1:nrow.rhc.model) %in% ind.tree)]
            rhc.model.est <- rhc.model[ind.est, ]
            
            # whether or not this iteration splits on the correct split point
            res <- (ind.surv2md1[i] %in% ind.surv2md1.75)
            
            trt.eff <- NULL
            for (boot.i in 1:boot.R) {
              
              set.seed(a[boot.i])
              rhc.model.est.used <- rhc.model.est[sample.int(nrow(rhc.model.est), replace = T), ]
              
              # left and right node observations
              rhc.model.est.used.l <- rhc.model.est.used[rhc.model.est.used$surv2md1 < splt.pt, ]
              rhc.model.est.used.r <- rhc.model.est.used[rhc.model.est.used$surv2md1 >= splt.pt, ]
              
              ################################################################################
              ############## fit the propensity score model on the whole dataset #############
              ################################################################################
              rhc.model.est.used.noy <- rhc.model.est.used[, !colnames(rhc.model.est.used) %in% c("Y")]
              
              whl.propsc <- gen.fullrank.ipw(df.noy           = rhc.model.est.used.noy, 
                                             propsc.form.true = NULL)
              whl.propsc <- withWarnings(est.prop.sc(df.noy    = whl.propsc$df.noy.fullrank,
                                                     method    = propsc.mthd,
                                                     form.true = whl.propsc$propsc.form.true.updated))
              
              ################################################################################
              ################## fit the outcome model on the whole dataset ##################
              ################################################################################
              whl.g <- gen.fullrank.g(df            = rhc.model.est.used,
                                      adj.form.true = NULL)
              whl.g <- withWarnings(est.cond.eff(df        = whl.g$df.fullrank,
                                                 method    = adj.mthd,
                                                 form.true = whl.g$adj.form.true.updated,
                                                 type.var  = type.var))
              
              ################################################################################
              ############## fit the propensity score model in two terminal nodes ############
              ################################################################################
              rhc.model.est.used.l.noy     <- rhc.model.est.used.l[, !colnames(rhc.model.est.used.l) %in% c("Y")]
              
              left.propsc <- gen.fullrank.ipw(df.noy           = rhc.model.est.used.l.noy, 
                                              propsc.form.true = NULL)
              left.propsc <- withWarnings(est.prop.sc(df.noy    = left.propsc$df.noy.fullrank,
                                                      method    = propsc.mthd,
                                                      form.true = left.propsc$propsc.form.true.updated))
              
              rhc.model.est.used.r.noy     <- rhc.model.est.used.r[, !colnames(rhc.model.est.used.r) %in% c("Y")]
              
              right.propsc <- gen.fullrank.ipw(df.noy           = rhc.model.est.used.r.noy, 
                                               propsc.form.true = NULL)
              right.propsc <- withWarnings(est.prop.sc(df.noy    = right.propsc$df.noy.fullrank,
                                                       method    = propsc.mthd,
                                                       form.true = right.propsc$propsc.form.true.updated))
              
              ################################################################################
              ################## fit the outcome model in two terminal nodes #################
              ################################################################################
              left.g <- gen.fullrank.g(df            = rhc.model.est.used.l,
                                       adj.form.true = NULL)
              left.g <- withWarnings(est.cond.eff(df        = left.g$df.fullrank,
                                                  method    = adj.mthd,
                                                  form.true = left.g$adj.form.true.updated,
                                                  type.var  = type.var))
              
              right.g <- gen.fullrank.g(df            = rhc.model.est.used.r,
                                        adj.form.true = NULL)
              right.g <- withWarnings(est.cond.eff(df        = right.g$df.fullrank,
                                                   method    = adj.mthd,
                                                   form.true = right.g$adj.form.true.updated,
                                                   type.var  = type.var))
              
              ################################################################################
              ################## find propensity scores in two terminal nodes ################
              ################################################################################
              # if 1) propensity score model fitting produces warning or 2) only produce thresholded propensity scores
              # use outside model fitting results
              if (!is.null(left.propsc$warnings)) {
                
                cond.warnings <- T
                for (warnings.i in 1:length(left.propsc$warnings)) {
                  cond.warnings <- cond.warnings & (left.propsc$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                # when there is only rank-deficient warning, use the local prediction
                if (cond.warnings) {
                  propsc.l <- left.propsc$value$prop.sc
                } else {
                  propsc.l <- whl.propsc$value$prop.sc[rhc.model.est.used$surv2md1 < splt.pt]
                }
                
              } else if (identical(sort(unique(left.propsc$value$prop.sc)), c(0.1, 0.9))) {
                propsc.l <- whl.propsc$value$prop.sc[rhc.model.est.used$surv2md1 < splt.pt]
                
              } else {
                propsc.l <- left.propsc$value$prop.sc
              }
              
              if (!is.null(right.propsc$warnings)) {
                
                cond.warnings <- T
                for (warnings.i in 1:length(right.propsc$warnings)) {
                  cond.warnings <- cond.warnings & (right.propsc$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                # when there is only rank-deficient warning, use the local prediction
                if (cond.warnings) {
                  propsc.r <- right.propsc$value$prop.sc
                } else {
                  propsc.r <- whl.propsc$value$prop.sc[rhc.model.est.used$surv2md1 >= splt.pt]
                }
                
              } else if (identical(sort(unique(right.propsc$value$prop.sc)), c(0.1, 0.9))) {
                propsc.r <- whl.propsc$value$prop.sc[rhc.model.est.used$surv2md1 >= splt.pt]
                
              } else {
                propsc.r <- right.propsc$value$prop.sc
              }
              
              ################################################################################
              ################# find potential outcomes in two terminal nodes ################
              ################################################################################
              # when the outcome model fitting produces warning
              # use outside model fitting results for conditional means
              if (!is.null(left.g$warnings)) {
                
                # when there is only rank-deficient warning, use the local prediction
                cond.warnings <- T
                for (warnings.i in 1:length(left.g$warnings)) {
                  cond.warnings <- cond.warnings & (left.g$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                if (cond.warnings) {
                  est.cond.eff.0.l <- left.g$value$pred.A.0
                  est.cond.eff.1.l <- left.g$value$pred.A.1
                } else {
                  est.cond.eff.0.l <- whl.g$value$pred.A.0[rhc.model.est.used$surv2md1 < splt.pt]
                  est.cond.eff.1.l <- whl.g$value$pred.A.1[rhc.model.est.used$surv2md1 < splt.pt]
                }
                
              } else { # if (!is.null(left.g$warnings)) {
                est.cond.eff.0.l <- left.g$value$pred.A.0
                est.cond.eff.1.l <- left.g$value$pred.A.1
              }
              
              if (!is.null(right.g$warnings)) {
                
                # when there is only rank-deficient warning, use the local prediction
                cond.warnings <- T
                for (warnings.i in 1:length(right.g$warnings)) {
                  cond.warnings <- cond.warnings & (right.g$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                if (cond.warnings) {
                  est.cond.eff.0.r <- right.g$value$pred.A.0
                  est.cond.eff.1.r <- right.g$value$pred.A.1
                } else {
                  est.cond.eff.0.r <- whl.g$value$pred.A.0[rhc.model.est.used$surv2md1 >= splt.pt]
                  est.cond.eff.1.r <- whl.g$value$pred.A.1[rhc.model.est.used$surv2md1 >= splt.pt]
                }
                
              } else { # if (!is.null(right.g$warnings)) {
                est.cond.eff.0.r <- right.g$value$pred.A.0
                est.cond.eff.1.r <- right.g$value$pred.A.1
              }
              
              ################################################################################
              ################## Find treatment effects in two terminal nodes ################
              ################################################################################
              ipw.trt.eff.l <- mean(rhc.model.est.used.l$Y * rhc.model.est.used.l$A / propsc.l) - 
                mean(rhc.model.est.used.l$Y * (1 - rhc.model.est.used.l$A) / (1 - propsc.l))
              ipw.trt.eff.r <- mean(rhc.model.est.used.r$Y * rhc.model.est.used.r$A / propsc.r) - 
                mean(rhc.model.est.used.r$Y * (1 - rhc.model.est.used.r$A) / (1 - propsc.r))
              
              g.trt.eff.l <- mean(est.cond.eff.1.l - est.cond.eff.0.l)
              g.trt.eff.r <- mean(est.cond.eff.1.r - est.cond.eff.0.r)
              
              dr.trt.eff.l <- mean(rhc.model.est.used.l$A * (rhc.model.est.used.l$Y - est.cond.eff.1.l) / propsc.l + est.cond.eff.1.l) - 
                mean((1 - rhc.model.est.used.l$A) * (rhc.model.est.used.l$Y - est.cond.eff.0.l) / (1 - propsc.l) + est.cond.eff.0.l)
              dr.trt.eff.r <- mean(rhc.model.est.used.r$A * (rhc.model.est.used.r$Y - est.cond.eff.1.r) / propsc.r + est.cond.eff.1.r) - 
                mean((1 - rhc.model.est.used.r$A) * (rhc.model.est.used.r$Y - est.cond.eff.0.r) / (1 - propsc.r) + est.cond.eff.0.r)
              
              trt.eff <- rbind(trt.eff, c(ipw.trt.eff.l, ipw.trt.eff.r, g.trt.eff.l, g.trt.eff.r, dr.trt.eff.l, dr.trt.eff.r))
              
              if ((boot.i %% 100) == 0) {print(i)}
              
            } # for boot.i loop
            
            res <- c(res, as.vector(apply(trt.eff, 2, quantile, probs = c(0.025, 0.5, 0.975))))
            
            res
            
          }

save(ci.frst.splt.cit, file = "../Data/AppendixD3/RhcSpltFrstSpltCi.RData")
