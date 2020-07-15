#!/usr/bin/env Rscript
setwd("../")
folder <- paste(getwd(), "/Functions/", sep="")
functions <- list.files(folder)
functions <- grep("main.R", functions, invert = T, value = T)
functions <- paste(folder, functions, sep = "")
for (i in functions){
  source(i)
}

setwd("Appendix/")

registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))

# Read in the arguments from command line
option_list = list(
  make_option(c("-s", "--start"), type="integer"),
  make_option(c("-e", "--end"), type = "integer"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
start <- opt$start
end   <- opt$end

load("../seed1000.rda")
# load("../Data/Revision/BinMixed/BinMixedDr20200524_1_1000.RData")
# # dim(performance.drInnd)
# performance.dr <- performance.drInnd
# performance.dr <- data.frame(rownames = rownames(performance.dr),
#                              performance.dr)
# 
# # create NA's in the dataframe
# all.rownames <- data.frame(rownames = paste0("result.", 1:1000))
# performance.dr <- all.rownames %>% 
#   left_join(performance.dr, by = "rownames")
# performance.dr <- performance.dr[, -1]
# load("../Data/Revision/BinMixed/BinMixedDrMissIters20200530.RData")
# # dim(performance.drInnd)
# performance.dr[is.na(performance.dr$hetero.dr.adjTGlmInnd.propscTGlmInnd.cv1.mse), ] <- performance.drInnd
# # dim(performance.dr)
# performance.dr <- as.matrix(performance.dr)
# 
# performance.dr <- rbind(performance.dr[, 8 * 0 + 1:8],
#                         performance.dr[, 8 * 1 + 1:8],
#                         performance.dr[, 8 * 2 + 1:8],
#                         performance.dr[, 8 * 3 + 1:8],
#                         performance.dr[, 8 * 4 + 1:8],
#                         performance.dr[, 8 * 5 + 1:8],
#                         performance.dr[, 8 * 6 + 1:8],
#                         performance.dr[, 8 * 7 + 1:8],
#                         performance.dr[, 8 * 8 + 1:8],
#                         performance.dr[, 8 * 9 + 1:8],
#                         cbind(performance.dr[, 8 * 10 + 7 * 0 + 1:7], NA),
#                         cbind(performance.dr[, 8 * 10 + 7 * 1 + 1:7], NA),
#                         cbind(performance.dr[, 8 * 10 + 7 * 2 + 1:7], NA),
#                         cbind(performance.dr[, 8 * 10 + 7 * 3 + 1:7], NA),
#                         cbind(performance.dr[, 8 * 10 + 7 * 4 + 1:7], NA),
#                         cbind(performance.dr[, 8 * 10 + 7 * 5 + 1:7], NA),
#                         cbind(performance.dr[, 8 * 10 + 7 * 6 + 1:7], NA),
#                         cbind(performance.dr[, 8 * 10 + 7 * 7 + 1:7], NA),
#                         cbind(performance.dr[, 8 * 10 + 7 * 8 + 1:7], NA),
#                         cbind(performance.dr[, 8 * 10 + 7 * 9 + 1:7], NA))
# rm(performance.drInnd)
# colnames(performance.dr)
# colnames(performance.dr) <- gsub("hetero.dr.adjTGlmInnd.propscTGlmInnd.cv1.", "", colnames(performance.dr))
# 
# algorithm <- c("Main FTS", "Alternative FTS")
# setting   <- c("Heterogeneous", "Homogeneous")
# Method   <- c("Both True DR-CIT", "True Out Mis Func Prop DR-CIT", "True Prop Mis Func Out DR-CIT", "Both Mis Func DR-CIT", "Both Unmeasured Cov DR-CIT")
# scnrs.cit <- expand.grid(Method, setting, algorithm)
# scnrs.cit <- scnrs.cit[c(1:5, 1:5 + 10, 1:5 + 5, 1:5 + 15), ]
# colnames(scnrs.cit) <- c("Method", "setting", "algorithm")
# 
# # scnrs <- rbind(scnrs.ct, scnrs.cit)
# performance.dr <- cbind(scnrs.cit[rep(1:nrow(scnrs.cit), each = 10^3), ],
#                         performance.dr)

boot.R <- 10^3
propsc.mthd <- "GLM"
adj.mthd    <- "GLM"
type.var    <- "bin"

ci.dr <- 
  foreach(i = start:end,
          .combine  = "rbind") %dopar% {
            
            # newly generated data should have seed different from before
            set.seed(a[i + 1000])  
            
            data.bin.mixed            <- makeData.bin.eff.mixed(N             = 1000, 
                                                                n.test        = 1000, 
                                                                p.cont        = 3, 
                                                                p.cate        = 3, 
                                                                n.cate        = 4:6, 
                                                                coeff.prop.sc = 0.3,
                                                                seed          = a[i + 1000])
            data.used.full.bin.mixed  <- data.bin.mixed$data.used
            data.used.bin.mixed       <- data.used.full.bin.mixed[1:800, ]
            # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
            data.validation.bin.mixed <- data.used.full.bin.mixed[801:1000, ]  
            
            # / 4 because only hetero and both final tree selection method should have the same treatment effect estimate
            # * 2 because two terminal nodes
            hetero.trt.eff <- matrix(NA, nrow = boot.R, ncol = 20 / 4 * 2)
            
            # Calculate the treatment effect and CI later anyway since might be used for a correct tree in either final tree selection method
            # will only calculate coverage using those iterations that get correct trees
            for (boot.i in 1:boot.R) {
              
              # seet seed for each bootstrap sample so can be replicated for each parameter setting
              set.seed(a[boot.i])
              data.boot <- data.used.full.bin.mixed[sample.int(1000, size = 1000, replace = T), ]
              data.boot.mis <- data.boot %>%
                select(-X2)
              
              data.l <- data.boot %>%
                filter(X4 %in% c("A", "C"))
              data.l.mis <- data.boot.mis %>%
                filter(X4 %in% c("A", "C"))
              data.r <- data.boot %>%
                filter(X4 %in% c("B", "D"))
              data.r.mis <- data.boot.mis %>%
                filter(X4 %in% c("B", "D"))
              
              ################################################################################
              ############## fit the propensity score model on the whole dataset #############
              ################################################################################
              data.boot.noy     <- data.boot[, !colnames(data.boot) %in% c("Y")]
              data.boot.mis.noy <- data.boot.mis[, !colnames(data.boot.mis) %in% c("Y")]
              
              whl.propsc.T <- gen.fullrank.ipw(df.noy           = data.boot.noy, 
                                               propsc.form.true = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))")
              whl.propsc.T <- withWarnings(est.prop.sc(df.noy    = whl.propsc.T$df.noy.fullrank,
                                                       method    = propsc.mthd,
                                                       form.true = whl.propsc.T$propsc.form.true.updated))
              
              whl.propsc.Nois <- gen.fullrank.ipw(df.noy           = data.boot.noy, 
                                                  propsc.form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6")
              whl.propsc.Nois <- withWarnings(est.prop.sc(df.noy    = whl.propsc.Nois$df.noy.fullrank,
                                                          method    = propsc.mthd,
                                                          form.true = whl.propsc.Nois$propsc.form.true.updated))
              
              whl.propsc.F <- gen.fullrank.ipw(df.noy           = data.boot.mis.noy, 
                                               propsc.form.true = NULL)
              whl.propsc.F <- withWarnings(est.prop.sc(df.noy    = whl.propsc.F$df.noy.fullrank,
                                                       method    = propsc.mthd,
                                                       form.true = whl.propsc.F$propsc.form.true.updated))
              
              ################################################################################
              ################## fit the outcome model on the whole dataset ##################
              ################################################################################
              whl.g.T <- gen.fullrank.g(df            = data.boot,
                                        adj.form.true = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))")
              whl.g.T <- withWarnings(est.cond.eff(df        = whl.g.T$df.fullrank,
                                                   method    = adj.mthd,
                                                   form.true = whl.g.T$adj.form.true.updated,
                                                   type.var  = type.var))
              
              whl.g.Nois <- gen.fullrank.g(df            = data.boot,
                                           adj.form.true = NULL)
              whl.g.Nois <- withWarnings(est.cond.eff(df        = whl.g.Nois$df.fullrank,
                                                      method    = adj.mthd,
                                                      form.true = whl.g.Nois$adj.form.true.updated,
                                                      type.var  = type.var))
              
              whl.g.F <- gen.fullrank.g(df            = data.boot.mis,
                                        adj.form.true = NULL)
              whl.g.F <- withWarnings(est.cond.eff(df        = whl.g.F$df.fullrank,
                                                   method    = adj.mthd,
                                                   form.true = whl.g.F$adj.form.true.updated,
                                                   type.var  = type.var))
              
              ################################################################################
              ############## fit the propensity score model in left terminal node ############
              ################################################################################
              data.l.noy     <- data.l[, !colnames(data.l) %in% c("Y")]
              data.l.mis.noy <- data.l.mis[, !colnames(data.l.mis) %in% c("Y")]
              
              left.propsc.T <- gen.fullrank.ipw(df.noy           = data.l.noy, 
                                                propsc.form.true = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))")
              left.propsc.T <- withWarnings(est.prop.sc(df.noy    = left.propsc.T$df.noy.fullrank,
                                                        method    = propsc.mthd,
                                                        form.true = left.propsc.T$propsc.form.true.updated))
              
              left.propsc.Nois <- gen.fullrank.ipw(df.noy           = data.l.noy, 
                                                   propsc.form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6")
              left.propsc.Nois <- withWarnings(est.prop.sc(df.noy    = left.propsc.Nois$df.noy.fullrank,
                                                           method    = propsc.mthd,
                                                           form.true = left.propsc.Nois$propsc.form.true.updated))
              
              left.propsc.F <- gen.fullrank.ipw(df.noy           = data.l.mis.noy, 
                                                propsc.form.true = NULL)
              left.propsc.F <- withWarnings(est.prop.sc(df.noy    = left.propsc.F$df.noy.fullrank,
                                                        method    = propsc.mthd,
                                                        form.true = left.propsc.F$propsc.form.true.updated))
              
              ################################################################################
              ################## fit the outcome model in left terminal node #################
              ################################################################################
              left.g.T <- gen.fullrank.g(df            = data.l,
                                         adj.form.true = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))")
              left.g.T <- withWarnings(est.cond.eff(df        = left.g.T$df.fullrank,
                                                    method    = adj.mthd,
                                                    form.true = left.g.T$adj.form.true.updated,
                                                    type.var  = type.var))
              
              left.g.Nois <- gen.fullrank.g(df            = data.l,
                                            adj.form.true = NULL)
              left.g.Nois <- withWarnings(est.cond.eff(df        = left.g.Nois$df.fullrank,
                                                       method    = adj.mthd,
                                                       form.true = left.g.Nois$adj.form.true.updated,
                                                       type.var  = type.var))
              
              left.g.F <- gen.fullrank.g(df            = data.l.mis,
                                         adj.form.true = NULL)
              left.g.F <- withWarnings(est.cond.eff(df        = left.g.F$df.fullrank,
                                                    method    = adj.mthd,
                                                    form.true = left.g.F$adj.form.true.updated,
                                                    type.var  = type.var))
              
              ################################################################################
              ############## fit the propensity score model in right terminal node ###########
              ################################################################################
              data.r.noy     <- data.r[, !colnames(data.r) %in% c("Y")]
              data.r.mis.noy <- data.r.mis[, !colnames(data.r.mis) %in% c("Y")]
              
              right.propsc.T <- gen.fullrank.ipw(df.noy           = data.r.noy, 
                                                 propsc.form.true = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))")
              right.propsc.T <- withWarnings(est.prop.sc(df.noy    = right.propsc.T$df.noy.fullrank,
                                                         method    = propsc.mthd,
                                                         form.true = right.propsc.T$propsc.form.true.updated))
              
              right.propsc.Nois <- gen.fullrank.ipw(df.noy           = data.r.noy, 
                                                    propsc.form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6")
              right.propsc.Nois <- withWarnings(est.prop.sc(df.noy    = right.propsc.Nois$df.noy.fullrank,
                                                            method    = propsc.mthd,
                                                            form.true = right.propsc.Nois$propsc.form.true.updated))
              
              right.propsc.F <- gen.fullrank.ipw(df.noy           = data.r.mis.noy, 
                                                 propsc.form.true = NULL)
              right.propsc.F <- withWarnings(est.prop.sc(df.noy    = right.propsc.F$df.noy.fullrank,
                                                         method    = propsc.mthd,
                                                         form.true = right.propsc.F$propsc.form.true.updated))
              
              ################################################################################
              ################## fit the outcome model in right terminal node ################
              ################################################################################
              right.g.T <- gen.fullrank.g(df            = data.r,
                                          adj.form.true = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))")
              right.g.T <- withWarnings(est.cond.eff(df        = right.g.T$df.fullrank,
                                                     method    = adj.mthd,
                                                     form.true = right.g.T$adj.form.true.updated,
                                                     type.var  = type.var))
              
              right.g.Nois <- gen.fullrank.g(df            = data.r,
                                             adj.form.true = NULL)
              right.g.Nois <- withWarnings(est.cond.eff(df        = right.g.Nois$df.fullrank,
                                                        method    = adj.mthd,
                                                        form.true = right.g.Nois$adj.form.true.updated,
                                                        type.var  = type.var))
              
              right.g.F <- gen.fullrank.g(df            = data.r.mis,
                                          adj.form.true = NULL)
              right.g.F <- withWarnings(est.cond.eff(df        = right.g.F$df.fullrank,
                                                     method    = adj.mthd,
                                                     form.true = right.g.F$adj.form.true.updated,
                                                     type.var  = type.var))
              
              ################################################################################
              ################## find propensity scores in left terminal node ################
              ################################################################################
              # if 1) propensity score model fitting produces warning or 2) only produce thresholded propensity scores
              # use outside model fitting results
              if (!is.null(left.propsc.T$warnings)) {
                
                cond.warnings <- T
                for (warnings.i in 1:length(left.propsc.T$warnings)) {
                  cond.warnings <- cond.warnings & (left.propsc.T$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                # when there is only rank-deficient warning, use the local prediction
                if (cond.warnings) {
                  propsc.T.l <- left.propsc.T$value$prop.sc
                } else {
                  propsc.T.l <- whl.propsc.T$value$prop.sc[data.boot$X4 %in% c("A", "C")]
                }
                
              } else if (identical(sort(unique(left.propsc.T$value$prop.sc)), c(0.1, 0.9))) {
                propsc.T.l <- whl.propsc.T$value$prop.sc[data.boot$X4 %in% c("A", "C")]
                
              } else {
                propsc.T.l <- left.propsc.T$value$prop.sc
              }
              
              if (!is.null(left.propsc.Nois$warnings)) {
                
                cond.warnings <- T
                for (warnings.i in 1:length(left.propsc.Nois$warnings)) {
                  cond.warnings <- cond.warnings & (left.propsc.Nois$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                # when there is only rank-deficient warning, use the local prediction
                if (cond.warnings) {
                  propsc.Nois.l <- left.propsc.Nois$value$prop.sc
                } else {
                  propsc.Nois.l <- whl.propsc.Nois$value$prop.sc[data.boot$X4 %in% c("A", "C")]
                }
                
              } else if (identical(sort(unique(left.propsc.Nois$value$prop.sc)), c(0.1, 0.9))) {
                propsc.Nois.l <- whl.propsc.Nois$value$prop.sc[data.boot$X4 %in% c("A", "C")]
                
              } else {
                propsc.Nois.l <- left.propsc.Nois$value$prop.sc
              }
              
              if (!is.null(left.propsc.F$warnings)) {
                
                cond.warnings <- T
                for (warnings.i in 1:length(left.propsc.F$warnings)) {
                  cond.warnings <- cond.warnings & (left.propsc.F$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                # when there is only rank-deficient warning, use the local prediction
                if (cond.warnings) {
                  propsc.F.l <- left.propsc.F$value$prop.sc
                } else {
                  propsc.F.l <- whl.propsc.F$value$prop.sc[data.boot$X4 %in% c("A", "C")]
                }
                
              } else if (identical(sort(unique(left.propsc.F$value$prop.sc)), c(0.1, 0.9))) {
                propsc.F.l <- whl.propsc.F$value$prop.sc[data.boot$X4 %in% c("A", "C")]
                
              } else {
                propsc.F.l <- left.propsc.F$value$prop.sc
              }
              
              ################################################################################
              ################# find potential outcomes in left terminal node ################
              ################################################################################
              # when the outcome model fitting produces warning
              # use outside model fitting results for conditional means
              if (!is.null(left.g.T$warnings)) {
                
                # when there is only rank-deficient warning, use the local prediction
                cond.warnings <- T
                for (warnings.i in 1:length(left.g.T$warnings)) {
                  cond.warnings <- cond.warnings & (left.g.T$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                if (cond.warnings) {
                  est.cond.eff.T.0.l <- left.g.T$value$pred.A.0
                  est.cond.eff.T.1.l <- left.g.T$value$pred.A.1
                } else {
                  est.cond.eff.T.0.l <- whl.g.T$value$pred.A.0[data.boot$X4 %in% c("A", "C")]
                  est.cond.eff.T.1.l <- whl.g.T$value$pred.A.1[data.boot$X4 %in% c("A", "C")]
                }
                
              } else { # if (!is.null(left.g.T$warnings)) {
                est.cond.eff.T.0.l <- left.g.T$value$pred.A.0
                est.cond.eff.T.1.l <- left.g.T$value$pred.A.1
              }
              
              if (!is.null(left.g.Nois$warnings)) {
                
                # when there is only rank-deficient warning, use the local prediction
                cond.warnings <- T
                for (warnings.i in 1:length(left.g.Nois$warnings)) {
                  cond.warnings <- cond.warnings & (left.g.Nois$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                if (cond.warnings) {
                  est.cond.eff.Nois.0.l <- left.g.Nois$value$pred.A.0
                  est.cond.eff.Nois.1.l <- left.g.Nois$value$pred.A.1
                } else {
                  est.cond.eff.Nois.0.l <- whl.g.Nois$value$pred.A.0[data.boot$X4 %in% c("A", "C")]
                  est.cond.eff.Nois.1.l <- whl.g.Nois$value$pred.A.1[data.boot$X4 %in% c("A", "C")]
                }
                
              } else { # if (!is.null(left.g.Nois$warnings)) {
                est.cond.eff.Nois.0.l <- left.g.Nois$value$pred.A.0
                est.cond.eff.Nois.1.l <- left.g.Nois$value$pred.A.1
              }
              
              if (!is.null(left.g.F$warnings)) {
                
                # when there is only rank-deficient warning, use the local prediction
                cond.warnings <- T
                for (warnings.i in 1:length(left.g.F$warnings)) {
                  cond.warnings <- cond.warnings & (left.g.F$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                if (cond.warnings) {
                  est.cond.eff.F.0.l <- left.g.F$value$pred.A.0
                  est.cond.eff.F.1.l <- left.g.F$value$pred.A.1
                } else {
                  est.cond.eff.F.0.l <- whl.g.F$value$pred.A.0[data.boot$X4 %in% c("A", "C")]
                  est.cond.eff.F.1.l <- whl.g.F$value$pred.A.1[data.boot$X4 %in% c("A", "C")]
                }
                
              } else { # if (!is.null(left.g.F$warnings)) {
                est.cond.eff.F.0.l <- left.g.F$value$pred.A.0
                est.cond.eff.F.1.l <- left.g.F$value$pred.A.1
              }
              
              ################################################################################
              ################# find propensity scores in right terminal node ################
              ################################################################################
              if (!is.null(right.propsc.T$warnings)) {
                
                cond.warnings <- T
                for (warnings.i in 1:length(right.propsc.T$warnings)) {
                  cond.warnings <- cond.warnings & (right.propsc.T$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                # when there is only rank-deficient warning, use the local prediction
                if (cond.warnings) {
                  propsc.T.r <- right.propsc.T$value$prop.sc
                } else {
                  propsc.T.r <- whl.propsc.T$value$prop.sc[data.boot$X4 %in% c("B", "D")]
                }
                
              } else if (identical(sort(unique(right.propsc.T$value$prop.sc)), c(0.1, 0.9))) {
                propsc.T.r <- whl.propsc.T$value$prop.sc[data.boot$X4 %in% c("B", "D")]
                
              } else {
                propsc.T.r <- right.propsc.T$value$prop.sc
              }
              
              if (!is.null(right.propsc.Nois$warnings)) {
                
                cond.warnings <- T
                for (warnings.i in 1:length(right.propsc.Nois$warnings)) {
                  cond.warnings <- cond.warnings & (right.propsc.Nois$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                # when there is only rank-deficient warning, use the local prediction
                if (cond.warnings) {
                  propsc.Nois.r <- right.propsc.Nois$value$prop.sc
                } else {
                  propsc.Nois.r <- whl.propsc.Nois$value$prop.sc[data.boot$X4 %in% c("B", "D")]
                }
                
              } else if (identical(sort(unique(right.propsc.Nois$value$prop.sc)), c(0.1, 0.9))) {
                propsc.Nois.r <- whl.propsc.Nois$value$prop.sc[data.boot$X4 %in% c("B", "D")]
                
              } else {
                propsc.Nois.r <- right.propsc.Nois$value$prop.sc
              }
              
              if (!is.null(right.propsc.F$warnings)) {
                
                cond.warnings <- T
                for (warnings.i in 1:length(right.propsc.F$warnings)) {
                  cond.warnings <- cond.warnings & (right.propsc.F$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                # when there is only rank-deficient warning, use the local prediction
                if (cond.warnings) {
                  propsc.F.r <- right.propsc.F$value$prop.sc
                } else {
                  propsc.F.r <- whl.propsc.F$value$prop.sc[data.boot$X4 %in% c("B", "D")]
                }
                
              } else if (identical(sort(unique(right.propsc.F$value$prop.sc)), c(0.1, 0.9))) {
                propsc.F.r <- whl.propsc.F$value$prop.sc[data.boot$X4 %in% c("B", "D")]
                
              } else {
                propsc.F.r <- right.propsc.F$value$prop.sc
              }
              
              ################################################################################
              ################ find potential outcomes in right terminal node ################
              ################################################################################
              # when the outcome model fitting produces warning
              # use outside model fitting results for conditional means
              if (!is.null(right.g.T$warnings)) {
                
                # when there is only rank-deficient warning, use the local prediction
                cond.warnings <- T
                for (warnings.i in 1:length(right.g.T$warnings)) {
                  cond.warnings <- cond.warnings & (right.g.T$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                if (cond.warnings) {
                  est.cond.eff.T.0.r <- right.g.T$value$pred.A.0
                  est.cond.eff.T.1.r <- right.g.T$value$pred.A.1
                } else {
                  est.cond.eff.T.0.r <- whl.g.T$value$pred.A.0[data.boot$X4 %in% c("B", "D")]
                  est.cond.eff.T.1.r <- whl.g.T$value$pred.A.1[data.boot$X4 %in% c("B", "D")]
                }
                
              } else { # if (!is.null(right.g.T$warnings)) {
                est.cond.eff.T.0.r <- right.g.T$value$pred.A.0
                est.cond.eff.T.1.r <- right.g.T$value$pred.A.1
              }
              
              if (!is.null(right.g.Nois$warnings)) {
                
                # when there is only rank-deficient warning, use the local prediction
                cond.warnings <- T
                for (warnings.i in 1:length(right.g.Nois$warnings)) {
                  cond.warnings <- cond.warnings & (right.g.Nois$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                if (cond.warnings) {
                  est.cond.eff.Nois.0.r <- right.g.Nois$value$pred.A.0
                  est.cond.eff.Nois.1.r <- right.g.Nois$value$pred.A.1
                } else {
                  est.cond.eff.Nois.0.r <- whl.g.Nois$value$pred.A.0[data.boot$X4 %in% c("B", "D")]
                  est.cond.eff.Nois.1.r <- whl.g.Nois$value$pred.A.1[data.boot$X4 %in% c("B", "D")]
                }
                
              } else { # if (!is.null(right.g.Nois$warnings)) {
                est.cond.eff.Nois.0.r <- right.g.Nois$value$pred.A.0
                est.cond.eff.Nois.1.r <- right.g.Nois$value$pred.A.1
              }
              
              if (!is.null(right.g.F$warnings)) {
                
                # when there is only rank-deficient warning, use the local prediction
                cond.warnings <- T
                for (warnings.i in 1:length(right.g.F$warnings)) {
                  cond.warnings <- cond.warnings & (right.g.F$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
                }
                
                if (cond.warnings) {
                  est.cond.eff.F.0.r <- right.g.F$value$pred.A.0
                  est.cond.eff.F.1.r <- right.g.F$value$pred.A.1
                } else {
                  est.cond.eff.F.0.r <- whl.g.F$value$pred.A.0[data.boot$X4 %in% c("B", "D")]
                  est.cond.eff.F.1.r <- whl.g.F$value$pred.A.1[data.boot$X4 %in% c("B", "D")]
                }
                
              } else { # if (!is.null(right.g.F$warnings)) {
                est.cond.eff.F.0.r <- right.g.F$value$pred.A.0
                est.cond.eff.F.1.r <- right.g.F$value$pred.A.1
              }
              
              ################################################################################
              ################## Find treatment effects in two terminal nodes ################
              ################################################################################
              # Both True
              hetero.trt.eff[boot.i, 1] <- mean(data.l$A * (data.l$Y - est.cond.eff.T.1.l) / propsc.T.l + est.cond.eff.T.1.l) - 
                mean((1 - data.l$A) * (data.l$Y - est.cond.eff.T.0.l) / (1 - propsc.T.l) + est.cond.eff.T.0.l)
              hetero.trt.eff[boot.i, 2] <- mean(data.r$A * (data.r$Y - est.cond.eff.T.1.r) / propsc.T.r + est.cond.eff.T.1.r) - 
                mean((1 - data.r$A) * (data.r$Y - est.cond.eff.T.0.r) / (1 - propsc.T.r) + est.cond.eff.T.0.r)
              
              # True outcome, Noisy propensity score
              hetero.trt.eff[boot.i, 3] <- mean(data.l$A * (data.l$Y - est.cond.eff.T.1.l) / propsc.Nois.l + est.cond.eff.T.1.l) - 
                mean((1 - data.l$A) * (data.l$Y - est.cond.eff.T.0.l) / (1 - propsc.Nois.l) + est.cond.eff.T.0.l)
              hetero.trt.eff[boot.i, 4] <- mean(data.r$A * (data.r$Y - est.cond.eff.T.1.r) / propsc.Nois.r + est.cond.eff.T.1.r) - 
                mean((1 - data.r$A) * (data.r$Y - est.cond.eff.T.0.r) / (1 - propsc.Nois.r) + est.cond.eff.T.0.r)
              
              # Noisy outcome, True propensity score
              hetero.trt.eff[boot.i, 5] <- mean(data.l$A * (data.l$Y - est.cond.eff.Nois.1.l) / propsc.T.l + est.cond.eff.Nois.1.l) - 
                mean((1 - data.l$A) * (data.l$Y - est.cond.eff.Nois.0.l) / (1 - propsc.T.l) + est.cond.eff.Nois.0.l)
              hetero.trt.eff[boot.i, 6] <- mean(data.r$A * (data.r$Y - est.cond.eff.Nois.1.r) / propsc.T.r + est.cond.eff.Nois.1.r) - 
                mean((1 - data.r$A) * (data.r$Y - est.cond.eff.Nois.0.r) / (1 - propsc.T.r) + est.cond.eff.Nois.0.r)
              
              # Noisy outcome, Noisy propensity score
              hetero.trt.eff[boot.i, 7] <- mean(data.l$A * (data.l$Y - est.cond.eff.Nois.1.l) / propsc.Nois.l + est.cond.eff.Nois.1.l) - 
                mean((1 - data.l$A) * (data.l$Y - est.cond.eff.Nois.0.l) / (1 - propsc.Nois.l) + est.cond.eff.Nois.0.l)
              hetero.trt.eff[boot.i, 8] <- mean(data.r$A * (data.r$Y - est.cond.eff.Nois.1.r) / propsc.Nois.r + est.cond.eff.Nois.1.r) - 
                mean((1 - data.r$A) * (data.r$Y - est.cond.eff.Nois.0.r) / (1 - propsc.Nois.r) + est.cond.eff.Nois.0.r)
              
              # Mis outcome, Mis propensity score
              hetero.trt.eff[boot.i, 9] <- mean(data.l$A * (data.l$Y - est.cond.eff.F.1.l) / propsc.F.l + est.cond.eff.F.1.l) - 
                mean((1 - data.l$A) * (data.l$Y - est.cond.eff.F.0.l) / (1 - propsc.F.l) + est.cond.eff.F.0.l)
              hetero.trt.eff[boot.i, 10] <- mean(data.r$A * (data.r$Y - est.cond.eff.F.1.r) / propsc.F.r + est.cond.eff.F.1.r) - 
                mean((1 - data.r$A) * (data.r$Y - est.cond.eff.F.0.r) / (1 - propsc.F.r) + est.cond.eff.F.0.r)
              
            } # for (boot.i in 1:boot.R) {
            
            res <- as.vector(apply(hetero.trt.eff, 2, quantile, probs = c(0.025, 0.5, 0.975)))
            
            #####################################################################################################################
            ################################################### Homogeneous #####################################################
            #####################################################################################################################
            # newly generated data should have seed different from before
            set.seed(a[i + 1000])  
            
            data.bin.mixed            <- makeData.bin.noeff.mixed(N             = 1000, 
                                                                  n.test        = 1000, 
                                                                  p.cont        = 3, 
                                                                  p.cate        = 3, 
                                                                  n.cate        = 4:6, 
                                                                  coeff.prop.sc = 0.3,
                                                                  seed          = a[i + 1000])
            data.used.full.bin.mixed  <- data.bin.mixed$data.used
            data.used.bin.mixed       <- data.used.full.bin.mixed[1:800, ]
            # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
            data.validation.bin.mixed <- data.used.full.bin.mixed[801:1000, ]  
            
            # / 4 because only hetero and both final tree selection method should have the same treatment effect estimate
            homo.trt.eff <- matrix(NA, nrow = boot.R, ncol = 20 / 4)
            
            # Calculate the treatment effect and CI later anyway since might be used for a correct tree in either final tree selection method
            # will only calculate coverage using those iterations that get correct trees
            for (boot.i in 1:boot.R) {
              
              # seet seed for each bootstrap sample so can be replicated for each parameter setting
              set.seed(a[boot.i])
              
              data.boot <- data.used.full.bin.mixed[sample.int(1000, size = 1000, replace = T), ]
              data.boot.mis <- data.boot %>%
                select(-X2)
              
              ################################################################################
              ############## fit the propensity score model on the whole dataset #############
              ################################################################################
              data.boot.noy     <- data.boot[, !colnames(data.boot) %in% c("Y")]
              data.boot.mis.noy <- data.boot.mis[, !colnames(data.boot.mis) %in% c("Y")]
              
              whl.propsc.T <- gen.fullrank.ipw(df.noy           = data.boot.noy, 
                                               propsc.form.true = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))")
              whl.propsc.T <- withWarnings(est.prop.sc(df.noy    = whl.propsc.T$df.noy.fullrank,
                                                       method    = propsc.mthd,
                                                       form.true = whl.propsc.T$propsc.form.true.updated))
              
              whl.propsc.Nois <- gen.fullrank.ipw(df.noy           = data.boot.noy, 
                                                  propsc.form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6")
              whl.propsc.Nois <- withWarnings(est.prop.sc(df.noy    = whl.propsc.Nois$df.noy.fullrank,
                                                          method    = propsc.mthd,
                                                          form.true = whl.propsc.Nois$propsc.form.true.updated))
              
              whl.propsc.F <- gen.fullrank.ipw(df.noy           = data.boot.mis.noy, 
                                               propsc.form.true = NULL)
              whl.propsc.F <- withWarnings(est.prop.sc(df.noy    = whl.propsc.F$df.noy.fullrank,
                                                       method    = propsc.mthd,
                                                       form.true = whl.propsc.F$propsc.form.true.updated))
              
              ################################################################################
              ################## fit the outcome model on the whole dataset ##################
              ################################################################################
              whl.g.T <- gen.fullrank.g(df            = data.boot,
                                        adj.form.true = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))")
              whl.g.T <- withWarnings(est.cond.eff(df        = whl.g.T$df.fullrank,
                                                   method    = adj.mthd,
                                                   form.true = whl.g.T$adj.form.true.updated,
                                                   type.var  = type.var))
              
              whl.g.Nois <- gen.fullrank.g(df            = data.boot,
                                           adj.form.true = NULL)
              whl.g.Nois <- withWarnings(est.cond.eff(df        = whl.g.Nois$df.fullrank,
                                                      method    = adj.mthd,
                                                      form.true = whl.g.Nois$adj.form.true.updated,
                                                      type.var  = type.var))
              
              whl.g.F <- gen.fullrank.g(df            = data.boot.mis,
                                        adj.form.true = NULL)
              whl.g.F <- withWarnings(est.cond.eff(df        = whl.g.F$df.fullrank,
                                                   method    = adj.mthd,
                                                   form.true = whl.g.F$adj.form.true.updated,
                                                   type.var  = type.var))
              
              ################################################################################
              ################## Find treatment effects in two terminal nodes ################
              ################################################################################
              # Both True
              homo.trt.eff[boot.i, 1] <- mean(data.boot$A * (data.boot$Y - whl.g.T$value$pred.A.1) / whl.propsc.T$value$prop.sc + whl.g.T$value$pred.A.1) - 
                mean((1 - data.boot$A) * (data.boot$Y - whl.g.T$value$pred.A.0) / (1 - whl.propsc.T$value$prop.sc) + whl.g.T$value$pred.A.0)
              
              # True outcome, Noisy propensity score
              homo.trt.eff[boot.i, 2] <- mean(data.boot$A * (data.boot$Y - whl.g.T$value$pred.A.1) / whl.propsc.Nois$value$prop.sc + whl.g.T$value$pred.A.1) - 
                mean((1 - data.boot$A) * (data.boot$Y - whl.g.T$value$pred.A.0) / (1 - whl.propsc.Nois$value$prop.sc) + whl.g.T$value$pred.A.0)
              
              # Noisy outcome, True propensity score
              homo.trt.eff[boot.i, 3] <- mean(data.boot$A * (data.boot$Y - whl.g.Nois$value$pred.A.1) / whl.propsc.T$value$prop.sc + whl.g.Nois$value$pred.A.1) - 
                mean((1 - data.boot$A) * (data.boot$Y - whl.g.Nois$value$pred.A.0) / (1 - whl.propsc.T$value$prop.sc) + whl.g.Nois$value$pred.A.0)
              
              # Noisy outcome, Noisy propensity score
              homo.trt.eff[boot.i, 4] <- mean(data.boot$A * (data.boot$Y - whl.g.Nois$value$pred.A.1) / whl.propsc.Nois$value$prop.sc + whl.g.Nois$value$pred.A.1) - 
                mean((1 - data.boot$A) * (data.boot$Y - whl.g.Nois$value$pred.A.0) / (1 - whl.propsc.Nois$value$prop.sc) + whl.g.Nois$value$pred.A.0)
              
              # Mis outcome, Mis propensity score
              homo.trt.eff[boot.i, 5] <- mean(data.boot$A * (data.boot$Y - whl.g.F$value$pred.A.1) / whl.propsc.F$value$prop.sc + whl.g.F$value$pred.A.1) - 
                mean((1 - data.boot$A) * (data.boot$Y - whl.g.F$value$pred.A.0) / (1 - whl.propsc.F$value$prop.sc) + whl.g.F$value$pred.A.0)
              
            } # for (boot.i in 1:boot.R) {
            
            
            res <- c(res, as.vector(apply(homo.trt.eff, 2, quantile, probs = c(0.025, 0.5, 0.975))))
            
            # print must be put before output, otherwise the output will be print output
            if (i%%30 == 0) {print(i)}
            
            res
            
          }

save(ci.dr, file = paste0("../Data/AppendixC7/BinMixedDrCover.RData"))



