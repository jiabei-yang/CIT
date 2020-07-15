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

boot.R <- 10^3
adj.mthd    <- "GLM"
type.var    <- "bin"

ci.g <- 
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
            hetero.trt.eff <- matrix(NA, nrow = boot.R, ncol = 12 / 4 * 2)
            
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
              # True
              hetero.trt.eff[boot.i, 1] <- mean(est.cond.eff.T.1.l - est.cond.eff.T.0.l)
              hetero.trt.eff[boot.i, 2] <- mean(est.cond.eff.T.1.r - est.cond.eff.T.0.r)
              
              # Noisy
              hetero.trt.eff[boot.i, 3] <- mean(est.cond.eff.Nois.1.l - est.cond.eff.Nois.0.l)
              hetero.trt.eff[boot.i, 4] <- mean(est.cond.eff.Nois.1.r - est.cond.eff.Nois.0.r)
              
              # Mis 
              hetero.trt.eff[boot.i, 5] <- mean(est.cond.eff.F.1.l - est.cond.eff.F.0.l)
              hetero.trt.eff[boot.i, 6] <- mean(est.cond.eff.F.1.r - est.cond.eff.F.0.r)
              
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
            homo.trt.eff <- matrix(NA, nrow = boot.R, ncol = 12 / 4)
            
            # Calculate the treatment effect and CI later anyway since might be used for a correct tree in either final tree selection method
            # will only calculate coverage using those iterations that get correct trees
            for (boot.i in 1:boot.R) {
              
              # seet seed for each bootstrap sample so can be replicated for each parameter setting
              set.seed(a[boot.i])
              
              data.boot <- data.used.full.bin.mixed[sample.int(1000, size = 1000, replace = T), ]
              data.boot.mis <- data.boot %>%
                select(-X2)
              
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
              #################### Find treatment effects in two root node ###################
              ################################################################################
              # True
              homo.trt.eff[boot.i, 1] <- mean(whl.g.T$value$pred.A.1 - whl.g.T$value$pred.A.0)
              
              # Noisy
              homo.trt.eff[boot.i, 2] <- mean(whl.g.Nois$value$pred.A.1 - whl.g.Nois$value$pred.A.0)
              
              # Mis 
              homo.trt.eff[boot.i, 3] <- mean(whl.g.F$value$pred.A.1 - whl.g.F$value$pred.A.0)
              
            } # for (boot.i in 1:boot.R) {
            
            res <- c(res, as.vector(apply(homo.trt.eff, 2, quantile, probs = c(0.025, 0.5, 0.975))))
            
            # print must be put before output, otherwise the output will be print output
            if (i%%30 == 0) {print(i)}
            
            res
            
          }

save(ci.g, file = paste0("../Data/AppendixC7/BinMixedGCover.RData"))


