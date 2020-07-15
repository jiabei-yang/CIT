#!/usr/bin/env Rscript
library(causalTree)

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

performance.ct <- foreach(i = start:end,
                          .combine  = "rbind") %dopar% {
                            
                            set.seed(a[i])
                            
                            #####################################################################################################################
                            ################################################# Heterogeneous #####################################################
                            #####################################################################################################################
                            data.bin.mixed            <- makeData.bin.eff.mixed(N             = 1000, 
                                                                                n.test        = 1000, 
                                                                                p.cont        = 3, 
                                                                                p.cate        = 3, 
                                                                                n.cate        = 4:6, 
                                                                                coeff.prop.sc = 0.3,
                                                                                seed          = a[i])
                            data.used.full.bin.mixed  <- data.bin.mixed$data.used
                            
                            trtIdx  <- which(data.used.full.bin.mixed$A == 1)
                            ctrlIdx <- which(data.used.full.bin.mixed$A == 0)
                            train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
                                           sample(ctrlIdx, length(ctrlIdx) / 2))
                            train.data <- data.used.full.bin.mixed[train.idx, ]
                            est.data   <- data.used.full.bin.mixed[-train.idx, ]
                            
                            #####################################################################################################################
                            #################################### 1. True propensity score model, no honest ######################################
                            #####################################################################################################################
                            # In help document: Unit-specific propensity scores are not supported
                            t0 <- Sys.time()
                            tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                                                      method    = "GLM",
                                                      form.true = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))")
                            tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
                            
                            ct.propsc.true.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
                                                                  data         = data.used.full.bin.mixed,
                                                                  weights      = 1 / tmp.propsc$prop.sc,
                                                                  treatment    = data.used.full.bin.mixed$A,
                                                                  split.Rule   = "CT", 
                                                                  split.Honest = T, 
                                                                  cv.option    = "CT", 
                                                                  cv.Honest    = T,
                                                                  split.Bucket = F, 
                                                                  xval         = 5, 
                                                                  cp           = 0, 
                                                                  minsize      = 20)
                            cptable.propsc.true.nohonest <- ct.propsc.true.nohonest$cptable[,1][which.min(ct.propsc.true.nohonest$cptable[,4])]
                            final.tree.propsc.true.nohonest <- prune(ct.propsc.true.nohonest, cptable.propsc.true.nohonest)
                            t1 <- Sys.time()
                            
                            eval.ct.propsc.true.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.true.nohonest,
                                                                              test.data    = data.bin.mixed$test.data,
                                                                              true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                              noise.var    = data.bin.mixed$noise.var,
                                                                              corr.split   = data.bin.mixed$corr.split,
                                                                              where.split  = data.bin.mixed$where.split,
                                                                              dir.split    = data.bin.mixed$dir.split,
                                                                              split.cate   = data.bin.mixed$split.cate,
                                                                              CT           = T)
                            eval.ct.propsc.true.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
                            eval.ct.propsc.true.nohonest$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = ct.propsc.true.nohonest,
                                                                                                    corr.split = data.bin.mixed$corr.split,
                                                                                                    split.cate = data.bin.mixed$split.cate)
                            eval.ct.propsc.true.nohonest                <- unlist(eval.ct.propsc.true.nohonest)
                            names(eval.ct.propsc.true.nohonest)         <- paste0("hetero.propsc.true.nohonest.", names(eval.ct.propsc.true.nohonest))
                            print("hetero.true.nohonest")   
                            
                            #####################################################################################################################
                            ################################## 2. True propensity score model, honest Estimation ################################
                            #####################################################################################################################
                            t0 <- Sys.time()
                            tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                                                      method    = "GLM",
                                                      form.true = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))")
                            tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
                            
                            ct.propsc.true.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                                                       data             = train.data,
                                                                       weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                                                       treatment        = train.data$A,
                                                                       est_data         = est.data,
                                                                       est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                                                       est_treatment    = est.data$A,
                                                                       split.Rule       = "CT", 
                                                                       split.Honest     = T,
                                                                       HonestSampleSize = nrow(est.data),
                                                                       split.Bucket     = F,
                                                                       cv.option        = "CT",
                                                                       cv.Honest        = T)
                            cptable.propsc.true.honest <- ct.propsc.true.honest$cptable[,1][which.min(ct.propsc.true.honest$cptable[,4])]
                            final.tree.propsc.true.honest <- prune(ct.propsc.true.honest, cptable.propsc.true.honest)
                            t1 <- Sys.time()
                            
                            eval.ct.propsc.true.honest <- eval.measures.eff(final.tree   = final.tree.propsc.true.honest,
                                                                            test.data    = data.bin.mixed$test.data,
                                                                            true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                            noise.var    = data.bin.mixed$noise.var,
                                                                            corr.split   = data.bin.mixed$corr.split,
                                                                            where.split  = data.bin.mixed$where.split,
                                                                            dir.split    = data.bin.mixed$dir.split,
                                                                            split.cate   = data.bin.mixed$split.cate,
                                                                            CT           = T)
                            eval.ct.propsc.true.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
                            eval.ct.propsc.true.honest$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = ct.propsc.true.honest,
                                                                                                  corr.split = data.bin.mixed$corr.split,
                                                                                                  split.cate = data.bin.mixed$split.cate)
                            eval.ct.propsc.true.honest                <- unlist(eval.ct.propsc.true.honest)
                            names(eval.ct.propsc.true.honest)         <- paste0("hetero.propsc.true.honest.", names(eval.ct.propsc.true.honest))
                            print("hetero.true.honest")  
                            
                            #####################################################################################################################
                            ##################################### 3. Noisy propensity score model, no honest ####################################
                            #####################################################################################################################
                            t0 <- Sys.time()
                            tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                                                      method    = "GLM",
                                                      form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6")
                            tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
                            
                            ct.propsc.nois.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
                                                                  data         = data.used.full.bin.mixed,
                                                                  weights      = 1 / tmp.propsc$prop.sc,
                                                                  treatment    = data.used.full.bin.mixed$A,
                                                                  split.Rule   = "CT", 
                                                                  cv.option    = "CT", 
                                                                  split.Honest = T, 
                                                                  split.Bucket = F,
                                                                  cv.Honest    = T,
                                                                  xval         = 5, 
                                                                  cp           = 0, 
                                                                  minsize      = 20)
                            
                            cptable.propsc.nois.nohonest <- ct.propsc.nois.nohonest$cptable[,1][which.min(ct.propsc.nois.nohonest$cptable[,4])]
                            final.tree.propsc.nois.nohonest <- prune(ct.propsc.nois.nohonest, cptable.propsc.nois.nohonest)
                            t1 <- Sys.time()
                            
                            eval.ct.propsc.nois.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.nois.nohonest,
                                                                              test.data    = data.bin.mixed$test.data,
                                                                              true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                              noise.var    = data.bin.mixed$noise.var,
                                                                              corr.split   = data.bin.mixed$corr.split,
                                                                              where.split  = data.bin.mixed$where.split,
                                                                              dir.split    = data.bin.mixed$dir.split,
                                                                              split.cate   = data.bin.mixed$split.cate,
                                                                              CT           = T)
                            eval.ct.propsc.nois.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
                            eval.ct.propsc.nois.nohonest$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = ct.propsc.nois.nohonest,
                                                                                                    corr.split = data.bin.mixed$corr.split,
                                                                                                    split.cate = data.bin.mixed$split.cate)
                            eval.ct.propsc.nois.nohonest                <- unlist(eval.ct.propsc.nois.nohonest)
                            names(eval.ct.propsc.nois.nohonest)         <- paste0("hetero.propsc.nois.nohonest.", names(eval.ct.propsc.nois.nohonest))
                            print("hetero.nois.nohonest")  
                            
                            #####################################################################################################################
                            ################################## 4. Noisy propensity score model, honest estimation ###############################
                            #####################################################################################################################
                            t0 <- Sys.time()
                            tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                                                      method    = "GLM",
                                                      form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6")
                            tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
                            
                            ct.propsc.nois.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                                                       data             = train.data,
                                                                       weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                                                       treatment        = train.data$A,
                                                                       est_data         = est.data,
                                                                       est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                                                       est_treatment    = est.data$A,
                                                                       split.Rule       = "CT", 
                                                                       split.Honest     = T,
                                                                       HonestSampleSize = nrow(est.data),
                                                                       split.Bucket     = F,
                                                                       cv.option        = "CT",
                                                                       cv.Honest        = T)
                            cptable.propsc.nois.honest <- ct.propsc.nois.honest$cptable[,1][which.min(ct.propsc.nois.honest$cptable[,4])]
                            final.tree.propsc.nois.honest <- prune(ct.propsc.nois.honest, cptable.propsc.nois.honest)
                            t1 <- Sys.time()
                            
                            eval.ct.propsc.nois.honest <- eval.measures.eff(final.tree   = final.tree.propsc.nois.honest,
                                                                            test.data    = data.bin.mixed$test.data,
                                                                            true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                            noise.var    = data.bin.mixed$noise.var,
                                                                            corr.split   = data.bin.mixed$corr.split,
                                                                            where.split  = data.bin.mixed$where.split,
                                                                            dir.split    = data.bin.mixed$dir.split,
                                                                            split.cate   = data.bin.mixed$split.cate,
                                                                            CT           = T)
                            eval.ct.propsc.nois.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
                            eval.ct.propsc.nois.honest$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = ct.propsc.nois.honest,
                                                                                                  corr.split = data.bin.mixed$corr.split,
                                                                                                  split.cate = data.bin.mixed$split.cate)
                            eval.ct.propsc.nois.honest                <- unlist(eval.ct.propsc.nois.honest)
                            names(eval.ct.propsc.nois.honest)         <- paste0("hetero.propsc.nois.honest.", names(eval.ct.propsc.nois.honest))
                            print("hetero.nois.honest")   
                            
                            #####################################################################################################################
                            ################################## 5. Misspecified propensity score model, no honest ################################
                            #####################################################################################################################
                            data.used.full.bin.mixed.mis <- data.used.full.bin.mixed %>%
                              select(-X2)
                            
                            t0 <- Sys.time()
                            tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed.mis[, !colnames(data.used.full.bin.mixed.mis) %in% c("Y")],
                                                      method    = "GLM",
                                                      form.true = NULL)
                            tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
                            
                            ct.propsc.mis.nohonest <- causalTree(Y ~ X1 + X3 + X4 + X5 + X6, 
                                                                 data         = data.used.full.bin.mixed.mis,
                                                                 weights      = 1 / tmp.propsc$prop.sc,
                                                                 treatment    = data.used.full.bin.mixed.mis$A,
                                                                 split.Rule   = "CT", 
                                                                 cv.option    = "CT", 
                                                                 split.Honest = T, 
                                                                 split.Bucket = F, 
                                                                 cv.Honest    = T,
                                                                 xval         = 5, 
                                                                 cp           = 0, 
                                                                 minsize      = 20)
                            
                            cptable.propsc.mis.nohonest <- ct.propsc.mis.nohonest$cptable[,1][which.min(ct.propsc.mis.nohonest$cptable[,4])]
                            final.tree.propsc.mis.nohonest <- prune(ct.propsc.mis.nohonest, cptable.propsc.mis.nohonest)
                            t1 <- Sys.time()
                            
                            eval.ct.propsc.mis.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.nohonest,
                                                                             test.data    = data.bin.mixed$test.data,
                                                                             true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                             noise.var    = data.bin.mixed$noise.var,
                                                                             corr.split   = data.bin.mixed$corr.split,
                                                                             where.split  = data.bin.mixed$where.split,
                                                                             dir.split    = data.bin.mixed$dir.split,
                                                                             split.cate   = data.bin.mixed$split.cate,
                                                                             CT           = T)
                            eval.ct.propsc.mis.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
                            eval.ct.propsc.mis.nohonest$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = ct.propsc.mis.nohonest,
                                                                                                   corr.split = data.bin.mixed$corr.split,
                                                                                                   split.cate = data.bin.mixed$split.cate)
                            eval.ct.propsc.mis.nohonest                <- unlist(eval.ct.propsc.mis.nohonest)
                            names(eval.ct.propsc.mis.nohonest)         <- paste0("hetero.propsc.mis.nohonest.", names(eval.ct.propsc.mis.nohonest))
                            print("hetero.mis.nohonest") 
                            
                            #####################################################################################################################
                            ############################ 6. Misspecified propensity score model, honest estimation ##############################
                            #####################################################################################################################
                            train.data.mis <- train.data %>%
                              select(-X2)
                            est.data.mis   <- est.data %>%
                              select(-X2)
                            
                            t0 <- Sys.time()
                            tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed.mis[, !colnames(data.used.full.bin.mixed.mis) %in% c("Y")],
                                                      method    = "GLM",
                                                      form.true = NULL)
                            tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
                            
                            ct.propsc.mis.honest <- honest.causalTree(Y ~ X1 + X3 + X4 + X5 + X6,
                                                                      data             = train.data.mis,
                                                                      weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                                                      treatment        = train.data.mis$A,
                                                                      est_data         = est.data.mis,
                                                                      est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                                                      est_treatment    = est.data.mis$A,
                                                                      split.Rule       = "CT", 
                                                                      split.Honest     = T,
                                                                      HonestSampleSize = nrow(est.data.mis),
                                                                      split.Bucket     = F,
                                                                      cv.option        = "CT",
                                                                      cv.Honest        = T)
                            cptable.propsc.mis.honest <- ct.propsc.mis.honest$cptable[,1][which.min(ct.propsc.mis.honest$cptable[,4])]
                            final.tree.propsc.mis.honest <- prune(ct.propsc.mis.honest, cptable.propsc.mis.honest)
                            t1 <- Sys.time()
                            
                            eval.ct.propsc.mis.honest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.honest,
                                                                           test.data    = data.bin.mixed$test.data,
                                                                           true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                           noise.var    = data.bin.mixed$noise.var,
                                                                           corr.split   = data.bin.mixed$corr.split,
                                                                           where.split  = data.bin.mixed$where.split,
                                                                           dir.split    = data.bin.mixed$dir.split,
                                                                           split.cate   = data.bin.mixed$split.cate,
                                                                           CT           = T)
                            eval.ct.propsc.mis.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
                            eval.ct.propsc.mis.honest$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = ct.propsc.mis.honest,
                                                                                                 corr.split = data.bin.mixed$corr.split,
                                                                                                 split.cate = data.bin.mixed$split.cate)
                            eval.ct.propsc.mis.honest                <- unlist(eval.ct.propsc.mis.honest)
                            names(eval.ct.propsc.mis.honest)         <- paste0("hetero.propsc.mis.honest.", names(eval.ct.propsc.mis.honest))
                            print("hetero.mis.nohonest")   
                            
                            performance.hetero.ct <- c(eval.ct.propsc.true.nohonest,
                                                       eval.ct.propsc.true.honest,
                                                       eval.ct.propsc.nois.nohonest,
                                                       eval.ct.propsc.nois.honest,
                                                       eval.ct.propsc.mis.nohonest,
                                                       eval.ct.propsc.mis.honest)
                            
                            #####################################################################################################################
                            ################################################### Homogeneous #####################################################
                            #####################################################################################################################
                            data.bin.mixed            <- makeData.bin.noeff.mixed(N             = 1000, 
                                                                                  n.test        = 1000, 
                                                                                  p.cont        = 3, 
                                                                                  p.cate        = 3, 
                                                                                  n.cate        = 4:6, 
                                                                                  coeff.prop.sc = 0.3,
                                                                                  seed          = a[i])
                            data.used.full.bin.mixed  <- data.bin.mixed$data.used
                            
                            trtIdx  <- which(data.used.full.bin.mixed$A == 1)
                            ctrlIdx <- which(data.used.full.bin.mixed$A == 0)
                            train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
                                           sample(ctrlIdx, length(ctrlIdx) / 2))
                            train.data <- data.used.full.bin.mixed[train.idx, ]
                            est.data   <- data.used.full.bin.mixed[-train.idx, ]
                            
                            #####################################################################################################################
                            #################################### 7. True propensity score model, no honest ######################################
                            #####################################################################################################################
                            # In help document: Unit-specific propensity scores are not supported
                            t0 <- Sys.time()
                            tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                                                      method    = "GLM",
                                                      form.true = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))")
                            tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
                            
                            ct.propsc.true.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
                                                                  data         = data.used.full.bin.mixed,
                                                                  weights      = 1 / tmp.propsc$prop.sc,
                                                                  treatment    = data.used.full.bin.mixed$A,
                                                                  split.Rule   = "CT", 
                                                                  split.Honest = T, 
                                                                  cv.option    = "CT", 
                                                                  cv.Honest    = T,
                                                                  split.Bucket = F, 
                                                                  xval         = 5, 
                                                                  cp           = 0, 
                                                                  minsize      = 20)
                            cptable.propsc.true.nohonest <- ct.propsc.true.nohonest$cptable[,1][which.min(ct.propsc.true.nohonest$cptable[,4])]
                            final.tree.propsc.true.nohonest <- prune(ct.propsc.true.nohonest, cptable.propsc.true.nohonest)
                            t1 <- Sys.time()
                            
                            eval.ct.propsc.true.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.true.nohonest,
                                                                              test.data    = data.bin.mixed$test.data,
                                                                              true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                              noise.var    = data.bin.mixed$noise.var,
                                                                              corr.split   = data.bin.mixed$corr.split,
                                                                              where.split  = data.bin.mixed$where.split,
                                                                              dir.split    = data.bin.mixed$dir.split,
                                                                              split.cate   = data.bin.mixed$split.cate,
                                                                              CT           = T)
                            eval.ct.propsc.true.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
                            eval.ct.propsc.true.nohonest                <- unlist(eval.ct.propsc.true.nohonest)
                            names(eval.ct.propsc.true.nohonest)         <- paste0("homo.propsc.true.nohonest.", names(eval.ct.propsc.true.nohonest))
                            print("homo.true.nohonest")   
                            
                            #####################################################################################################################
                            ################################## 8. True propensity score model, honest Estimation ################################
                            #####################################################################################################################
                            t0 <- Sys.time()
                            tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                                                      method    = "GLM",
                                                      form.true = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))")
                            tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
                            
                            ct.propsc.true.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                                                       data             = train.data,
                                                                       weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                                                       treatment        = train.data$A,
                                                                       est_data         = est.data,
                                                                       est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                                                       est_treatment    = est.data$A,
                                                                       split.Rule       = "CT", 
                                                                       split.Honest     = T,
                                                                       HonestSampleSize = nrow(est.data),
                                                                       split.Bucket     = F,
                                                                       cv.option        = "CT",
                                                                       cv.Honest        = T)
                            cptable.propsc.true.honest <- ct.propsc.true.honest$cptable[,1][which.min(ct.propsc.true.honest$cptable[,4])]
                            final.tree.propsc.true.honest <- prune(ct.propsc.true.honest, cptable.propsc.true.honest)
                            t1 <- Sys.time()
                            
                            eval.ct.propsc.true.honest <- eval.measures.eff(final.tree   = final.tree.propsc.true.honest,
                                                                            test.data    = data.bin.mixed$test.data,
                                                                            true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                            noise.var    = data.bin.mixed$noise.var,
                                                                            corr.split   = data.bin.mixed$corr.split,
                                                                            where.split  = data.bin.mixed$where.split,
                                                                            dir.split    = data.bin.mixed$dir.split,
                                                                            split.cate   = data.bin.mixed$split.cate,
                                                                            CT           = T)
                            eval.ct.propsc.true.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
                            eval.ct.propsc.true.honest                <- unlist(eval.ct.propsc.true.honest)
                            names(eval.ct.propsc.true.honest)         <- paste0("homo.propsc.true.honest.", names(eval.ct.propsc.true.honest))
                            print("homo.true.honest")  
                            
                            #####################################################################################################################
                            ##################################### 9. Noisy propensity score model, no honest ####################################
                            #####################################################################################################################
                            t0 <- Sys.time()
                            tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                                                      method    = "GLM",
                                                      form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6")
                            tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
                            
                            ct.propsc.nois.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
                                                                  data         = data.used.full.bin.mixed,
                                                                  weights      = 1 / tmp.propsc$prop.sc,
                                                                  treatment    = data.used.full.bin.mixed$A,
                                                                  split.Rule   = "CT", 
                                                                  cv.option    = "CT", 
                                                                  split.Honest = T, 
                                                                  split.Bucket = F,
                                                                  cv.Honest    = T,
                                                                  xval         = 5, 
                                                                  cp           = 0, 
                                                                  minsize      = 20)
                            
                            cptable.propsc.nois.nohonest <- ct.propsc.nois.nohonest$cptable[,1][which.min(ct.propsc.nois.nohonest$cptable[,4])]
                            final.tree.propsc.nois.nohonest <- prune(ct.propsc.nois.nohonest, cptable.propsc.nois.nohonest)
                            t1 <- Sys.time()
                            
                            eval.ct.propsc.nois.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.nois.nohonest,
                                                                              test.data    = data.bin.mixed$test.data,
                                                                              true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                              noise.var    = data.bin.mixed$noise.var,
                                                                              corr.split   = data.bin.mixed$corr.split,
                                                                              where.split  = data.bin.mixed$where.split,
                                                                              dir.split    = data.bin.mixed$dir.split,
                                                                              split.cate   = data.bin.mixed$split.cate,
                                                                              CT           = T)
                            eval.ct.propsc.nois.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
                            eval.ct.propsc.nois.nohonest                <- unlist(eval.ct.propsc.nois.nohonest)
                            names(eval.ct.propsc.nois.nohonest)         <- paste0("homo.propsc.nois.nohonest.", names(eval.ct.propsc.nois.nohonest))
                            print("homo.nois.nohonest")  
                            
                            #####################################################################################################################
                            ################################## 10. Noisy propensity score model, honest estimation ##############################
                            #####################################################################################################################
                            t0 <- Sys.time()
                            tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                                                      method    = "GLM",
                                                      form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6")
                            tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
                            
                            ct.propsc.nois.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                                                       data             = train.data,
                                                                       weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                                                       treatment        = train.data$A,
                                                                       est_data         = est.data,
                                                                       est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                                                       est_treatment    = est.data$A,
                                                                       split.Rule       = "CT", 
                                                                       split.Honest     = T,
                                                                       HonestSampleSize = nrow(est.data),
                                                                       split.Bucket     = F,
                                                                       cv.option        = "CT",
                                                                       cv.Honest        = T)
                            cptable.propsc.nois.honest <- ct.propsc.nois.honest$cptable[,1][which.min(ct.propsc.nois.honest$cptable[,4])]
                            final.tree.propsc.nois.honest <- prune(ct.propsc.nois.honest, cptable.propsc.nois.honest)
                            t1 <- Sys.time()
                            
                            eval.ct.propsc.nois.honest <- eval.measures.eff(final.tree   = final.tree.propsc.nois.honest,
                                                                            test.data    = data.bin.mixed$test.data,
                                                                            true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                            noise.var    = data.bin.mixed$noise.var,
                                                                            corr.split   = data.bin.mixed$corr.split,
                                                                            where.split  = data.bin.mixed$where.split,
                                                                            dir.split    = data.bin.mixed$dir.split,
                                                                            split.cate   = data.bin.mixed$split.cate,
                                                                            CT           = T)
                            eval.ct.propsc.nois.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
                            eval.ct.propsc.nois.honest                <- unlist(eval.ct.propsc.nois.honest)
                            names(eval.ct.propsc.nois.honest)         <- paste0("homo.propsc.nois.honest.", names(eval.ct.propsc.nois.honest))
                            print("homo.nois.honest")   
                            
                            #####################################################################################################################
                            ################################## 11. Misspecified propensity score model, no honest ###############################
                            #####################################################################################################################
                            data.used.full.bin.mixed.mis <- data.used.full.bin.mixed %>%
                              select(-X2)
                            
                            t0 <- Sys.time()
                            tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed.mis[, !colnames(data.used.full.bin.mixed.mis) %in% c("Y")],
                                                      method    = "GLM",
                                                      form.true = NULL)
                            tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
                            
                            ct.propsc.mis.nohonest <- causalTree(Y ~ X1 + X3 + X4 + X5 + X6, 
                                                                 data         = data.used.full.bin.mixed.mis,
                                                                 weights      = 1 / tmp.propsc$prop.sc,
                                                                 treatment    = data.used.full.bin.mixed.mis$A,
                                                                 split.Rule   = "CT", 
                                                                 cv.option    = "CT", 
                                                                 split.Honest = T, 
                                                                 split.Bucket = F, 
                                                                 cv.Honest    = T,
                                                                 xval         = 5, 
                                                                 cp           = 0, 
                                                                 minsize      = 20)
                            
                            cptable.propsc.mis.nohonest <- ct.propsc.mis.nohonest$cptable[,1][which.min(ct.propsc.mis.nohonest$cptable[,4])]
                            final.tree.propsc.mis.nohonest <- prune(ct.propsc.mis.nohonest, cptable.propsc.mis.nohonest)
                            t1 <- Sys.time()
                            
                            eval.ct.propsc.mis.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.nohonest,
                                                                             test.data    = data.bin.mixed$test.data,
                                                                             true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                             noise.var    = data.bin.mixed$noise.var,
                                                                             corr.split   = data.bin.mixed$corr.split,
                                                                             where.split  = data.bin.mixed$where.split,
                                                                             dir.split    = data.bin.mixed$dir.split,
                                                                             split.cate   = data.bin.mixed$split.cate,
                                                                             CT           = T)
                            eval.ct.propsc.mis.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
                            eval.ct.propsc.mis.nohonest                <- unlist(eval.ct.propsc.mis.nohonest)
                            names(eval.ct.propsc.mis.nohonest)         <- paste0("homo.propsc.mis.nohonest.", names(eval.ct.propsc.mis.nohonest))
                            print("homo.mis.nohonest") 
                            
                            #####################################################################################################################
                            ############################ 12. Misspecified propensity score model, honest estimation #############################
                            #####################################################################################################################
                            train.data.mis <- train.data %>%
                              select(-X2)
                            est.data.mis   <- est.data %>%
                              select(-X2)
                            
                            t0 <- Sys.time()
                            tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed.mis[, !colnames(data.used.full.bin.mixed.mis) %in% c("Y")],
                                                      method    = "GLM",
                                                      form.true = NULL)
                            tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
                            
                            ct.propsc.mis.honest <- honest.causalTree(Y ~ X1 + X3 + X4 + X5 + X6,
                                                                      data             = train.data.mis,
                                                                      weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                                                      treatment        = train.data.mis$A,
                                                                      est_data         = est.data.mis,
                                                                      est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                                                      est_treatment    = est.data.mis$A,
                                                                      split.Rule       = "CT", 
                                                                      split.Honest     = T,
                                                                      HonestSampleSize = nrow(est.data.mis),
                                                                      split.Bucket     = F,
                                                                      cv.option        = "CT",
                                                                      cv.Honest        = T)
                            cptable.propsc.mis.honest <- ct.propsc.mis.honest$cptable[,1][which.min(ct.propsc.mis.honest$cptable[,4])]
                            final.tree.propsc.mis.honest <- prune(ct.propsc.mis.honest, cptable.propsc.mis.honest)
                            t1 <- Sys.time()
                            
                            eval.ct.propsc.mis.honest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.honest,
                                                                           test.data    = data.bin.mixed$test.data,
                                                                           true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                           noise.var    = data.bin.mixed$noise.var,
                                                                           corr.split   = data.bin.mixed$corr.split,
                                                                           where.split  = data.bin.mixed$where.split,
                                                                           dir.split    = data.bin.mixed$dir.split,
                                                                           split.cate   = data.bin.mixed$split.cate,
                                                                           CT           = T)
                            eval.ct.propsc.mis.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
                            eval.ct.propsc.mis.honest                <- unlist(eval.ct.propsc.mis.honest)
                            names(eval.ct.propsc.mis.honest)         <- paste0("homo.propsc.mis.honest.", names(eval.ct.propsc.mis.honest))
                            print("homo.mis.nohonest")   
                            
                            performance.homo.ct <- c(eval.ct.propsc.true.nohonest,
                                                     eval.ct.propsc.true.honest,
                                                     eval.ct.propsc.nois.nohonest,
                                                     eval.ct.propsc.nois.honest,
                                                     eval.ct.propsc.mis.nohonest,
                                                     eval.ct.propsc.mis.honest)
                            
                            # print must be put before output, otherwise the output will be print output
                            if (i%%60 == 0) {print(i)}
                            
                            c(performance.hetero.ct, performance.homo.ct)
                            
                          }

save(performance.ct, file = paste0("../Data/AppendixC7/BinCausalTreeOriginal.RData"))
