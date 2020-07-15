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

performance.ct <- foreach(i = start:end,
                          .combine  = "rbind") %dopar% {
                            
                            set.seed(a[i])
                            
                            #####################################################################################################################
                            ################################################# Heterogeneous #####################################################
                            #####################################################################################################################
                            performance.hetero.ct     <- NULL
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
                            
                            t0 <- Sys.time()
                            # Estimate propensity scores
                            tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                                                      method    = "GLM",
                                                      form.true = NULL)
                            tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
                            
                            for (comb.scnr.i in 1:nrow(comb.scnr)) {
                              
                              if (is.na(comb.scnr$split_honest[comb.scnr.i])) {
                                
                                if (is.na(comb.scnr$cv_honest[comb.scnr.i])) {
                                  ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                                                 data             = train.data,
                                                                 weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                                                 treatment        = train.data$A,
                                                                 est_data         = est.data,
                                                                 est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                                                 est_treatment    = est.data$A,
                                                                 split.Rule       = comb.scnr$split_rule[comb.scnr.i], 
                                                                 HonestSampleSize = nrow(est.data),
                                                                 split.Bucket     = F,
                                                                 cv.option        = comb.scnr$cv_option[comb.scnr.i])
                                  
                                  name.i <- paste(tolower(comb.scnr$split_rule[comb.scnr.i]), tolower(comb.scnr$cv_option[comb.scnr.i]), sep = ".")
                                } else {
                                  ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                                                 data             = train.data,
                                                                 weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                                                 treatment        = train.data$A,
                                                                 est_data         = est.data,
                                                                 est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                                                 est_treatment    = est.data$A,
                                                                 split.Rule       = comb.scnr$split_rule[comb.scnr.i], 
                                                                 HonestSampleSize = nrow(est.data),
                                                                 split.Bucket     = F,
                                                                 cv.option        = comb.scnr$cv_option[comb.scnr.i],
                                                                 cv.Honest        = comb.scnr$cv_honest[comb.scnr.i])
                                  name.i <- paste(tolower(comb.scnr$split_rule[comb.scnr.i]), paste0(tolower(comb.scnr$cv_option[comb.scnr.i]), substring(comb.scnr$cv_honest[comb.scnr.i], 1, 1)), sep = ".")
                                }
                              } else {
                                
                                if (is.na(comb.scnr$cv_honest[comb.scnr.i])) {
                                  ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                                                 data             = train.data,
                                                                 weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                                                 treatment        = train.data$A,
                                                                 est_data         = est.data,
                                                                 est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                                                 est_treatment    = est.data$A,
                                                                 split.Rule       = comb.scnr$split_rule[comb.scnr.i], 
                                                                 split.Honest     = comb.scnr$split_honest[comb.scnr.i],
                                                                 HonestSampleSize = nrow(est.data),
                                                                 split.Bucket     = F,
                                                                 cv.option        = comb.scnr$cv_option[comb.scnr.i])
                                  name.i <- paste(paste0(tolower(comb.scnr$split_rule[comb.scnr.i]), 
                                                         substring(comb.scnr$split_honest[comb.scnr.i], 1, 1)), 
                                                  tolower(comb.scnr$cv_option[comb.scnr.i]), sep = ".")
                                  
                                } else {
                                  ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                                                 data             = train.data,
                                                                 weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                                                 treatment        = train.data$A,
                                                                 est_data         = est.data,
                                                                 est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                                                 est_treatment    = est.data$A,
                                                                 split.Rule       = comb.scnr$split_rule[comb.scnr.i], 
                                                                 split.Honest     = comb.scnr$split_honest[comb.scnr.i],
                                                                 HonestSampleSize = nrow(est.data),
                                                                 split.Bucket     = F,
                                                                 cv.option        = comb.scnr$cv_option[comb.scnr.i],
                                                                 cv.Honest        = comb.scnr$cv_honest[comb.scnr.i])
                                  name.i <- paste(paste0(tolower(comb.scnr$split_rule[comb.scnr.i]), 
                                                         substring(comb.scnr$split_honest[comb.scnr.i], 1, 1)), 
                                                  paste0(tolower(comb.scnr$cv_option[comb.scnr.i]), 
                                                         substring(comb.scnr$cv_honest[comb.scnr.i], 1, 1)), sep = ".")
                                  
                                }
                                
                              }
                              
                              cptable.propsc <- ct.propsc$cptable[,1][which.min(ct.propsc$cptable[,4])]
                              final.tree.propsc <- prune(ct.propsc, cptable.propsc)
                              t1 <- Sys.time()
                              
                              eval.ct.propsc <- eval.measures.eff(final.tree   = final.tree.propsc,
                                                                  test.data    = data.bin.mixed$test.data,
                                                                  true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                  noise.var    = data.bin.mixed$noise.var,
                                                                  corr.split   = data.bin.mixed$corr.split,
                                                                  where.split  = data.bin.mixed$where.split,
                                                                  dir.split    = data.bin.mixed$dir.split,
                                                                  split.cate   = data.bin.mixed$split.cate,
                                                                  CT           = T)
                              eval.ct.propsc$t              <- as.numeric(difftime(t1, t0, units = "secs"))
                              eval.ct.propsc$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = ct.propsc,
                                                                                        corr.split = data.bin.mixed$corr.split,
                                                                                        split.cate = data.bin.mixed$split.cate)
                              
                              eval.ct.propsc                <- unlist(eval.ct.propsc)
                              names(eval.ct.propsc)         <- paste0("hetero.", name.i, ".", names(eval.ct.propsc))
                              
                              performance.hetero.ct <- c(performance.hetero.ct, eval.ct.propsc)
                              print(name.i)
                              
                            } # for (comb.scnr.i in 1:nrow(comb.scnr)) {
                            
                            #####################################################################################################################
                            ################################################### Homogeneous #####################################################
                            #####################################################################################################################
                            performance.homo.ct     <- NULL
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
                            
                            t0 <- Sys.time()
                            # Estimate propensity scores
                            tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                                                      method    = "GLM",
                                                      form.true = NULL)
                            tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
                            
                            for (comb.scnr.i in 1:nrow(comb.scnr)) {
                              
                              if (is.na(comb.scnr$split_honest[comb.scnr.i])) {
                                
                                if (is.na(comb.scnr$cv_honest[comb.scnr.i])) {
                                  ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                                                 data             = train.data,
                                                                 weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                                                 treatment        = train.data$A,
                                                                 est_data         = est.data,
                                                                 est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                                                 est_treatment    = est.data$A,
                                                                 split.Rule       = comb.scnr$split_rule[comb.scnr.i], 
                                                                 HonestSampleSize = nrow(est.data),
                                                                 split.Bucket     = F,
                                                                 cv.option        = comb.scnr$cv_option[comb.scnr.i])
                                  
                                  name.i <- paste(tolower(comb.scnr$split_rule[comb.scnr.i]), tolower(comb.scnr$cv_option[comb.scnr.i]), sep = ".")
                                } else {
                                  ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                                                 data             = train.data,
                                                                 weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                                                 treatment        = train.data$A,
                                                                 est_data         = est.data,
                                                                 est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                                                 est_treatment    = est.data$A,
                                                                 split.Rule       = comb.scnr$split_rule[comb.scnr.i], 
                                                                 HonestSampleSize = nrow(est.data),
                                                                 split.Bucket     = F,
                                                                 cv.option        = comb.scnr$cv_option[comb.scnr.i],
                                                                 cv.Honest        = comb.scnr$cv_honest[comb.scnr.i])
                                  name.i <- paste(tolower(comb.scnr$split_rule[comb.scnr.i]), paste0(tolower(comb.scnr$cv_option[comb.scnr.i]), substring(comb.scnr$cv_honest[comb.scnr.i], 1, 1)), sep = ".")
                                }
                              } else {
                                
                                if (is.na(comb.scnr$cv_honest[comb.scnr.i])) {
                                  ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                                                 data             = train.data,
                                                                 weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                                                 treatment        = train.data$A,
                                                                 est_data         = est.data,
                                                                 est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                                                 est_treatment    = est.data$A,
                                                                 split.Rule       = comb.scnr$split_rule[comb.scnr.i], 
                                                                 split.Honest     = comb.scnr$split_honest[comb.scnr.i],
                                                                 HonestSampleSize = nrow(est.data),
                                                                 split.Bucket     = F,
                                                                 cv.option        = comb.scnr$cv_option[comb.scnr.i])
                                  name.i <- paste(paste0(tolower(comb.scnr$split_rule[comb.scnr.i]), 
                                                         substring(comb.scnr$split_honest[comb.scnr.i], 1, 1)), 
                                                  tolower(comb.scnr$cv_option[comb.scnr.i]), sep = ".")
                                  
                                } else {
                                  ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                                                 data             = train.data,
                                                                 weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                                                 treatment        = train.data$A,
                                                                 est_data         = est.data,
                                                                 est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                                                 est_treatment    = est.data$A,
                                                                 split.Rule       = comb.scnr$split_rule[comb.scnr.i], 
                                                                 split.Honest     = comb.scnr$split_honest[comb.scnr.i],
                                                                 HonestSampleSize = nrow(est.data),
                                                                 split.Bucket     = F,
                                                                 cv.option        = comb.scnr$cv_option[comb.scnr.i],
                                                                 cv.Honest        = comb.scnr$cv_honest[comb.scnr.i])
                                  name.i <- paste(paste0(tolower(comb.scnr$split_rule[comb.scnr.i]), 
                                                         substring(comb.scnr$split_honest[comb.scnr.i], 1, 1)), 
                                                  paste0(tolower(comb.scnr$cv_option[comb.scnr.i]), 
                                                         substring(comb.scnr$cv_honest[comb.scnr.i], 1, 1)), sep = ".")
                                  
                                }
                                
                              }
                              
                              cptable.propsc <- ct.propsc$cptable[,1][which.min(ct.propsc$cptable[,4])]
                              final.tree.propsc <- prune(ct.propsc, cptable.propsc)
                              t1 <- Sys.time()
                              
                              eval.ct.propsc <- eval.measures.eff(final.tree   = final.tree.propsc,
                                                                  test.data    = data.bin.mixed$test.data,
                                                                  true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                  noise.var    = data.bin.mixed$noise.var,
                                                                  corr.split   = data.bin.mixed$corr.split,
                                                                  where.split  = data.bin.mixed$where.split,
                                                                  dir.split    = data.bin.mixed$dir.split,
                                                                  split.cate   = data.bin.mixed$split.cate,
                                                                  CT           = T)
                              eval.ct.propsc$t              <- as.numeric(difftime(t1, t0, units = "secs"))
                              eval.ct.propsc                <- unlist(eval.ct.propsc)
                              names(eval.ct.propsc)         <- paste0("homo.", name.i, ".", names(eval.ct.propsc))
                              
                              performance.homo.ct <- c(performance.homo.ct, eval.ct.propsc)
                              print(name.i)
                              
                            } # for (comb.scnr.i in 1:nrow(comb.scnr)) {
                            
                            # print must be put before output, otherwise the output will be print output
                            if (i %% 60 == 0) {print(i)}
                            
                            c(performance.hetero.ct, performance.homo.ct)
                            
                          }

save(performance.ct, file = paste0("../Data/AppendixC7/BinCausalTreeSettings.RData"))
