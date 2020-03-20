job.number <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
load("../seed1000.rda")
set.seed(a[job.number])

setwd("..")
folder <- paste(getwd(), "/Functions/", sep="")
functions <- list.files(folder)
functions <- paste(folder, functions, sep = "")
for (i in functions){
  source(i)
}

setwd("Appendix/")

# Causal Tree
# library(devtools) 
# install_github("susanathey/causalTree")
library(causalTree)

coeff.prop.sc <- 0.6
N.training    <- 10^3
N.testing     <- 10^3

#####################################################################################################################
################################################# Heterogeneous #####################################################
#####################################################################################################################
data.cont.cont            <- makeData.cont.eff.cont(N             = N.training, 
                                                    n.test        = N.testing, 
                                                    coeff.prop.sc = coeff.prop.sc)
data.used.full.cont.cont  <- data.cont.cont$data.used
data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
# val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ]  

trtIdx  <- which(data.used.full.cont.cont$A == 1)
ctrlIdx <- which(data.used.full.cont.cont$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.cont.cont[train.idx, ]
est.data   <- data.used.full.cont.cont[-train.idx, ]

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

performance.hetero.ct <- list()
for (i in 1:nrow(comb.scnr)) {
  
  # Use noise model to evaluate different settings
  t0 <- Sys.time()
  tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                            method    = "GLM",
                            form.true = NULL)
  tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
  
  if (is.na(comb.scnr$split_honest[i])) {
    
     if (is.na(comb.scnr$cv_honest[i])) {
       ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                      data             = train.data,
                                      weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                      treatment        = train.data$A,
                                      est_data         = est.data,
                                      est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                      est_treatment    = est.data$A,
                                      split.Rule       = comb.scnr$split_rule[i], 
                                      HonestSampleSize = nrow(est.data),
                                      split.Bucket     = F,
                                      cv.option        = comb.scnr$cv_option[i])
       
       name.i <- paste(tolower(comb.scnr$split_rule[i]), tolower(comb.scnr$cv_option[i]), sep = ".")
     } else {
       ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                      data             = train.data,
                                      weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                      treatment        = train.data$A,
                                      est_data         = est.data,
                                      est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                      est_treatment    = est.data$A,
                                      split.Rule       = comb.scnr$split_rule[i], 
                                      HonestSampleSize = nrow(est.data),
                                      split.Bucket     = F,
                                      cv.option        = comb.scnr$cv_option[i],
                                      cv.Honest        = comb.scnr$cv_honest[i])
       name.i <- paste(tolower(comb.scnr$split_rule[i]), paste0(tolower(comb.scnr$cv_option[i]), substring(comb.scnr$cv_honest[i], 1, 1)), sep = ".")
     }
    
  } else {
    
    if (is.na(comb.scnr$cv_honest[i])) {
      ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                     data             = train.data,
                                     weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                     treatment        = train.data$A,
                                     est_data         = est.data,
                                     est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                     est_treatment    = est.data$A,
                                     split.Rule       = comb.scnr$split_rule[i], 
                                     split.Honest     = comb.scnr$split_honest[i],
                                     HonestSampleSize = nrow(est.data),
                                     split.Bucket     = F,
                                     cv.option        = comb.scnr$cv_option[i])
      name.i <- paste(paste0(tolower(comb.scnr$split_rule[i]), 
                             substring(comb.scnr$split_honest[i], 1, 1)), 
                      tolower(comb.scnr$cv_option[i]), sep = ".")
      
    } else {
      ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                     data             = train.data,
                                     weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                     treatment        = train.data$A,
                                     est_data         = est.data,
                                     est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                     est_treatment    = est.data$A,
                                     split.Rule       = comb.scnr$split_rule[i], 
                                     split.Honest     = comb.scnr$split_honest[i],
                                     HonestSampleSize = nrow(est.data),
                                     split.Bucket     = F,
                                     cv.option        = comb.scnr$cv_option[i],
                                     cv.Honest        = comb.scnr$cv_honest[i])
      name.i <- paste(paste0(tolower(comb.scnr$split_rule[i]), 
                             substring(comb.scnr$split_honest[i], 1, 1)), 
                      paste0(tolower(comb.scnr$cv_option[i]), 
                             substring(comb.scnr$cv_honest[i], 1, 1)), sep = ".")
      
    }
    
  }
  
  cptable.propsc <- ct.propsc$cptable[,1][which.min(ct.propsc$cptable[,4])]
  final.tree.propsc <- prune(ct.propsc, cptable.propsc)
  t1 <- Sys.time()
  
  eval.ct.propsc <- eval.measures.eff(final.tree   = final.tree.propsc,
                                      test.data    = data.cont.cont$test.data,
                                      true.trt.eff = data.cont.cont$true.trt.eff,
                                      noise.var    = data.cont.cont$noise.var,
                                      corr.split   = data.cont.cont$corr.split,
                                      where.split  = data.cont.cont$where.split,
                                      dir.split    = data.cont.cont$dir.split)
  eval.ct.propsc$t <- as.numeric(difftime(t1, t0, units = "secs"))
  eval.ct.propsc$corr.frst.splt <- as.character(ct.propsc$frame$var[1]) == data.cont.cont$corr.split
  
  performance.hetero.ct <- c(performance.hetero.ct, list(eval.ct.propsc))
  names(performance.hetero.ct)[i] <- name.i
  print(name.i)
}

#####################################################################################################################
################################################### Homogeneous #####################################################
#####################################################################################################################
data.cont.cont            <- makeData.cont.noeff.cont(N             = N.training, 
                                                      n.test        = N.testing, 
                                                      coeff.prop.sc = coeff.prop.sc)
data.used.full.cont.cont  <- data.cont.cont$data.used
data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
# val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ]  

trtIdx  <- which(data.used.full.cont.cont$A == 1)
ctrlIdx <- which(data.used.full.cont.cont$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.cont.cont[train.idx, ]
est.data   <- data.used.full.cont.cont[-train.idx, ]

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

performance.homo.ct <- list()
for (i in 1:nrow(comb.scnr)) {
  
  t0 <- Sys.time()
  tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                            method    = "GLM",
                            form.true = NULL)
  tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
  
  if (is.na(comb.scnr$split_honest[i])) {
    
    if (is.na(comb.scnr$cv_honest[i])) {
      ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                     data             = train.data,
                                     weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                     treatment        = train.data$A,
                                     est_data         = est.data,
                                     est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                     est_treatment    = est.data$A,
                                     split.Rule       = comb.scnr$split_rule[i], 
                                     HonestSampleSize = nrow(est.data),
                                     split.Bucket     = F,
                                     cv.option        = comb.scnr$cv_option[i])
      name.i <- paste(tolower(comb.scnr$split_rule[i]), tolower(comb.scnr$cv_option[i]), sep = ".")
    } else {
      ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                     data             = train.data,
                                     weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                     treatment        = train.data$A,
                                     est_data         = est.data,
                                     est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                     est_treatment    = est.data$A,
                                     split.Rule       = comb.scnr$split_rule[i], 
                                     HonestSampleSize = nrow(est.data),
                                     split.Bucket     = F,
                                     cv.option        = comb.scnr$cv_option[i],
                                     cv.Honest        = comb.scnr$cv_honest[i])
      name.i <- paste(tolower(comb.scnr$split_rule[i]), paste0(tolower(comb.scnr$cv_option[i]), substring(comb.scnr$cv_honest[i], 1, 1)), sep = ".")
    }
    
  } else {
    
    if (is.na(comb.scnr$cv_honest[i])) {
      ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                     data             = train.data,
                                     weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                     treatment        = train.data$A,
                                     est_data         = est.data,
                                     est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                     est_treatment    = est.data$A,
                                     split.Rule       = comb.scnr$split_rule[i], 
                                     split.Honest     = comb.scnr$split_honest[i],
                                     HonestSampleSize = nrow(est.data),
                                     split.Bucket     = F,
                                     cv.option        = comb.scnr$cv_option[i])
      name.i <- paste(paste0(tolower(comb.scnr$split_rule[i]), 
                             substring(comb.scnr$split_honest[i], 1, 1)), 
                      tolower(comb.scnr$cv_option[i]), sep = ".")
      
    } else {
      ct.propsc <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                     data             = train.data,
                                     weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                     treatment        = train.data$A,
                                     est_data         = est.data,
                                     est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                     est_treatment    = est.data$A,
                                     split.Rule       = comb.scnr$split_rule[i], 
                                     split.Honest     = comb.scnr$split_honest[i],
                                     HonestSampleSize = nrow(est.data),
                                     split.Bucket     = F,
                                     cv.option        = comb.scnr$cv_option[i],
                                     cv.Honest        = comb.scnr$cv_honest[i])
      name.i <- paste(paste0(tolower(comb.scnr$split_rule[i]), 
                             substring(comb.scnr$split_honest[i], 1, 1)), 
                      paste0(tolower(comb.scnr$cv_option[i]), 
                             substring(comb.scnr$cv_honest[i], 1, 1)), sep = ".")
      
    }
    
  }
  
  cptable.propsc <- ct.propsc$cptable[,1][which.min(ct.propsc$cptable[,4])]
  final.tree.propsc <- prune(ct.propsc, cptable.propsc)
  t1 <- Sys.time()
  
  eval.ct.propsc <- eval.measures.eff(final.tree   = final.tree.propsc,
                                      test.data    = data.cont.cont$test.data,
                                      true.trt.eff = data.cont.cont$true.trt.eff,
                                      noise.var    = data.cont.cont$noise.var,
                                      corr.split   = data.cont.cont$corr.split,
                                      where.split  = data.cont.cont$where.split,
                                      dir.split    = data.cont.cont$dir.split)
  eval.ct.propsc$t <- as.numeric(difftime(t1, t0, units = "secs"))
  
  performance.homo.ct <- c(performance.homo.ct, list(eval.ct.propsc))
  names(performance.homo.ct)[i] <- name.i
  print(name.i)
}

file.name = paste("../Data/AppendixC1/", toString(job.number), ".RData", sep = "")
save(performance.homo.ct, performance.hetero.ct,
     file = file.name)
