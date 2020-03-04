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

#####################################################################################################################
################################################# Heterogeneous #####################################################
#####################################################################################################################
data.bin.mixed            <- makeData.bin.eff.mixed(N             = 1000, 
                                                    n.test        = 1000, 
                                                    p.cont        = 3, 
                                                    p.cate        = 3, 
                                                    n.cate        = 4:6, 
                                                    coeff.prop.sc = 0.3,
                                                    seed          = a[job.number])
data.used.full.bin.mixed  <- data.bin.mixed$data.used

trtIdx  <- which(data.used.full.bin.mixed$A == 1)
ctrlIdx <- which(data.used.full.bin.mixed$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.bin.mixed[train.idx, ]
est.data   <- data.used.full.bin.mixed[-train.idx, ]

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
  tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                            method    = "GLM",
                            form.true = NULL)
  tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
  
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
                                      test.data    = data.bin.mixed$test.data,
                                      true.trt.eff = data.bin.mixed$true.trt.eff,
                                      noise.var    = data.bin.mixed$noise.var,
                                      corr.split   = data.bin.mixed$corr.split,
                                      where.split  = data.bin.mixed$where.split,
                                      dir.split    = data.bin.mixed$dir.split,
                                      split.cate   = data.bin.mixed$split.cate,
                                      CT           = T)
  eval.ct.propsc$t <- as.numeric(difftime(t1, t0, units = "secs"))
  eval.ct.propsc$corr.frst.splt <- as.character(ct.propsc$frame$var[1]) == data.bin.mixed$corr.split
  
  performance.hetero.ct <- c(performance.hetero.ct, list(eval.ct.propsc))
  names(performance.hetero.ct)[i] <- name.i
  print(name.i)
}

#####################################################################################################################
################################################## Homogeneous ######################################################
#####################################################################################################################
data.bin.mixed            <- makeData.bin.noeff.mixed(N             = 1000, 
                                                      n.test        = 1000, 
                                                      p.cont        = 3, 
                                                      p.cate        = 3, 
                                                      n.cate        = 4:6, 
                                                      coeff.prop.sc = 0.3,
                                                      seed          = a[job.number])
data.used.full.bin.mixed  <- data.bin.mixed$data.used

trtIdx  <- which(data.used.full.bin.mixed$A == 1)
ctrlIdx <- which(data.used.full.bin.mixed$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.bin.mixed[train.idx, ]
est.data   <- data.used.full.bin.mixed[-train.idx, ]

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
  tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                            method    = "GLM",
                            form.true = NULL)
  tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)
  
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
                                      test.data    = data.bin.mixed$test.data,
                                      true.trt.eff = data.bin.mixed$true.trt.eff,
                                      noise.var    = data.bin.mixed$noise.var,
                                      corr.split   = data.bin.mixed$corr.split,
                                      where.split  = data.bin.mixed$where.split,
                                      dir.split    = data.bin.mixed$dir.split,
                                      split.cate   = data.bin.mixed$split.cate,
                                      CT           = T)
  eval.ct.propsc$t <- as.numeric(difftime(t1, t0, units = "secs"))
  
  performance.homo.ct <- c(performance.homo.ct, list(eval.ct.propsc))
  names(performance.homo.ct)[i] <- name.i
  print(name.i)
}


#####################################################################################################################
################################################# Heterogeneous #####################################################
#####################################################################################################################
data.bin.mixed            <- makeData.bin.eff.mixed(N             = 1000, 
                                                    n.test        = 1000, 
                                                    p.cont        = 3, 
                                                    p.cate        = 3, 
                                                    n.cate        = 4:6, 
                                                    coeff.prop.sc = 0.3,
                                                    seed          = a[job.number])
data.used.full.bin.mixed  <- data.bin.mixed$data.used
data.used.bin.mixed       <- data.used.full.bin.mixed[1:800, ]
# val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
data.validation.bin.mixed <- data.used.full.bin.mixed[801:1000, ]  

# Pretend X2 is unmeasured for unmeasured cov
data.used.full.bin.mixed.mis <- data.used.full.bin.mixed %>%
  select(-X2)
test.data.mis <- data.bin.mixed$test.data %>%
  select(-X2)

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
print("hetero-true-nohonest")

#####################################################################################################################
################################## 2. True propensity score model, honest Estimation ################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.bin.mixed$A == 1)
ctrlIdx <- which(data.used.full.bin.mixed$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.bin.mixed[train.idx, ]
est.data   <- data.used.full.bin.mixed[-train.idx, ]

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
print("hetero-true-honest")

#####################################################################################################################
#################################### 3. Mis func propensity score model, no honest ##################################
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
print("hetero-nois-nohonest")

#####################################################################################################################
################################ 4. Mis func propensity score model, honest estimation ##############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.bin.mixed$A == 1)
ctrlIdx <- which(data.used.full.bin.mixed$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.bin.mixed[train.idx, ]
est.data   <- data.used.full.bin.mixed[-train.idx, ]

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
print("hetero-nois-honest")

#####################################################################################################################
################################## 5. Unmeasured cov propensity score model, no honest ##############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed.mis[, !colnames(data.used.full.bin.mixed.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.mis.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
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
                                                 test.data    = test.data.mis,
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
print("hetero-mis-nohonest")

#####################################################################################################################
############################ 6. Unmeasured cov propensity score model, honest estimation ############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed.mis[, !colnames(data.used.full.bin.mixed.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.bin.mixed.mis$A == 1)
ctrlIdx <- which(data.used.full.bin.mixed.mis$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.bin.mixed.mis[train.idx, ]
est.data   <- data.used.full.bin.mixed.mis[-train.idx, ]

ct.propsc.mis.honest <- honest.causalTree(Y ~ X1 + X3 + X4 + X5 + X6,
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
cptable.propsc.mis.honest <- ct.propsc.mis.honest$cptable[,1][which.min(ct.propsc.mis.honest$cptable[,4])]
final.tree.propsc.mis.honest <- prune(ct.propsc.mis.honest, cptable.propsc.mis.honest)
t1 <- Sys.time()

eval.ct.propsc.mis.honest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.honest,
                                               test.data    = test.data.mis,
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
print("hetero-mis-honest")

performance.hetero.ct <- list(true.nohonest = eval.ct.propsc.true.nohonest,
                              true.honest   = eval.ct.propsc.true.honest,
                              nois.nohonest = eval.ct.propsc.nois.nohonest,
                              nois.honest   = eval.ct.propsc.nois.honest,
                              mis.nohonest  = eval.ct.propsc.mis.nohonest,
                              mis.honest    = eval.ct.propsc.mis.honest)

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
                                      split.Rule   = "tstats", 
                                      split.Honest = F, 
                                      cv.option    = "fit", 
                                      cv.Honest    = F,
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
print("best-hetero-true-nohonest")

#####################################################################################################################
################################## 8. True propensity score model, honest Estimation ################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.bin.mixed$A == 1)
ctrlIdx <- which(data.used.full.bin.mixed$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.bin.mixed[train.idx, ]
est.data   <- data.used.full.bin.mixed[-train.idx, ]

ct.propsc.true.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                           data             = train.data,
                                           weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                           treatment        = train.data$A,
                                           est_data         = est.data,
                                           est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                           est_treatment    = est.data$A,
                                           split.Rule       = "tstats", 
                                           split.Honest     = F,
                                           HonestSampleSize = nrow(est.data),
                                           split.Bucket     = F,
                                           cv.option        = "fit",
                                           cv.Honest        = F)
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
print("best-hetero-true-honest")

#####################################################################################################################
#################################### 9. Mis func propensity score model, no honest ##################################
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
                                      split.Rule   = "tstats", 
                                      cv.option    = "fit", 
                                      split.Honest = F, 
                                      split.Bucket = F,
                                      cv.Honest    = F,
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
print("best-hetero-nois-nohonest")

#####################################################################################################################
################################ 10. Mis func propensity score model, honest estimation #############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.bin.mixed$A == 1)
ctrlIdx <- which(data.used.full.bin.mixed$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.bin.mixed[train.idx, ]
est.data   <- data.used.full.bin.mixed[-train.idx, ]

ct.propsc.nois.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                           data             = train.data,
                                           weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                           treatment        = train.data$A,
                                           est_data         = est.data,
                                           est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                           est_treatment    = est.data$A,
                                           split.Rule       = "tstats", 
                                           split.Honest     = F,
                                           HonestSampleSize = nrow(est.data),
                                           split.Bucket     = F,
                                           cv.option        = "fit",
                                           cv.Honest        = F)
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
print("best-hetero-nois-honest")

#####################################################################################################################
############################### 11. Unmeasured cov propensity score model, no honest ################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed.mis[, !colnames(data.used.full.bin.mixed.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.mis.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
                                     data         = data.used.full.bin.mixed.mis,
                                     weights      = 1 / tmp.propsc$prop.sc,
                                     treatment    = data.used.full.bin.mixed.mis$A,
                                     split.Rule   = "tstats", 
                                     cv.option    = "fit", 
                                     split.Honest = F, 
                                     split.Bucket = F, 
                                     cv.Honest    = F,
                                     xval         = 5, 
                                     cp           = 0, 
                                     minsize      = 20)

cptable.propsc.mis.nohonest <- ct.propsc.mis.nohonest$cptable[,1][which.min(ct.propsc.mis.nohonest$cptable[,4])]
final.tree.propsc.mis.nohonest <- prune(ct.propsc.mis.nohonest, cptable.propsc.mis.nohonest)
t1 <- Sys.time()

eval.ct.propsc.mis.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.nohonest,
                                                 test.data    = test.data.mis,
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
print("best-hetero-mis-nohonest")

#####################################################################################################################
############################ 12. Unmeasured cov propensity score model, honest estimation ###########################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed.mis[, !colnames(data.used.full.bin.mixed.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.bin.mixed.mis$A == 1)
ctrlIdx <- which(data.used.full.bin.mixed.mis$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.bin.mixed.mis[train.idx, ]
est.data   <- data.used.full.bin.mixed.mis[-train.idx, ]

ct.propsc.mis.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                          data             = train.data,
                                          weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                          treatment        = train.data$A,
                                          est_data         = est.data,
                                          est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                          est_treatment    = est.data$A,
                                          split.Rule       = "tstats", 
                                          split.Honest     = F,
                                          HonestSampleSize = nrow(est.data),
                                          split.Bucket     = F,
                                          cv.option        = "fit",
                                          cv.Honest        = F)
cptable.propsc.mis.honest <- ct.propsc.mis.honest$cptable[,1][which.min(ct.propsc.mis.honest$cptable[,4])]
final.tree.propsc.mis.honest <- prune(ct.propsc.mis.honest, cptable.propsc.mis.honest)
t1 <- Sys.time()

eval.ct.propsc.mis.honest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.honest,
                                               test.data    = test.data.mis,
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
print("best-hetero-mis-honest")

performance.hetero.ct.final <- list(true.nohonest = eval.ct.propsc.true.nohonest,
                                    true.honest   = eval.ct.propsc.true.honest,
                                    nois.nohonest = eval.ct.propsc.nois.nohonest,
                                    nois.honest   = eval.ct.propsc.nois.honest,
                                    mis.nohonest  = eval.ct.propsc.mis.nohonest,
                                    mis.honest    = eval.ct.propsc.mis.honest)

#####################################################################################################################
################################################## Homogeneous ######################################################
#####################################################################################################################
data.bin.mixed            <- makeData.bin.noeff.mixed(N             = 1000, 
                                                      n.test        = 1000, 
                                                      p.cont        = 3, 
                                                      p.cate        = 3, 
                                                      n.cate        = 4:6, 
                                                      coeff.prop.sc = 0.3,
                                                      seed          = a[job.number])
data.used.full.bin.mixed  <- data.bin.mixed$data.used
data.used.bin.mixed       <- data.used.full.bin.mixed[1:800, ]
# val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
data.validation.bin.mixed <- data.used.full.bin.mixed[801:1000, ]  

# Pretend X2 is unmeasured for unmeasured cov
data.used.full.bin.mixed.mis <- data.used.full.bin.mixed %>%
  select(-X2)
test.data.mis <- data.bin.mixed$test.data %>%
  select(-X2)

#####################################################################################################################
#################################### 13. True propensity score model, no honest #####################################
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
print("homo-true-nohonest")

#####################################################################################################################
################################## 14. True propensity score model, honest Estimation ###############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.bin.mixed$A == 1)
ctrlIdx <- which(data.used.full.bin.mixed$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.bin.mixed[train.idx, ]
est.data   <- data.used.full.bin.mixed[-train.idx, ]

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
print("homo-true-honest")

#####################################################################################################################
#################################### 15. Mis func propensity score model, no honest #################################
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
print("homo-nois-nohonest")

#####################################################################################################################
############################### 16. Mis func propensity score model, honest estimation ##############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.bin.mixed$A == 1)
ctrlIdx <- which(data.used.full.bin.mixed$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.bin.mixed[train.idx, ]
est.data   <- data.used.full.bin.mixed[-train.idx, ]

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
print("homo-nois-honest")

#####################################################################################################################
################################# 17. Unmeasured cov propensity score model, no honest ##############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed.mis[, !colnames(data.used.full.bin.mixed.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.mis.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
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
                                                 test.data    = test.data.mis,
                                                 true.trt.eff = data.bin.mixed$true.trt.eff,
                                                 noise.var    = data.bin.mixed$noise.var,
                                                 corr.split   = data.bin.mixed$corr.split,
                                                 where.split  = data.bin.mixed$where.split,
                                                 dir.split    = data.bin.mixed$dir.split,
                                                 split.cate   = data.bin.mixed$split.cate,
                                                 CT           = T)
eval.ct.propsc.mis.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("homo-mis-nohonest")

#####################################################################################################################
############################ 18. Unmeasured cov propensity score model, honest estimation ###########################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed.mis[, !colnames(data.used.full.bin.mixed.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.bin.mixed.mis$A == 1)
ctrlIdx <- which(data.used.full.bin.mixed.mis$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.bin.mixed.mis[train.idx, ]
est.data   <- data.used.full.bin.mixed.mis[-train.idx, ]

ct.propsc.mis.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
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
cptable.propsc.mis.honest <- ct.propsc.mis.honest$cptable[,1][which.min(ct.propsc.mis.honest$cptable[,4])]
final.tree.propsc.mis.honest <- prune(ct.propsc.mis.honest, cptable.propsc.mis.honest)
t1 <- Sys.time()

eval.ct.propsc.mis.honest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.honest,
                                               test.data    = test.data.mis,
                                               true.trt.eff = data.bin.mixed$true.trt.eff,
                                               noise.var    = data.bin.mixed$noise.var,
                                               corr.split   = data.bin.mixed$corr.split,
                                               where.split  = data.bin.mixed$where.split,
                                               dir.split    = data.bin.mixed$dir.split,
                                               split.cate   = data.bin.mixed$split.cate,
                                               CT           = T)
eval.ct.propsc.mis.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("homo-mis-honest")

performance.homo.ct <- list(true.nohonest = eval.ct.propsc.true.nohonest,
                            true.honest   = eval.ct.propsc.true.honest,
                            nois.nohonest = eval.ct.propsc.nois.nohonest,
                            nois.honest   = eval.ct.propsc.nois.honest,
                            mis.nohonest  = eval.ct.propsc.mis.nohonest,
                            mis.honest    = eval.ct.propsc.mis.honest)

#####################################################################################################################
#################################### 19. True propensity score model, no honest #####################################
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
                                      split.Rule   = "tstats", 
                                      split.Honest = F, 
                                      cv.option    = "fit", 
                                      cv.Honest    = F,
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
print("best-homo-true-nohonest")

#####################################################################################################################
################################## 20. True propensity score model, honest Estimation ###############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.bin.mixed$A == 1)
ctrlIdx <- which(data.used.full.bin.mixed$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.bin.mixed[train.idx, ]
est.data   <- data.used.full.bin.mixed[-train.idx, ]

ct.propsc.true.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                           data             = train.data,
                                           weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                           treatment        = train.data$A,
                                           est_data         = est.data,
                                           est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                           est_treatment    = est.data$A,
                                           split.Rule       = "tstats", 
                                           split.Honest     = F,
                                           HonestSampleSize = nrow(est.data),
                                           split.Bucket     = F,
                                           cv.option        = "fit",
                                           cv.Honest        = F)
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
print("best-homo-true-honest")

#####################################################################################################################
#################################### 21. Mis func propensity score model, no honest #################################
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
                                      split.Rule   = "tstats", 
                                      cv.option    = "fit", 
                                      split.Honest = F, 
                                      split.Bucket = F,
                                      cv.Honest    = F,
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
print("best-homo-nois-nohonest")

#####################################################################################################################
############################### 22. Mis func propensity score model, honest estimation ##############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed[, !colnames(data.used.full.bin.mixed) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.bin.mixed$A == 1)
ctrlIdx <- which(data.used.full.bin.mixed$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.bin.mixed[train.idx, ]
est.data   <- data.used.full.bin.mixed[-train.idx, ]

ct.propsc.nois.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                           data             = train.data,
                                           weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                           treatment        = train.data$A,
                                           est_data         = est.data,
                                           est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                           est_treatment    = est.data$A,
                                           split.Rule       = "tstats", 
                                           split.Honest     = F,
                                           HonestSampleSize = nrow(est.data),
                                           split.Bucket     = F,
                                           cv.option        = "fit",
                                           cv.Honest        = F)
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
print("best-homo-nois-honest")

#####################################################################################################################
################################# 23. Unmeasured cov propensity score model, no honest ##############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed.mis[, !colnames(data.used.full.bin.mixed.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.mis.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
                                     data         = data.used.full.bin.mixed.mis,
                                     weights      = 1 / tmp.propsc$prop.sc,
                                     treatment    = data.used.full.bin.mixed.mis$A,
                                     split.Rule   = "tstats", 
                                     cv.option    = "fit", 
                                     split.Honest = F, 
                                     split.Bucket = F, 
                                     cv.Honest    = F,
                                     xval         = 5, 
                                     cp           = 0, 
                                     minsize      = 20)

cptable.propsc.mis.nohonest <- ct.propsc.mis.nohonest$cptable[,1][which.min(ct.propsc.mis.nohonest$cptable[,4])]
final.tree.propsc.mis.nohonest <- prune(ct.propsc.mis.nohonest, cptable.propsc.mis.nohonest)
t1 <- Sys.time()

eval.ct.propsc.mis.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.nohonest,
                                                 test.data    = test.data.mis,
                                                 true.trt.eff = data.bin.mixed$true.trt.eff,
                                                 noise.var    = data.bin.mixed$noise.var,
                                                 corr.split   = data.bin.mixed$corr.split,
                                                 where.split  = data.bin.mixed$where.split,
                                                 dir.split    = data.bin.mixed$dir.split,
                                                 split.cate   = data.bin.mixed$split.cate,
                                                 CT           = T)
eval.ct.propsc.mis.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("best-homo-mis-nohonest")

#####################################################################################################################
############################ 24. Unmeasured cov propensity score model, honest estimation ###########################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.bin.mixed.mis[, !colnames(data.used.full.bin.mixed.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.bin.mixed.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.bin.mixed.mis$A == 1)
ctrlIdx <- which(data.used.full.bin.mixed.mis$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.bin.mixed.mis[train.idx, ]
est.data   <- data.used.full.bin.mixed.mis[-train.idx, ]

ct.propsc.mis.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                          data             = train.data,
                                          weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                          treatment        = train.data$A,
                                          est_data         = est.data,
                                          est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                          est_treatment    = est.data$A,
                                          split.Rule       = "tstats", 
                                          split.Honest     = F,
                                          HonestSampleSize = nrow(est.data),
                                          split.Bucket     = F,
                                          cv.option        = "fit",
                                          cv.Honest        = F)
cptable.propsc.mis.honest <- ct.propsc.mis.honest$cptable[,1][which.min(ct.propsc.mis.honest$cptable[,4])]
final.tree.propsc.mis.honest <- prune(ct.propsc.mis.honest, cptable.propsc.mis.honest)
t1 <- Sys.time()

eval.ct.propsc.mis.honest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.honest,
                                               test.data    = test.data.mis,
                                               true.trt.eff = data.bin.mixed$true.trt.eff,
                                               noise.var    = data.bin.mixed$noise.var,
                                               corr.split   = data.bin.mixed$corr.split,
                                               where.split  = data.bin.mixed$where.split,
                                               dir.split    = data.bin.mixed$dir.split,
                                               split.cate   = data.bin.mixed$split.cate,
                                               CT           = T)
eval.ct.propsc.mis.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("best-homo-mis-honest")

performance.homo.ct.final <- list(true.nohonest = eval.ct.propsc.true.nohonest,
                                  true.honest   = eval.ct.propsc.true.honest,
                                  nois.nohonest = eval.ct.propsc.nois.nohonest,
                                  nois.honest   = eval.ct.propsc.nois.honest,
                                  mis.nohonest  = eval.ct.propsc.mis.nohonest,
                                  mis.honest    = eval.ct.propsc.mis.honest)


