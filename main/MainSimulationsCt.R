#!/usr/bin/env Rscript
setwd("..")
folder <- paste(getwd(), "/Functions/", sep="")
functions <- list.files(folder)
functions <- paste(folder, functions, sep = "")
for (i in functions){
  source(i)
}

setwd("main/")

# Read in the arguments from command line
option_list = list(
  make_option(c("-a", "--array"), type="integer"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
iters <- opt$array

job.number <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) + iters
load("../seed1000.rda")
set.seed(a[job.number])
print(job.number)

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

# Pretend X2 is unmeasured for unmeasured cov
data.used.full.cont.cont.mis <- data.used.full.cont.cont %>%
  select(-X2)
test.data.mis <- data.cont.cont$test.data %>%
  select(-X2)

#####################################################################################################################
###################################### True propensity score model, no honest #######################################
#####################################################################################################################
# In help document: Unit-specific propensity scores are not supported
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X2 + X3")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.true.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
                                      data         = data.used.full.cont.cont,
                                      weights      = 1 / tmp.propsc$prop.sc,
                                      treatment    = data.used.full.cont.cont$A,
                                      split.Rule   = "CT", 
                                      cv.option    = "CT", 
                                      split.Honest = T, 
                                      cv.Honest    = T, 
                                      split.Bucket = F, 
                                      xval         = 5, 
                                      cp           = 0, 
                                      minsize      = 20)
cptable.propsc.true.nohonest <- ct.propsc.true.nohonest$cptable[,1][which.min(ct.propsc.true.nohonest$cptable[,4])]
final.tree.propsc.true.nohonest <- prune(ct.propsc.true.nohonest, cptable.propsc.true.nohonest)
t1 <- Sys.time()

eval.ct.propsc.true.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.true.nohonest,
                                                  test.data    = data.cont.cont$test.data,
                                                  true.trt.eff = data.cont.cont$true.trt.eff,
                                                  noise.var    = data.cont.cont$noise.var,
                                                  corr.split   = data.cont.cont$corr.split,
                                                  where.split  = data.cont.cont$where.split,
                                                  dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.true.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.ct.propsc.true.nohonest$corr.frst.splt <- as.character(ct.propsc.true.nohonest$frame$var[1]) == data.cont.cont$corr.split
# mean((predict(final.tree.propsc.true.nohonest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("1")

#####################################################################################################################
#################################### Mis Func propensity score model, no honest #####################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.nois.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
                                      data       = data.used.full.cont.cont,
                                      weights    = 1 / tmp.propsc$prop.sc,
                                      treatment  = data.used.full.cont.cont$A,
                                      split.Rule = "CT", 
                                      cv.option  = "CT", 
                                      split.Honest = T, 
                                      cv.Honest = T, 
                                      split.Bucket = F, 
                                      xval = 5, 
                                      cp = 0, 
                                      minsize = 20)

cptable.propsc.nois.nohonest <- ct.propsc.nois.nohonest$cptable[,1][which.min(ct.propsc.nois.nohonest$cptable[,4])]
final.tree.propsc.nois.nohonest <- prune(ct.propsc.nois.nohonest, cptable.propsc.nois.nohonest)
t1 <- Sys.time()

eval.ct.propsc.nois.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.nois.nohonest,
                                                  test.data    = data.cont.cont$test.data,
                                                  true.trt.eff = data.cont.cont$true.trt.eff,
                                                  noise.var    = data.cont.cont$noise.var,
                                                  corr.split   = data.cont.cont$corr.split,
                                                  where.split  = data.cont.cont$where.split,
                                                  dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.nois.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.ct.propsc.nois.nohonest$corr.frst.splt <- as.character(ct.propsc.nois.nohonest$frame$var[1]) == data.cont.cont$corr.split
# mean((predict(final.tree.propsc.nois.nohonest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("2")

#####################################################################################################################
################################## Unmeasured cov propensity score model, no honest #################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont.mis[, !colnames(data.used.full.cont.cont.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.mis.nohonest <- causalTree(Y ~ X1 + X3 + X4 + X5 + X6, 
                                     data       = data.used.full.cont.cont.mis,
                                     weights    = 1 / tmp.propsc$prop.sc,
                                     treatment  = data.used.full.cont.cont.mis$A,
                                     split.Rule = "CT", 
                                     cv.option  = "CT", 
                                     split.Honest = T, 
                                     cv.Honest = T, 
                                     split.Bucket = F, 
                                     xval = 5, 
                                     cp = 0, 
                                     minsize = 20)

cptable.propsc.mis.nohonest <- ct.propsc.mis.nohonest$cptable[,1][which.min(ct.propsc.mis.nohonest$cptable[,4])]
final.tree.propsc.mis.nohonest <- prune(ct.propsc.mis.nohonest, cptable.propsc.mis.nohonest)
t1 <- Sys.time()

eval.ct.propsc.mis.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.nohonest,
                                                 test.data    = test.data.mis,
                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                 noise.var    = data.cont.cont$noise.var,
                                                 corr.split   = data.cont.cont$corr.split,
                                                 where.split  = data.cont.cont$where.split,
                                                 dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.mis.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.ct.propsc.mis.nohonest$corr.frst.splt <- as.character(ct.propsc.mis.nohonest$frame$var[1]) == data.cont.cont$corr.split
# mean((predict(final.tree.propsc.mis.nohonest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("3")

performance.hetero.ct <- list(true.nohonest = eval.ct.propsc.true.nohonest,
                              nois.nohonest = eval.ct.propsc.nois.nohonest,
                              mis.nohonest  = eval.ct.propsc.mis.nohonest)

#####################################################################################################################
###################################### True propensity score model, no honest #######################################
#####################################################################################################################
# In help document: Unit-specific propensity scores are not supported
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X2 + X3")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.true.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
                                      data         = data.used.full.cont.cont,
                                      weights      = 1 / tmp.propsc$prop.sc,
                                      treatment    = data.used.full.cont.cont$A,
                                      split.Rule   = "tstats", 
                                      split.Honest = T, 
                                      cv.option    = "matching", 
                                      split.Bucket = F, 
                                      xval         = 5, 
                                      cp           = 0, 
                                      minsize      = 20)
cptable.propsc.true.nohonest <- ct.propsc.true.nohonest$cptable[,1][which.min(ct.propsc.true.nohonest$cptable[,4])]
final.tree.propsc.true.nohonest <- prune(ct.propsc.true.nohonest, cptable.propsc.true.nohonest)
t1 <- Sys.time()

eval.ct.propsc.true.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.true.nohonest,
                                                  test.data    = data.cont.cont$test.data,
                                                  true.trt.eff = data.cont.cont$true.trt.eff,
                                                  noise.var    = data.cont.cont$noise.var,
                                                  corr.split   = data.cont.cont$corr.split,
                                                  where.split  = data.cont.cont$where.split,
                                                  dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.true.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.ct.propsc.true.nohonest$corr.frst.splt <- as.character(ct.propsc.true.nohonest$frame$var[1]) == data.cont.cont$corr.split
# mean((predict(final.tree.propsc.true.nohonest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("4")

#####################################################################################################################
#################################### Mis Func propensity score model, no honest #####################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.nois.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
                                      data         = data.used.full.cont.cont,
                                      weights      = 1 / tmp.propsc$prop.sc,
                                      treatment    = data.used.full.cont.cont$A,
                                      split.Rule   = "tstats", 
                                      cv.option    = "matching", 
                                      split.Honest = T, 
                                      split.Bucket = F, 
                                      xval         = 5, 
                                      cp           = 0, 
                                      minsize      = 20)

cptable.propsc.nois.nohonest <- ct.propsc.nois.nohonest$cptable[,1][which.min(ct.propsc.nois.nohonest$cptable[,4])]
final.tree.propsc.nois.nohonest <- prune(ct.propsc.nois.nohonest, cptable.propsc.nois.nohonest)
t1 <- Sys.time()

eval.ct.propsc.nois.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.nois.nohonest,
                                                  test.data    = data.cont.cont$test.data,
                                                  true.trt.eff = data.cont.cont$true.trt.eff,
                                                  noise.var    = data.cont.cont$noise.var,
                                                  corr.split   = data.cont.cont$corr.split,
                                                  where.split  = data.cont.cont$where.split,
                                                  dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.nois.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.ct.propsc.nois.nohonest$corr.frst.splt <- as.character(ct.propsc.nois.nohonest$frame$var[1]) == data.cont.cont$corr.split
# mean((predict(final.tree.propsc.nois.nohonest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("5")

#####################################################################################################################
################################# Unmeasured cov propensity score model, no honest ##################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont.mis[, !colnames(data.used.full.cont.cont.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.mis.nohonest <- causalTree(Y ~ X1 + X3 + X4 + X5 + X6, 
                                     data         = data.used.full.cont.cont.mis,
                                     weights      = 1 / tmp.propsc$prop.sc,
                                     treatment    = data.used.full.cont.cont.mis$A,
                                     split.Rule   = "tstats", 
                                     cv.option    = "matching", 
                                     split.Honest = T, 
                                     split.Bucket = F, 
                                     xval         = 5, 
                                     cp           = 0, 
                                     minsize      = 20)

cptable.propsc.mis.nohonest <- ct.propsc.mis.nohonest$cptable[,1][which.min(ct.propsc.mis.nohonest$cptable[,4])]
final.tree.propsc.mis.nohonest <- prune(ct.propsc.mis.nohonest, cptable.propsc.mis.nohonest)
t1 <- Sys.time()

eval.ct.propsc.mis.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.nohonest,
                                                 test.data    = test.data.mis,
                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                 noise.var    = data.cont.cont$noise.var,
                                                 corr.split   = data.cont.cont$corr.split,
                                                 where.split  = data.cont.cont$where.split,
                                                 dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.mis.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.ct.propsc.mis.nohonest$corr.frst.splt <- as.character(ct.propsc.mis.nohonest$frame$var[1]) == data.cont.cont$corr.split
# mean((predict(final.tree.propsc.mis.nohonest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("6")

performance.hetero.ct.final <- list(true.nohonest = eval.ct.propsc.true.nohonest,
                                    nois.nohonest = eval.ct.propsc.nois.nohonest,
                                    mis.nohonest  = eval.ct.propsc.mis.nohonest)



#####################################################################################################################
################################################## Homogeneous ######################################################
#####################################################################################################################
data.cont.cont            <- makeData.cont.noeff.cont(N             = N.training, 
                                                      n.test        = N.testing, 
                                                      coeff.prop.sc = coeff.prop.sc)
data.used.full.cont.cont  <- data.cont.cont$data.used
data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
# val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ]  

# Pretend X2 is unmeasured for unmeasured cov
data.used.full.cont.cont.mis <- data.used.full.cont.cont %>%
  select(-X2)
test.data.mis <- data.cont.cont$test.data %>%
  select(-X2)

#####################################################################################################################
###################################### True propensity score model, no honest #######################################
#####################################################################################################################
# In help document: Unit-specific propensity scores are not supported
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X2 + X3")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.true.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
                                      data         = data.used.full.cont.cont,
                                      weights      = 1 / tmp.propsc$prop.sc,
                                      treatment    = data.used.full.cont.cont$A,
                                      split.Rule   = "CT", 
                                      cv.option    = "CT", 
                                      split.Honest = T, 
                                      cv.Honest    = T, 
                                      split.Bucket = F, 
                                      xval         = 5, 
                                      cp           = 0, 
                                      minsize      = 20)
cptable.propsc.true.nohonest <- ct.propsc.true.nohonest$cptable[,1][which.min(ct.propsc.true.nohonest$cptable[,4])]
final.tree.propsc.true.nohonest <- prune(ct.propsc.true.nohonest, cptable.propsc.true.nohonest)
t1 <- Sys.time()

eval.ct.propsc.true.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.true.nohonest,
                                                  test.data    = data.cont.cont$test.data,
                                                  true.trt.eff = data.cont.cont$true.trt.eff,
                                                  noise.var    = data.cont.cont$noise.var,
                                                  corr.split   = data.cont.cont$corr.split,
                                                  where.split  = data.cont.cont$where.split,
                                                  dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.true.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
# mean((predict(final.tree.propsc.true.nohonest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("7")

#####################################################################################################################
#################################### Mis Func propensity score model, no honest #####################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.nois.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
                                      data       = data.used.full.cont.cont,
                                      weights    = 1 / tmp.propsc$prop.sc,
                                      treatment  = data.used.full.cont.cont$A,
                                      split.Rule = "CT", 
                                      cv.option  = "CT", 
                                      split.Honest = T, 
                                      cv.Honest = T, 
                                      split.Bucket = F, 
                                      xval = 5, 
                                      cp = 0, 
                                      minsize = 20)

cptable.propsc.nois.nohonest <- ct.propsc.nois.nohonest$cptable[,1][which.min(ct.propsc.nois.nohonest$cptable[,4])]
final.tree.propsc.nois.nohonest <- prune(ct.propsc.nois.nohonest, cptable.propsc.nois.nohonest)
t1 <- Sys.time()

eval.ct.propsc.nois.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.nois.nohonest,
                                                  test.data    = data.cont.cont$test.data,
                                                  true.trt.eff = data.cont.cont$true.trt.eff,
                                                  noise.var    = data.cont.cont$noise.var,
                                                  corr.split   = data.cont.cont$corr.split,
                                                  where.split  = data.cont.cont$where.split,
                                                  dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.nois.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
# mean((predict(final.tree.propsc.nois.nohonest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("8")

#####################################################################################################################
################################ Unmeasured cov propensity score model, no honest ###################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont.mis[, !colnames(data.used.full.cont.cont.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.mis.nohonest <- causalTree(Y ~ X1 + X3 + X4 + X5 + X6, 
                                     data       = data.used.full.cont.cont.mis,
                                     weights    = 1 / tmp.propsc$prop.sc,
                                     treatment  = data.used.full.cont.cont.mis$A,
                                     split.Rule = "CT", 
                                     cv.option  = "CT", 
                                     split.Honest = T, 
                                     cv.Honest = T, 
                                     split.Bucket = F, 
                                     xval = 5, 
                                     cp = 0, 
                                     minsize = 20)

cptable.propsc.mis.nohonest <- ct.propsc.mis.nohonest$cptable[,1][which.min(ct.propsc.mis.nohonest$cptable[,4])]
final.tree.propsc.mis.nohonest <- prune(ct.propsc.mis.nohonest, cptable.propsc.mis.nohonest)
t1 <- Sys.time()

eval.ct.propsc.mis.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.nohonest,
                                                 test.data    = test.data.mis,
                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                 noise.var    = data.cont.cont$noise.var,
                                                 corr.split   = data.cont.cont$corr.split,
                                                 where.split  = data.cont.cont$where.split,
                                                 dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.mis.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
# mean((predict(final.tree.propsc.mis.nohonest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("9")

performance.homo.ct <- list(true.nohonest = eval.ct.propsc.true.nohonest,
                            nois.nohonest = eval.ct.propsc.nois.nohonest,
                            mis.nohonest  = eval.ct.propsc.mis.nohonest)

#####################################################################################################################
###################################### True propensity score model, no honest #######################################
#####################################################################################################################
# In help document: Unit-specific propensity scores are not supported
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X2 + X3")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.true.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
                                      data         = data.used.full.cont.cont,
                                      weights      = 1 / tmp.propsc$prop.sc,
                                      treatment    = data.used.full.cont.cont$A,
                                      split.Rule   = "tstats", 
                                      cv.option    = "matching", 
                                      split.Honest = T, 
                                      split.Bucket = F, 
                                      xval         = 5, 
                                      cp           = 0, 
                                      minsize      = 20)
cptable.propsc.true.nohonest <- ct.propsc.true.nohonest$cptable[,1][which.min(ct.propsc.true.nohonest$cptable[,4])]
final.tree.propsc.true.nohonest <- prune(ct.propsc.true.nohonest, cptable.propsc.true.nohonest)
t1 <- Sys.time()

eval.ct.propsc.true.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.true.nohonest,
                                                  test.data    = data.cont.cont$test.data,
                                                  true.trt.eff = data.cont.cont$true.trt.eff,
                                                  noise.var    = data.cont.cont$noise.var,
                                                  corr.split   = data.cont.cont$corr.split,
                                                  where.split  = data.cont.cont$where.split,
                                                  dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.true.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
# mean((predict(final.tree.propsc.true.nohonest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("10")

#####################################################################################################################
##################################### Mis Func propensity score model, no honest ####################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.nois.nohonest <- causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6, 
                                      data         = data.used.full.cont.cont,
                                      weights      = 1 / tmp.propsc$prop.sc,
                                      treatment    = data.used.full.cont.cont$A,
                                      split.Rule   = "tstats", 
                                      cv.option    = "matching", 
                                      split.Honest = T, 
                                      split.Bucket = F, 
                                      xval         = 5, 
                                      cp           = 0, 
                                      minsize      = 20)

cptable.propsc.nois.nohonest <- ct.propsc.nois.nohonest$cptable[,1][which.min(ct.propsc.nois.nohonest$cptable[,4])]
final.tree.propsc.nois.nohonest <- prune(ct.propsc.nois.nohonest, cptable.propsc.nois.nohonest)
t1 <- Sys.time()

eval.ct.propsc.nois.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.nois.nohonest,
                                                  test.data    = data.cont.cont$test.data,
                                                  true.trt.eff = data.cont.cont$true.trt.eff,
                                                  noise.var    = data.cont.cont$noise.var,
                                                  corr.split   = data.cont.cont$corr.split,
                                                  where.split  = data.cont.cont$where.split,
                                                  dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.nois.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
# mean((predict(final.tree.propsc.nois.nohonest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("11")

#####################################################################################################################
################################## Unmeasured cov propensity score model, no honest #################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont.mis[, !colnames(data.used.full.cont.cont.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

ct.propsc.mis.nohonest <- causalTree(Y ~ X1 + X3 + X4 + X5 + X6, 
                                     data         = data.used.full.cont.cont.mis,
                                     weights      = 1 / tmp.propsc$prop.sc,
                                     treatment    = data.used.full.cont.cont.mis$A,
                                     split.Rule   = "tstats", 
                                     cv.option    = "matching", 
                                     split.Honest = T, 
                                     split.Bucket = F, 
                                     xval         = 5, 
                                     cp           = 0, 
                                     minsize      = 20)

cptable.propsc.mis.nohonest <- ct.propsc.mis.nohonest$cptable[,1][which.min(ct.propsc.mis.nohonest$cptable[,4])]
final.tree.propsc.mis.nohonest <- prune(ct.propsc.mis.nohonest, cptable.propsc.mis.nohonest)
t1 <- Sys.time()

eval.ct.propsc.mis.nohonest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.nohonest,
                                                 test.data    = test.data.mis,
                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                 noise.var    = data.cont.cont$noise.var,
                                                 corr.split   = data.cont.cont$corr.split,
                                                 where.split  = data.cont.cont$where.split,
                                                 dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.mis.nohonest$t <- as.numeric(difftime(t1, t0, units = "secs"))
# mean((predict(final.tree.propsc.mis.nohonest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("12")

performance.homo.ct.final <- list(true.nohonest = eval.ct.propsc.true.nohonest,
                                  nois.nohonest = eval.ct.propsc.nois.nohonest,
                                  mis.nohonest  = eval.ct.propsc.mis.nohonest)

file.name = paste("../Data/MainSimulationsCt/", toString(job.number), ".RData", sep = "")
save(performance.hetero.ct, performance.hetero.ct.final, 
     performance.homo.ct, performance.homo.ct.final,
     file = file.name)
