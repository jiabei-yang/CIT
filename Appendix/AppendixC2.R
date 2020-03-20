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

# Pretend X2 is unmeasured for unmeasured cov
data.used.full.cont.cont.mis <- data.used.full.cont.cont %>%
  select(-X2)
test.data.mis <- data.cont.cont$test.data %>%
  select(-X2)

#####################################################################################################################
################################### True propensity score model, honest Estimation ##################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X2 + X3")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.cont.cont$A == 1)
ctrlIdx <- which(data.used.full.cont.cont$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.cont.cont[train.idx, ]
est.data   <- data.used.full.cont.cont[-train.idx, ]

ct.propsc.true.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                           data = train.data,
                                           weights = 1 / tmp.propsc$prop.sc[train.idx],
                                           treatment = train.data$A,
                                           est_data = est.data,
                                           est_weights = 1 / tmp.propsc$prop.sc[-train.idx],
                                           est_treatment = est.data$A,
                                           split.Rule = "CT", 
                                           split.Honest = T,
                                           HonestSampleSize = nrow(est.data),
                                           split.Bucket = F,
                                           cv.option = "CT",
                                           cv.Honest = T)
cptable.propsc.true.honest <- ct.propsc.true.honest$cptable[,1][which.min(ct.propsc.true.honest$cptable[,4])]
final.tree.propsc.true.honest <- prune(ct.propsc.true.honest, cptable.propsc.true.honest)
t1 <- Sys.time()

eval.ct.propsc.true.honest <- eval.measures.eff(final.tree   = final.tree.propsc.true.honest,
                                                test.data    = data.cont.cont$test.data,
                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                noise.var    = data.cont.cont$noise.var,
                                                corr.split   = data.cont.cont$corr.split,
                                                where.split  = data.cont.cont$where.split,
                                                dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.true.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.ct.propsc.true.honest$corr.frst.splt <- as.character(ct.propsc.true.honest$frame$var[1]) == data.cont.cont$corr.split
# mean((predict(final.tree.propsc.true.honest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("1")

#####################################################################################################################
################################# Mis func propensity score model, honest estimation ################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.cont.cont$A == 1)
ctrlIdx <- which(data.used.full.cont.cont$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.cont.cont[train.idx, ]
est.data   <- data.used.full.cont.cont[-train.idx, ]

ct.propsc.nois.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                           data = train.data,
                                           weights = 1 / tmp.propsc$prop.sc[train.idx],
                                           treatment = train.data$A,
                                           est_data = est.data,
                                           est_weights = 1 / tmp.propsc$prop.sc[-train.idx],
                                           est_treatment = est.data$A,
                                           split.Rule = "CT", 
                                           split.Honest = T,
                                           HonestSampleSize = nrow(est.data),
                                           split.Bucket = F,
                                           cv.option = "CT",
                                           cv.Honest = T)
cptable.propsc.nois.honest <- ct.propsc.nois.honest$cptable[,1][which.min(ct.propsc.nois.honest$cptable[,4])]
final.tree.propsc.nois.honest <- prune(ct.propsc.nois.honest, cptable.propsc.nois.honest)
t1 <- Sys.time()

eval.ct.propsc.nois.honest <- eval.measures.eff(final.tree   = final.tree.propsc.nois.honest,
                                                test.data    = data.cont.cont$test.data,
                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                noise.var    = data.cont.cont$noise.var,
                                                corr.split   = data.cont.cont$corr.split,
                                                where.split  = data.cont.cont$where.split,
                                                dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.nois.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.ct.propsc.nois.honest$corr.frst.splt <- as.character(ct.propsc.nois.honest$frame$var[1]) == data.cont.cont$corr.split
# mse.propsc.nois.honest <- mean((predict(final.tree.propsc.nois.honest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("2")

#####################################################################################################################
############################# Unmeasured cov propensity score model, honest estimation ##############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont.mis[, !colnames(data.used.full.cont.cont.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.cont.cont.mis$A == 1)
ctrlIdx <- which(data.used.full.cont.cont.mis$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.cont.cont.mis[train.idx, ]
est.data   <- data.used.full.cont.cont.mis[-train.idx, ]

ct.propsc.mis.honest <- honest.causalTree(Y ~ X1 + X3 + X4 + X5 + X6,
                                          data = train.data,
                                          weights = 1 / tmp.propsc$prop.sc[train.idx],
                                          treatment = train.data$A,
                                          est_data = est.data,
                                          est_weights = 1 / tmp.propsc$prop.sc[-train.idx],
                                          est_treatment = est.data$A,
                                          split.Rule = "CT", 
                                          split.Honest = T,
                                          HonestSampleSize = nrow(est.data),
                                          split.Bucket = F,
                                          cv.option = "CT",
                                          cv.Honest = T)
cptable.propsc.mis.honest <- ct.propsc.mis.honest$cptable[,1][which.min(ct.propsc.mis.honest$cptable[,4])]
final.tree.propsc.mis.honest <- prune(ct.propsc.mis.honest, cptable.propsc.mis.honest)
t1 <- Sys.time()

eval.ct.propsc.mis.honest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.honest,
                                               test.data    = test.data.mis,
                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                               noise.var    = data.cont.cont$noise.var,
                                               corr.split   = data.cont.cont$corr.split,
                                               where.split  = data.cont.cont$where.split,
                                               dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.mis.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.ct.propsc.mis.honest$corr.frst.splt <- as.character(ct.propsc.mis.honest$frame$var[1]) == data.cont.cont$corr.split
# mse.propsc.mis.honest <- mean((predict(final.tree.propsc.mis.honest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("3")

performance.hetero.ct <- list(true.honest   = eval.ct.propsc.true.honest,
                              nois.honest   = eval.ct.propsc.nois.honest,
                              mis.honest    = eval.ct.propsc.mis.honest)

#####################################################################################################################
################################### True propensity score model, honest Estimation ##################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X2 + X3")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.cont.cont$A == 1)
ctrlIdx <- which(data.used.full.cont.cont$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.cont.cont[train.idx, ]
est.data   <- data.used.full.cont.cont[-train.idx, ]

ct.propsc.true.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                           data             = train.data,
                                           weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                           treatment        = train.data$A,
                                           est_data         = est.data,
                                           est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                           est_treatment    = est.data$A,
                                           split.Rule       = "tstats", 
                                           split.Honest     = T,
                                           HonestSampleSize = nrow(est.data),
                                           split.Bucket     = F,
                                           cv.option        = "matching")
cptable.propsc.true.honest <- ct.propsc.true.honest$cptable[,1][which.min(ct.propsc.true.honest$cptable[,4])]
final.tree.propsc.true.honest <- prune(ct.propsc.true.honest, cptable.propsc.true.honest)
t1 <- Sys.time()

eval.ct.propsc.true.honest <- eval.measures.eff(final.tree   = final.tree.propsc.true.honest,
                                                test.data    = data.cont.cont$test.data,
                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                noise.var    = data.cont.cont$noise.var,
                                                corr.split   = data.cont.cont$corr.split,
                                                where.split  = data.cont.cont$where.split,
                                                dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.true.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.ct.propsc.true.honest$corr.frst.splt <- as.character(ct.propsc.true.honest$frame$var[1]) == data.cont.cont$corr.split
# mean((predict(final.tree.propsc.true.honest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("4")

#####################################################################################################################
################################## Mis func propensity score model, honest estimation ###############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.cont.cont$A == 1)
ctrlIdx <- which(data.used.full.cont.cont$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.cont.cont[train.idx, ]
est.data   <- data.used.full.cont.cont[-train.idx, ]

ct.propsc.nois.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                           data             = train.data,
                                           weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                           treatment        = train.data$A,
                                           est_data         = est.data,
                                           est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                           est_treatment    = est.data$A,
                                           split.Rule       = "tstats", 
                                           split.Honest     = T,
                                           HonestSampleSize = nrow(est.data),
                                           split.Bucket     = F,
                                           cv.option        = "matching")
cptable.propsc.nois.honest <- ct.propsc.nois.honest$cptable[,1][which.min(ct.propsc.nois.honest$cptable[,4])]
final.tree.propsc.nois.honest <- prune(ct.propsc.nois.honest, cptable.propsc.nois.honest)
t1 <- Sys.time()

eval.ct.propsc.nois.honest <- eval.measures.eff(final.tree   = final.tree.propsc.nois.honest,
                                                test.data    = data.cont.cont$test.data,
                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                noise.var    = data.cont.cont$noise.var,
                                                corr.split   = data.cont.cont$corr.split,
                                                where.split  = data.cont.cont$where.split,
                                                dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.nois.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.ct.propsc.nois.honest$corr.frst.splt <- as.character(ct.propsc.nois.honest$frame$var[1]) == data.cont.cont$corr.split
# mse.propsc.nois.honest <- mean((predict(final.tree.propsc.nois.honest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("5")

#####################################################################################################################
############################ Unmeasured cov propensity score model, honest estimation ###############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont.mis[, !colnames(data.used.full.cont.cont.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.cont.cont.mis$A == 1)
ctrlIdx <- which(data.used.full.cont.cont.mis$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.cont.cont.mis[train.idx, ]
est.data   <- data.used.full.cont.cont.mis[-train.idx, ]

ct.propsc.mis.honest <- honest.causalTree(Y ~ X1 + X3 + X4 + X5 + X6,
                                          data             = train.data,
                                          weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                          treatment        = train.data$A,
                                          est_data         = est.data,
                                          est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                          est_treatment    = est.data$A,
                                          split.Rule       = "tstats", 
                                          split.Honest     = T,
                                          HonestSampleSize = nrow(est.data),
                                          split.Bucket     = F,
                                          cv.option        = "matching")
cptable.propsc.mis.honest <- ct.propsc.mis.honest$cptable[,1][which.min(ct.propsc.mis.honest$cptable[,4])]
final.tree.propsc.mis.honest <- prune(ct.propsc.mis.honest, cptable.propsc.mis.honest)
t1 <- Sys.time()

eval.ct.propsc.mis.honest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.honest,
                                               test.data    = test.data.mis,
                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                               noise.var    = data.cont.cont$noise.var,
                                               corr.split   = data.cont.cont$corr.split,
                                               where.split  = data.cont.cont$where.split,
                                               dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.mis.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.ct.propsc.mis.honest$corr.frst.splt <- as.character(ct.propsc.mis.honest$frame$var[1]) == data.cont.cont$corr.split
# mse.propsc.mis.honest <- mean((predict(final.tree.propsc.mis.honest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("6")

performance.hetero.ct.final <- list(true.honest   = eval.ct.propsc.true.honest,
                                    nois.honest   = eval.ct.propsc.nois.honest,
                                    mis.honest    = eval.ct.propsc.mis.honest)



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
################################### True propensity score model, honest Estimation ##################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X2 + X3")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.cont.cont$A == 1)
ctrlIdx <- which(data.used.full.cont.cont$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.cont.cont[train.idx, ]
est.data   <- data.used.full.cont.cont[-train.idx, ]

ct.propsc.true.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                           data = train.data,
                                           weights = 1 / tmp.propsc$prop.sc[train.idx],
                                           treatment = train.data$A,
                                           est_data = est.data,
                                           est_weights = 1 / tmp.propsc$prop.sc[-train.idx],
                                           est_treatment = est.data$A,
                                           split.Rule = "CT", 
                                           split.Honest = T,
                                           HonestSampleSize = nrow(est.data),
                                           split.Bucket = F,
                                           cv.option = "CT",
                                           cv.Honest = T)
cptable.propsc.true.honest <- ct.propsc.true.honest$cptable[,1][which.min(ct.propsc.true.honest$cptable[,4])]
final.tree.propsc.true.honest <- prune(ct.propsc.true.honest, cptable.propsc.true.honest)
t1 <- Sys.time()

eval.ct.propsc.true.honest <- eval.measures.eff(final.tree   = final.tree.propsc.true.honest,
                                                test.data    = data.cont.cont$test.data,
                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                noise.var    = data.cont.cont$noise.var,
                                                corr.split   = data.cont.cont$corr.split,
                                                where.split  = data.cont.cont$where.split,
                                                dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.true.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
# mean((predict(final.tree.propsc.true.honest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("7")

#####################################################################################################################
################################ Mis func propensity score model, honest estimation #################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.cont.cont$A == 1)
ctrlIdx <- which(data.used.full.cont.cont$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.cont.cont[train.idx, ]
est.data   <- data.used.full.cont.cont[-train.idx, ]

ct.propsc.nois.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                           data = train.data,
                                           weights = 1 / tmp.propsc$prop.sc[train.idx],
                                           treatment = train.data$A,
                                           est_data = est.data,
                                           est_weights = 1 / tmp.propsc$prop.sc[-train.idx],
                                           est_treatment = est.data$A,
                                           split.Rule = "CT", 
                                           split.Honest = T,
                                           HonestSampleSize = nrow(est.data),
                                           split.Bucket = F,
                                           cv.option = "CT",
                                           cv.Honest = T)
cptable.propsc.nois.honest <- ct.propsc.nois.honest$cptable[,1][which.min(ct.propsc.nois.honest$cptable[,4])]
final.tree.propsc.nois.honest <- prune(ct.propsc.nois.honest, cptable.propsc.nois.honest)
t1 <- Sys.time()

eval.ct.propsc.nois.honest <- eval.measures.eff(final.tree   = final.tree.propsc.nois.honest,
                                                test.data    = data.cont.cont$test.data,
                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                noise.var    = data.cont.cont$noise.var,
                                                corr.split   = data.cont.cont$corr.split,
                                                where.split  = data.cont.cont$where.split,
                                                dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.nois.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
# mse.propsc.nois.honest <- mean((predict(final.tree.propsc.nois.honest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("8")

#####################################################################################################################
############################# Unmeasured cov propensity score model, honest estimation ##############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont.mis[, !colnames(data.used.full.cont.cont.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.cont.cont.mis$A == 1)
ctrlIdx <- which(data.used.full.cont.cont.mis$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.cont.cont.mis[train.idx, ]
est.data   <- data.used.full.cont.cont.mis[-train.idx, ]

ct.propsc.mis.honest <- honest.causalTree(Y ~ X1 + X3 + X4 + X5 + X6,
                                          data = train.data,
                                          weights = 1 / tmp.propsc$prop.sc[train.idx],
                                          treatment = train.data$A,
                                          est_data = est.data,
                                          est_weights = 1 / tmp.propsc$prop.sc[-train.idx],
                                          est_treatment = est.data$A,
                                          split.Rule = "CT", 
                                          split.Honest = T,
                                          HonestSampleSize = nrow(est.data),
                                          split.Bucket = F,
                                          cv.option = "CT",
                                          cv.Honest = T)
cptable.propsc.mis.honest <- ct.propsc.mis.honest$cptable[,1][which.min(ct.propsc.mis.honest$cptable[,4])]
final.tree.propsc.mis.honest <- prune(ct.propsc.mis.honest, cptable.propsc.mis.honest)
t1 <- Sys.time()

eval.ct.propsc.mis.honest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.honest,
                                               test.data    = test.data.mis,
                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                               noise.var    = data.cont.cont$noise.var,
                                               corr.split   = data.cont.cont$corr.split,
                                               where.split  = data.cont.cont$where.split,
                                               dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.mis.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
# mse.propsc.mis.honest <- mean((predict(final.tree.propsc.mis.honest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("9")

performance.homo.ct <- list(true.honest   = eval.ct.propsc.true.honest,
                            nois.honest   = eval.ct.propsc.nois.honest,
                            mis.honest    = eval.ct.propsc.mis.honest)

#####################################################################################################################
################################### True propensity score model, honest Estimation ##################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X2 + X3")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.cont.cont$A == 1)
ctrlIdx <- which(data.used.full.cont.cont$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.cont.cont[train.idx, ]
est.data   <- data.used.full.cont.cont[-train.idx, ]

ct.propsc.true.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                           data             = train.data,
                                           weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                           treatment        = train.data$A,
                                           est_data         = est.data,
                                           est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                           est_treatment    = est.data$A,
                                           split.Rule       = "tstats", 
                                           split.Honest     = T,
                                           HonestSampleSize = nrow(est.data),
                                           split.Bucket     = F,
                                           cv.option        = "matching")
cptable.propsc.true.honest <- ct.propsc.true.honest$cptable[,1][which.min(ct.propsc.true.honest$cptable[,4])]
final.tree.propsc.true.honest <- prune(ct.propsc.true.honest, cptable.propsc.true.honest)
t1 <- Sys.time()

eval.ct.propsc.true.honest <- eval.measures.eff(final.tree   = final.tree.propsc.true.honest,
                                                test.data    = data.cont.cont$test.data,
                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                noise.var    = data.cont.cont$noise.var,
                                                corr.split   = data.cont.cont$corr.split,
                                                where.split  = data.cont.cont$where.split,
                                                dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.true.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
# mean((predict(final.tree.propsc.true.honest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("10")

#####################################################################################################################
################################ Mis func propensity score model, honest estimation #################################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont[, !colnames(data.used.full.cont.cont) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.cont.cont$A == 1)
ctrlIdx <- which(data.used.full.cont.cont$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.cont.cont[train.idx, ]
est.data   <- data.used.full.cont.cont[-train.idx, ]

ct.propsc.nois.honest <- honest.causalTree(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                                           data             = train.data,
                                           weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                           treatment        = train.data$A,
                                           est_data         = est.data,
                                           est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                           est_treatment    = est.data$A,
                                           split.Rule       = "tstats", 
                                           split.Honest     = T,
                                           HonestSampleSize = nrow(est.data),
                                           split.Bucket     = F,
                                           cv.option        = "matching")
cptable.propsc.nois.honest <- ct.propsc.nois.honest$cptable[,1][which.min(ct.propsc.nois.honest$cptable[,4])]
final.tree.propsc.nois.honest <- prune(ct.propsc.nois.honest, cptable.propsc.nois.honest)
t1 <- Sys.time()

eval.ct.propsc.nois.honest <- eval.measures.eff(final.tree   = final.tree.propsc.nois.honest,
                                                test.data    = data.cont.cont$test.data,
                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                noise.var    = data.cont.cont$noise.var,
                                                corr.split   = data.cont.cont$corr.split,
                                                where.split  = data.cont.cont$where.split,
                                                dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.nois.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
# mse.propsc.nois.honest <- mean((predict(final.tree.propsc.nois.honest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("11")

#####################################################################################################################
############################# Unmeasured cov propensity score model, honest estimation ##############################
#####################################################################################################################
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = data.used.full.cont.cont.mis[, !colnames(data.used.full.cont.cont.mis) %in% c("Y")],
                          method    = "GLM",
                          form.true = "A ~ X1 + X3 + X4 + X5 + X6")
tmp.propsc$prop.sc <- ifelse(data.used.full.cont.cont.mis$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

trtIdx  <- which(data.used.full.cont.cont.mis$A == 1)
ctrlIdx <- which(data.used.full.cont.cont.mis$A == 0)
train.idx <- c(sample(trtIdx, length(trtIdx) / 2),
               sample(ctrlIdx, length(ctrlIdx) / 2))
train.data <- data.used.full.cont.cont.mis[train.idx, ]
est.data   <- data.used.full.cont.cont.mis[-train.idx, ]

ct.propsc.mis.honest <- honest.causalTree(Y ~ X1 + X3 + X4 + X5 + X6,
                                          data             = train.data,
                                          weights          = 1 / tmp.propsc$prop.sc[train.idx],
                                          treatment        = train.data$A,
                                          est_data         = est.data,
                                          est_weights      = 1 / tmp.propsc$prop.sc[-train.idx],
                                          est_treatment    = est.data$A,
                                          split.Rule       = "tstats", 
                                          split.Honest     = T,
                                          HonestSampleSize = nrow(est.data),
                                          split.Bucket     = F,
                                          cv.option        = "matching")
cptable.propsc.mis.honest <- ct.propsc.mis.honest$cptable[,1][which.min(ct.propsc.mis.honest$cptable[,4])]
final.tree.propsc.mis.honest <- prune(ct.propsc.mis.honest, cptable.propsc.mis.honest)
t1 <- Sys.time()

eval.ct.propsc.mis.honest <- eval.measures.eff(final.tree   = final.tree.propsc.mis.honest,
                                               test.data    = test.data.mis,
                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                               noise.var    = data.cont.cont$noise.var,
                                               corr.split   = data.cont.cont$corr.split,
                                               where.split  = data.cont.cont$where.split,
                                               dir.split    = data.cont.cont$dir.split)
eval.ct.propsc.mis.honest$t <- as.numeric(difftime(t1, t0, units = "secs"))
# mse.propsc.mis.honest <- mean((predict(final.tree.propsc.mis.honest, data.cont.cont$test.data) - data.cont.cont$true.trt.eff)^2)
print("12")

performance.homo.ct.final <- list(true.honest   = eval.ct.propsc.true.honest,
                                  nois.honest   = eval.ct.propsc.nois.honest,
                                  mis.honest    = eval.ct.propsc.mis.honest)

file.name = paste("../Data/AppendixC2/", toString(job.number), ".RData", sep = "")
save(performance.homo.ct, performance.hetero.ct,
     performance.homo.ct.final, performance.hetero.ct.final,
     file = file.name)


