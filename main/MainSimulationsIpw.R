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

setwd("main/")

#####################################################################################################################
################################################# Heterogeneous #####################################################
#####################################################################################################################
data.cont.cont            <- makeData.cont.eff.cont(1000, 1000, coeff.prop.sc = 0.6)
data.used.full.cont.cont  <- data.cont.cont$data.used
data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
# val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ]  

#####################################################################################################################
######################### 1. ipw: GLM Model, inside node, True propensity score model, cv1 ##########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.true.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                              est.used          = "IPW",
                                                              type.var          = "cont",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM", 
                                                              propsc.form.true  = "A ~ X1 + X2 + X3",
                                                              w                 = NULL,
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.true.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                              tree.list         = seq.created.estipw.glm.propscinnd.true.cv1$tree.list, 
                                                              lambda.used       = qchisq(0.95, 1), 
                                                              val.sample        = data.validation.cont.cont, 
                                                              type.var          = "cont",
                                                              propsc.mod.out    = F, 
                                                              propsc.mthd       = "GLM", 
                                                              propsc.form.true  = "A ~ X1 + X2 + X3",
                                                              val.w             = NULL,
                                                              propsc.mod.insplt = F, 
                                                              min.obs.mod       = 5)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.true.cv1[[1]], 
                                                               test.data    = data.cont.cont$test.data,
                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                               noise.var    = data.cont.cont$noise.var,
                                                               corr.split   = data.cont.cont$corr.split,
                                                               where.split  = data.cont.cont$where.split,
                                                               dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinnd.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinnd.true.cv1$corr.frst.splt <- as.character(seq.created.estipw.glm.propscinnd.true.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("1")

#####################################################################################################################
######################### 2. ipw: GLM Model, inside node, Noisy propensity score model, cv1 #########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.nois.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                              est.used          = "IPW",
                                                              type.var          = "cont",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM", 
                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)", 
                                                              w                 = NULL, 
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.nois.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                              tree.list         = seq.created.estipw.glm.propscinnd.nois.cv1$tree.list, 
                                                              lambda.used       = qchisq(0.95, 1), 
                                                              val.sample        = data.validation.cont.cont, 
                                                              type.var          = "cont",
                                                              propsc.mod.out    = F, 
                                                              propsc.mthd       = "GLM", 
                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)", 
                                                              val.w             = NULL,
                                                              propsc.mod.insplt = F, 
                                                              min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.nois.cv1[[1]], 
                                                               test.data    = data.cont.cont$test.data,
                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                               noise.var    = data.cont.cont$noise.var,
                                                               corr.split   = data.cont.cont$corr.split,
                                                               where.split  = data.cont.cont$where.split,
                                                               dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinnd.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinnd.nois.cv1$corr.frst.splt <- as.character(seq.created.estipw.glm.propscinnd.nois.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("2")

#####################################################################################################################
##################### 3. ipw: GLM Model, inside node, Misspecified propensity score model, cv1 ######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                             est.used          = "IPW",
                                                             type.var          = "cont",
                                                             propsc.mod.out    = F,
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6", 
                                                             w                 = NULL, 
                                                             propsc.mod.insplt = F,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscinnd.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                             tree.list         = seq.created.estipw.glm.propscinnd.mis.cv1$tree.list, 
                                                             lambda.used       = qchisq(0.95, 1), 
                                                             val.sample        = data.validation.cont.cont, 
                                                             type.var          = "cont",
                                                             propsc.mod.out    = F, 
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6", 
                                                             val.w             = NULL,
                                                             propsc.mod.insplt = F, 
                                                             min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.mis.cv1[[1]], 
                                                              test.data    = data.cont.cont$test.data,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split,
                                                              where.split  = data.cont.cont$where.split,
                                                              dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinnd.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinnd.mis.cv1$corr.frst.splt <- as.character(seq.created.estipw.glm.propscinnd.mis.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("3")

performance.hetero.ipw <- list(glm.propscinnd.true.cv1 = eval.final.estipw.glm.propscinnd.true.cv1,
                               glm.propscinnd.nois.cv1 = eval.final.estipw.glm.propscinnd.nois.cv1,
                               glm.propscinnd.mis.cv1  = eval.final.estipw.glm.propscinnd.mis.cv1)



#####################################################################################################################
################################################### Homogeneous #####################################################
#####################################################################################################################
data.cont.cont            <- makeData.cont.noeff.cont(1000, 1000, coeff.prop.sc = 0.6)
data.used.full.cont.cont  <- data.cont.cont$data.used
data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
# val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ]  

#####################################################################################################################
######################### 7. ipw: GLM Model, inside node, True propensity score model, cv1 ##########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.true.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                              est.used          = "IPW",
                                                              type.var          = "cont",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM", 
                                                              propsc.form.true  = "A ~ X1 + X2 + X3",
                                                              w                 = NULL,
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.true.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                              tree.list         = seq.created.estipw.glm.propscinnd.true.cv1$tree.list, 
                                                              lambda.used       = qchisq(0.95, 1), 
                                                              val.sample        = data.validation.cont.cont, 
                                                              type.var          = "cont",
                                                              propsc.mod.out    = F, 
                                                              propsc.mthd       = "GLM", 
                                                              propsc.form.true  = "A ~ X1 + X2 + X3",
                                                              val.w             = NULL,
                                                              propsc.mod.insplt = F, 
                                                              min.obs.mod       = 5)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.true.cv1[[1]], 
                                                               test.data    = data.cont.cont$test.data,
                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                               noise.var    = data.cont.cont$noise.var,
                                                               corr.split   = data.cont.cont$corr.split,
                                                               where.split  = data.cont.cont$where.split,
                                                               dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinnd.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("7")

#####################################################################################################################
######################### 8. ipw: GLM Model, inside node, Noisy propensity score model, cv1 #########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.nois.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                              est.used          = "IPW",
                                                              type.var          = "cont",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM", 
                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)", 
                                                              w                 = NULL, 
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.nois.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                              tree.list         = seq.created.estipw.glm.propscinnd.nois.cv1$tree.list, 
                                                              lambda.used       = qchisq(0.95, 1), 
                                                              val.sample        = data.validation.cont.cont, 
                                                              type.var          = "cont",
                                                              propsc.mod.out    = F, 
                                                              propsc.mthd       = "GLM", 
                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)", 
                                                              val.w             = NULL,
                                                              propsc.mod.insplt = F, 
                                                              min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.nois.cv1[[1]], 
                                                               test.data    = data.cont.cont$test.data,
                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                               noise.var    = data.cont.cont$noise.var,
                                                               corr.split   = data.cont.cont$corr.split,
                                                               where.split  = data.cont.cont$where.split,
                                                               dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinnd.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("8")

#####################################################################################################################
##################### 9. ipw: GLM Model, inside node, Misspecified propensity score model, cv1 ######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                             est.used          = "IPW",
                                                             type.var          = "cont",
                                                             propsc.mod.out    = F,
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6", 
                                                             w                 = NULL, 
                                                             propsc.mod.insplt = F,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscinnd.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                             tree.list         = seq.created.estipw.glm.propscinnd.mis.cv1$tree.list, 
                                                             lambda.used       = qchisq(0.95, 1), 
                                                             val.sample        = data.validation.cont.cont, 
                                                             type.var          = "cont",
                                                             propsc.mod.out    = F, 
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6", 
                                                             val.w             = NULL,
                                                             propsc.mod.insplt = F, 
                                                             min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.mis.cv1[[1]], 
                                                              test.data    = data.cont.cont$test.data,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split,
                                                              where.split  = data.cont.cont$where.split,
                                                              dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinnd.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("9")

performance.homo.ipw <- list(glm.propscinnd.true.cv1 = eval.final.estipw.glm.propscinnd.true.cv1,
                             glm.propscinnd.nois.cv1 = eval.final.estipw.glm.propscinnd.nois.cv1,
                             glm.propscinnd.mis.cv1  = eval.final.estipw.glm.propscinnd.mis.cv1)









