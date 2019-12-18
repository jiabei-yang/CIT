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

#####################################################################################################################
################################################# Heterogeneous #####################################################
#####################################################################################################################
data.cont.cont            <- makeData.cont.eff.cont(1000, 1000, coeff.prop.sc = 0.6)
data.used.full.cont.cont  <- data.cont.cont$data.used
data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
# val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ]  

#####################################################################################################################
######################## 1. ipw: GLM Model, inside split, True propensity score model, cv1 ##########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.true.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                est.used          = "IPW",
                                                                type.var          = "cont",
                                                                propsc.mod.out    = F,
                                                                propsc.mthd       = "GLM",
                                                                propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                w                 = NULL,
                                                                propsc.mod.insplt = T,
                                                                num.truc.obs      = 30,
                                                                min.node          = 20)

final.tree.estipw.glm.propscinsplt.true.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont,
                                                                tree.list         = seq.created.estipw.glm.propscinsplt.true.cv1$tree.list,
                                                                lambda.used       = qchisq(0.95, 1),
                                                                val.sample        = data.validation.cont.cont,
                                                                type.var          = "cont",
                                                                propsc.mod.out    = F,
                                                                propsc.mthd       = "GLM",
                                                                propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                val.w             = NULL,
                                                                propsc.mod.insplt = T,
                                                                min.obs.mod       = 5)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.true.cv1[[1]],
                                                                 test.data    = data.cont.cont$test.data,
                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                 noise.var    = data.cont.cont$noise.var,
                                                                 corr.split   = data.cont.cont$corr.split,
                                                                 where.split  = data.cont.cont$where.split,
                                                                 dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinsplt.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinsplt.true.cv1$corr.frst.splt <- as.character(seq.created.estipw.glm.propscinsplt.true.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("1")

#####################################################################################################################
###################### 2. ipw: GLM Model, inside split, Mis Func propensity score model, cv1 ########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.nois.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                est.used          = "IPW",
                                                                type.var          = "cont",
                                                                propsc.mod.out    = F,
                                                                propsc.mthd       = "GLM",
                                                                propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                w                 = NULL,
                                                                propsc.mod.insplt = T,
                                                                num.truc.obs      = 30,
                                                                min.node          = 20)

final.tree.estipw.glm.propscinsplt.nois.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont,
                                                                tree.list         = seq.created.estipw.glm.propscinsplt.nois.cv1$tree.list,
                                                                lambda.used       = qchisq(0.95, 1),
                                                                val.sample        = data.validation.cont.cont,
                                                                type.var          = "cont",
                                                                propsc.mod.out    = F,
                                                                propsc.mthd       = "GLM",
                                                                propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                val.w             = NULL,
                                                                propsc.mod.insplt = T,
                                                                min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.nois.cv1[[1]],
                                                                 test.data    = data.cont.cont$test.data,
                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                 noise.var    = data.cont.cont$noise.var,
                                                                 corr.split   = data.cont.cont$corr.split,
                                                                 where.split  = data.cont.cont$where.split,
                                                                 dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinsplt.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinsplt.nois.cv1$corr.frst.splt <- as.character(seq.created.estipw.glm.propscinsplt.nois.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("2")

#####################################################################################################################
####################### 3. ipw: GLM Model, inside split, Mis Cov propensity score model, cv1 ########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                               est.used          = "IPW",
                                                               type.var          = "cont",
                                                               propsc.mod.out    = F,
                                                               propsc.mthd       = "GLM",
                                                               propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                               w                 = NULL,
                                                               propsc.mod.insplt = T,
                                                               num.truc.obs      = 30,
                                                               min.node          = 20)

final.tree.estipw.glm.propscinsplt.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont,
                                                               tree.list         = seq.created.estipw.glm.propscinsplt.mis.cv1$tree.list,
                                                               lambda.used       = qchisq(0.95, 1),
                                                               val.sample        = data.validation.cont.cont,
                                                               type.var          = "cont",
                                                               propsc.mod.out    = F,
                                                               propsc.mthd       = "GLM",
                                                               propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                               val.w             = NULL,
                                                               propsc.mod.insplt = T,
                                                               min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.mis.cv1[[1]],
                                                                test.data    = data.cont.cont$test.data,
                                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                                noise.var    = data.cont.cont$noise.var,
                                                                corr.split   = data.cont.cont$corr.split,
                                                                where.split  = data.cont.cont$where.split,
                                                                dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinsplt.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinsplt.mis.cv1$corr.frst.splt <- as.character(seq.created.estipw.glm.propscinsplt.mis.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("3")

#####################################################################################################################
######################## 4. ipw: GLM Model, inside split, True propensity score model, Cv2 ##########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                est.used          = "IPW",
                                                                type.var          = "cont",
                                                                propsc.mod.out    = F,
                                                                propsc.mthd       = "GLM",
                                                                propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                w                 = NULL,
                                                                propsc.mod.insplt = T,
                                                                num.truc.obs      = 30,
                                                                min.node          = 20)

final.tree.estipw.glm.propscinsplt.true.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont,
                                                                tree.list        = seq.created.estipw.glm.propscinsplt.true.cv2$tree.list,
                                                                type.var         = "cont",
                                                                seed             = a[job.number],
                                                                n.cv             = 5,
                                                                propsc.mod.out   = F, 
                                                                propsc.mthd      = "GLM", 
                                                                propsc.form.true = "A ~ X1 + X2 + X3",
                                                                min.obs.mod      = 5)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.true.cv2[[1]],
                                                                 test.data    = data.cont.cont$test.data,
                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                 noise.var    = data.cont.cont$noise.var,
                                                                 corr.split   = data.cont.cont$corr.split,
                                                                 where.split  = data.cont.cont$where.split,
                                                                 dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinsplt.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinsplt.true.cv2$corr.frst.splt <- as.character(seq.created.estipw.glm.propscinsplt.true.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("4")

#####################################################################################################################
####################### 5. ipw: GLM Model, inside split, Mis Func propensity score model, Cv2 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                est.used          = "IPW",
                                                                type.var          = "cont",
                                                                propsc.mod.out    = F,
                                                                propsc.mthd       = "GLM",
                                                                propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                w                 = NULL,
                                                                propsc.mod.insplt = T,
                                                                num.truc.obs      = 30,
                                                                min.node          = 20)

final.tree.estipw.glm.propscinsplt.nois.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont,
                                                                tree.list        = seq.created.estipw.glm.propscinsplt.nois.cv2$tree.list,
                                                                type.var         = "cont",
                                                                seed             = a[job.number],
                                                                n.cv             = 5,
                                                                propsc.mod.out   = F, 
                                                                propsc.mthd      = "GLM", 
                                                                propsc.form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                min.obs.mod      = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.nois.cv2[[1]],
                                                                 test.data    = data.cont.cont$test.data,
                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                 noise.var    = data.cont.cont$noise.var,
                                                                 corr.split   = data.cont.cont$corr.split,
                                                                 where.split  = data.cont.cont$where.split,
                                                                 dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinsplt.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinsplt.nois.cv2$corr.frst.splt <- as.character(seq.created.estipw.glm.propscinsplt.nois.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("5")

#####################################################################################################################
######################## 6. ipw: GLM Model, inside split, Mis Cov propensity score model, Cv2 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                               est.used          = "IPW",
                                                               type.var          = "cont",
                                                               propsc.mod.out    = F,
                                                               propsc.mthd       = "GLM",
                                                               propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                               w                 = NULL,
                                                               propsc.mod.insplt = T,
                                                               num.truc.obs      = 30,
                                                               min.node          = 20)

final.tree.estipw.glm.propscinsplt.mis.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont,
                                                               tree.list        = seq.created.estipw.glm.propscinsplt.mis.cv2$tree.list,
                                                               type.var         = "cont",
                                                               seed             = a[job.number],
                                                               n.cv             = 5,
                                                               propsc.mod.out   = F, 
                                                               propsc.mthd      = "GLM", 
                                                               propsc.form.true = "A ~ X1 + X3 + X4 + X5 + X6",
                                                               min.obs.mod      = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.mis.cv2[[1]],
                                                                test.data    = data.cont.cont$test.data,
                                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                                noise.var    = data.cont.cont$noise.var,
                                                                corr.split   = data.cont.cont$corr.split,
                                                                where.split  = data.cont.cont$where.split,
                                                                dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinsplt.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinsplt.mis.cv2$corr.frst.splt <- as.character(seq.created.estipw.glm.propscinsplt.mis.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("6")

performance.hetero.ipw <- list(glm.propscinsplt.true.cv1 = eval.final.estipw.glm.propscinsplt.true.cv1,
                               glm.propscinsplt.nois.cv1 = eval.final.estipw.glm.propscinsplt.nois.cv1,
                               glm.propscinsplt.mis.cv1  = eval.final.estipw.glm.propscinsplt.mis.cv1,
                               glm.propscinsplt.true.cv2 = eval.final.estipw.glm.propscinsplt.true.cv2,
                               glm.propscinsplt.nois.cv2 = eval.final.estipw.glm.propscinsplt.nois.cv2,
                               glm.propscinsplt.mis.cv2  = eval.final.estipw.glm.propscinsplt.mis.cv2)

#####################################################################################################################
############################ 1. g: GLM Model, inside split, True adjustment model, Cv1 ##############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.true.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                           est.used          = "G",
                                                           type.var          = "cont",
                                                           adj.mod.out       = F,
                                                           adj.mthd          = "GLM",
                                                           adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                           w                 = NULL,
                                                           adj.mod.insplt    = T,
                                                           num.truc.obs      = 10,
                                                           min.node          = 10)

final.tree.estg.glm.modinsplt.true.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                         tree.list         = seq.created.estg.glm.modinsplt.true.cv1$tree.list,
                                                         lambda.used       = qchisq(0.95, 1),
                                                         val.sample        = data.validation.cont.cont,              
                                                         type.var          = "cont",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                         val.w             = NULL,
                                                         adj.mod.insplt    = T,
                                                         min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estg.glm.modinsplt.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinsplt.true.cv1[[1]],
                                                            test.data    = data.cont.cont$test.data,
                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                            noise.var    = data.cont.cont$noise.var,
                                                            corr.split   = data.cont.cont$corr.split,
                                                            where.split  = data.cont.cont$where.split,
                                                            dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modinsplt.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinsplt.true.cv1$corr.frst.splt <- as.character(seq.created.estg.glm.modinsplt.true.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("1")

#####################################################################################################################
############################# 2. g: GLM Model, inside split, Noisy adjustment model, cv1 ############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.nois.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                           est.used          = "G",
                                                           type.var          = "cont",
                                                           adj.mod.out       = F,
                                                           adj.mthd          = "GLM",
                                                           adj.form.true     = NULL,
                                                           w                 = NULL,
                                                           adj.mod.insplt    = T,
                                                           num.truc.obs      = 15,
                                                           min.node          = 15)

final.tree.estg.glm.modinsplt.nois.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                         tree.list         = seq.created.estg.glm.modinsplt.nois.cv1$tree.list,
                                                         lambda.used       = qchisq(0.95, 1),
                                                         val.sample        = data.validation.cont.cont,              
                                                         type.var          = "cont",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = NULL,
                                                         val.w             = NULL,
                                                         adj.mod.insplt    = T,
                                                         min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinsplt.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinsplt.nois.cv1[[1]],
                                                            test.data    = data.cont.cont$test.data,
                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                            noise.var    = data.cont.cont$noise.var,
                                                            corr.split   = data.cont.cont$corr.split,
                                                            where.split  = data.cont.cont$where.split,
                                                            dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modinsplt.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinsplt.nois.cv1$corr.frst.splt <- as.character(seq.created.estg.glm.modinsplt.nois.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("2")

#####################################################################################################################
########################## 3. g: GLM Model, inside split, Misspecified adjustment model, cv1 ########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                          est.used          = "G",
                                                          type.var          = "cont",
                                                          adj.mod.out       = F,
                                                          adj.mthd          = "GLM",
                                                          adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                          w                 = NULL,
                                                          adj.mod.insplt    = T,
                                                          num.truc.obs      = 15,
                                                          min.node          = 15)

final.tree.estg.glm.modinsplt.mis.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                        tree.list         = seq.created.estg.glm.modinsplt.mis.cv1$tree.list,
                                                        lambda.used       = qchisq(0.95, 1),
                                                        val.sample        = data.validation.cont.cont,              
                                                        type.var          = "cont",
                                                        adj.mod.out       = F,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                        val.w             = NULL,
                                                        adj.mod.insplt    = T,
                                                        min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinsplt.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinsplt.mis.cv1[[1]],
                                                           test.data    = data.cont.cont$test.data,
                                                           true.trt.eff = data.cont.cont$true.trt.eff,
                                                           noise.var    = data.cont.cont$noise.var,
                                                           corr.split   = data.cont.cont$corr.split,
                                                           where.split  = data.cont.cont$where.split,
                                                           dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modinsplt.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinsplt.mis.cv1$corr.frst.splt <- as.character(seq.created.estg.glm.modinsplt.mis.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("3")

#####################################################################################################################
############################ 4. g: GLM Model, inside split, True adjustment model, Cv2 ##############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                           est.used          = "G",
                                                           type.var          = "cont",
                                                           adj.mod.out       = F,
                                                           adj.mthd          = "GLM",
                                                           adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                           w                 = NULL,
                                                           adj.mod.insplt    = T,
                                                           num.truc.obs      = 10,
                                                           min.node          = 10)

final.tree.estg.glm.modinsplt.true.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                         tree.list         = seq.created.estg.glm.modinsplt.true.cv2$tree.list,             
                                                         type.var          = "cont",
                                                         seed              = a[job.number],
                                                         n.cv              = 5,
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                         min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estg.glm.modinsplt.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinsplt.true.cv2[[1]],
                                                            test.data    = data.cont.cont$test.data,
                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                            noise.var    = data.cont.cont$noise.var,
                                                            corr.split   = data.cont.cont$corr.split,
                                                            where.split  = data.cont.cont$where.split,
                                                            dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modinsplt.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinsplt.true.cv2$corr.frst.splt <- as.character(seq.created.estg.glm.modinsplt.true.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("4")

#####################################################################################################################
############################# 5. g: GLM Model, inside split, Noisy adjustment model, Cv2 ############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                           est.used          = "G",
                                                           type.var          = "cont",
                                                           adj.mod.out       = F,
                                                           adj.mthd          = "GLM",
                                                           adj.form.true     = NULL,
                                                           w                 = NULL,
                                                           adj.mod.insplt    = T,
                                                           num.truc.obs      = 15,
                                                           min.node          = 15)

final.tree.estg.glm.modinsplt.nois.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                         tree.list         = seq.created.estg.glm.modinsplt.nois.cv2$tree.list,             
                                                         type.var          = "cont",
                                                         seed              = a[job.number],
                                                         n.cv              = 5,
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = NULL,
                                                         min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinsplt.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinsplt.nois.cv2[[1]],
                                                            test.data    = data.cont.cont$test.data,
                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                            noise.var    = data.cont.cont$noise.var,
                                                            corr.split   = data.cont.cont$corr.split,
                                                            where.split  = data.cont.cont$where.split,
                                                            dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modinsplt.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinsplt.nois.cv2$corr.frst.splt <- as.character(seq.created.estg.glm.modinsplt.nois.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("5")

#####################################################################################################################
########################## 6. g: GLM Model, inside split, Misspecified adjustment model, Cv2 ########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                          est.used          = "G",
                                                          type.var          = "cont",
                                                          adj.mod.out       = F,
                                                          adj.mthd          = "GLM",
                                                          adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                          w                 = NULL,
                                                          adj.mod.insplt    = T,
                                                          num.truc.obs      = 15,
                                                          min.node          = 15)

final.tree.estg.glm.modinsplt.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                        tree.list         = seq.created.estg.glm.modinsplt.mis.cv2$tree.list,             
                                                        type.var          = "cont",
                                                        seed              = a[job.number],
                                                        n.cv              = 5,
                                                        adj.mod.out       = F,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                        min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinsplt.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinsplt.mis.cv2[[1]],
                                                           test.data    = data.cont.cont$test.data,
                                                           true.trt.eff = data.cont.cont$true.trt.eff,
                                                           noise.var    = data.cont.cont$noise.var,
                                                           corr.split   = data.cont.cont$corr.split,
                                                           where.split  = data.cont.cont$where.split,
                                                           dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modinsplt.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinsplt.mis.cv2$corr.frst.splt <- as.character(seq.created.estg.glm.modinsplt.mis.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("6")

performance.hetero.g <- list(glm.modinsplt.true.cv1 = eval.final.estg.glm.modinsplt.true.cv1,
                             glm.modinsplt.nois.cv1 = eval.final.estg.glm.modinsplt.nois.cv1,
                             glm.modinsplt.mis.cv1  = eval.final.estg.glm.modinsplt.mis.cv1,
                             glm.modinsplt.true.cv2 = eval.final.estg.glm.modinsplt.true.cv2,
                             glm.modinsplt.nois.cv2 = eval.final.estg.glm.modinsplt.nois.cv2,
                             glm.modinsplt.mis.cv2  = eval.final.estg.glm.modinsplt.mis.cv2)

#####################################################################################################################
#################### 1. dr: True GLM Model in split, True propensity score model in split, Cv1 ######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                        est.used          = "DR",
                                                                        type.var          = "cont",
                                                                        propsc.mod.out    = F,
                                                                        propsc.mthd       = "GLM",
                                                                        propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                        propsc.mod.insplt = T,
                                                                        adj.mod.out       = F, 
                                                                        adj.mthd          = "GLM", 
                                                                        adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                        adj.mod.insplt    = T, 
                                                                        num.truc.obs      = 30,
                                                                        min.node          = 20)

final.tree.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                       tree.list         = seq.created.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1$tree.list, 
                                                                       lambda.used       = qchisq(0.95, 1),
                                                                       val.sample        = data.validation.cont.cont,
                                                                       type.var          = "cont",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                       propsc.mod.insplt = T,
                                                                       adj.mod.out       = F,
                                                                       adj.mthd          = "GLM",
                                                                       adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                       adj.mod.insplt    = T,
                                                                       min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1[[1]],
                                                                         test.data    = data.cont.cont$test.data,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("1")

#####################################################################################################################
###################### 2. dr: True GLM Model in split, Noisy propensity score model in split, Cv1 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                           est.used          = "DR",
                                                                           type.var          = "cont",
                                                                           propsc.mod.out    = F,
                                                                           propsc.mthd       = "GLM",
                                                                           propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                           propsc.mod.insplt = T,
                                                                           adj.mod.out       = F, 
                                                                           adj.mthd          = "GLM", 
                                                                           adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                           adj.mod.insplt    = T, 
                                                                           num.truc.obs      = 30,
                                                                           min.node          = 20)

final.tree.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                          tree.list         = seq.created.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1$tree.list, 
                                                                          lambda.used       = qchisq(0.95, 1),
                                                                          val.sample        = data.validation.cont.cont,
                                                                          type.var          = "cont",
                                                                          propsc.mod.out    = F,
                                                                          propsc.mthd       = "GLM",
                                                                          propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                          propsc.mod.insplt = T,
                                                                          adj.mod.out       = F,
                                                                          adj.mthd          = "GLM",
                                                                          adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                          adj.mod.insplt    = T,
                                                                          min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1[[1]],
                                                                            test.data    = data.cont.cont$test.data,
                                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                                            noise.var    = data.cont.cont$noise.var,
                                                                            corr.split   = data.cont.cont$corr.split,
                                                                            where.split  = data.cont.cont$where.split,
                                                                            dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("2")

#####################################################################################################################
#################### 3. dr: Noisy GLM Model in split, True propensity score model in split, Cv1 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                           est.used          = "DR",
                                                                           type.var          = "cont",
                                                                           propsc.mod.out    = F,
                                                                           propsc.mthd       = "GLM",
                                                                           propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                           propsc.mod.insplt = T,
                                                                           adj.mod.out       = F, 
                                                                           adj.mthd          = "GLM", 
                                                                           adj.form.true     = NULL, 
                                                                           adj.mod.insplt    = T, 
                                                                           num.truc.obs      = 30,
                                                                           min.node          = 20)

final.tree.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                          tree.list         = seq.created.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1$tree.list, 
                                                                          lambda.used       = qchisq(0.95, 1),
                                                                          val.sample        = data.validation.cont.cont,
                                                                          type.var          = "cont",
                                                                          propsc.mod.out    = F,
                                                                          propsc.mthd       = "GLM",
                                                                          propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                          propsc.mod.insplt = T,
                                                                          adj.mod.out       = F,
                                                                          adj.mthd          = "GLM",
                                                                          adj.form.true     = NULL, 
                                                                          adj.mod.insplt    = T,
                                                                          min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1[[1]],
                                                                            test.data    = data.cont.cont$test.data,
                                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                                            noise.var    = data.cont.cont$noise.var,
                                                                            corr.split   = data.cont.cont$corr.split,
                                                                            where.split  = data.cont.cont$where.split,
                                                                            dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("3")

#####################################################################################################################
##################### 4. dr: Noisy GLM Model in split, Noisy propensity score model in split, Cv1 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                              est.used          = "DR",
                                                                              type.var          = "cont",
                                                                              propsc.mod.out    = F,
                                                                              propsc.mthd       = "GLM",
                                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                              propsc.mod.insplt = T,
                                                                              adj.mod.out       = F, 
                                                                              adj.mthd          = "GLM", 
                                                                              adj.form.true     = NULL, 
                                                                              adj.mod.insplt    = T, 
                                                                              num.truc.obs      = 30,
                                                                              min.node          = 20)

final.tree.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                             tree.list         = seq.created.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1$tree.list, 
                                                                             lambda.used       = qchisq(0.95, 1),
                                                                             val.sample        = data.validation.cont.cont,
                                                                             type.var          = "cont",
                                                                             propsc.mod.out    = F,
                                                                             propsc.mthd       = "GLM",
                                                                             propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                             propsc.mod.insplt = T,
                                                                             adj.mod.out       = F,
                                                                             adj.mthd          = "GLM",
                                                                             adj.form.true     = NULL, 
                                                                             adj.mod.insplt    = T,
                                                                             min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1[[1]],
                                                                               test.data    = data.cont.cont$test.data,
                                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                                               noise.var    = data.cont.cont$noise.var,
                                                                               corr.split   = data.cont.cont$corr.split,
                                                                               where.split  = data.cont.cont$where.split,
                                                                               dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("4")

#####################################################################################################################
############ 5. dr: Misspecified GLM Model in split, Misspecified propensity score model in split, Cv1 ##############
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                        est.used          = "DR",
                                                                        type.var          = "cont",
                                                                        propsc.mod.out    = F,
                                                                        propsc.mthd       = "GLM",
                                                                        propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                        propsc.mod.insplt = T,
                                                                        adj.mod.out       = F, 
                                                                        adj.mthd          = "GLM", 
                                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                        adj.mod.insplt    = T, 
                                                                        num.truc.obs      = 30,
                                                                        min.node          = 20)

final.tree.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                       tree.list         = seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1$tree.list, 
                                                                       lambda.used       = qchisq(0.95, 1),
                                                                       val.sample        = data.validation.cont.cont,
                                                                       type.var          = "cont",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                       propsc.mod.insplt = T,
                                                                       adj.mod.out       = F,
                                                                       adj.mthd          = "GLM",
                                                                       adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                       adj.mod.insplt    = T,
                                                                       min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1[[1]],
                                                                         test.data    = data.cont.cont$test.data,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("5")

#####################################################################################################################
###################### 6. dr: True GLM Model in split, True propensity score model in split, Cv2 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                        est.used          = "DR",
                                                                        type.var          = "cont",
                                                                        propsc.mod.out    = F,
                                                                        propsc.mthd       = "GLM",
                                                                        propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                        propsc.mod.insplt = T,
                                                                        adj.mod.out       = F, 
                                                                        adj.mthd          = "GLM", 
                                                                        adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                        adj.mod.insplt    = T, 
                                                                        num.truc.obs      = 30,
                                                                        min.node          = 20)

final.tree.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                       tree.list         = seq.created.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2$tree.list,
                                                                       type.var          = "cont",
                                                                       seed              = a[job.number], 
                                                                       n.cv              = 5,
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                       min.obs.propsc    = 5,
                                                                       adj.mod.out       = F,
                                                                       adj.mthd          = "GLM",
                                                                       adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                       min.obs.adj       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2[[1]],
                                                                         test.data    = data.cont.cont$test.data,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("6")

#####################################################################################################################
###################### 7. dr: True GLM Model in split, Noisy propensity score model in split, Cv2 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                           est.used          = "DR",
                                                                           type.var          = "cont",
                                                                           propsc.mod.out    = F,
                                                                           propsc.mthd       = "GLM",
                                                                           propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                           propsc.mod.insplt = T,
                                                                           adj.mod.out       = F, 
                                                                           adj.mthd          = "GLM", 
                                                                           adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                           adj.mod.insplt    = T, 
                                                                           num.truc.obs      = 30,
                                                                           min.node          = 20)

final.tree.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                          tree.list         = seq.created.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2$tree.list,
                                                                          type.var          = "cont",
                                                                          seed              = a[job.number], 
                                                                          n.cv              = 5,
                                                                          propsc.mod.out    = F,
                                                                          propsc.mthd       = "GLM",
                                                                          propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                          min.obs.propsc    = 10,
                                                                          adj.mod.out       = F,
                                                                          adj.mthd          = "GLM",
                                                                          adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                          min.obs.adj       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2[[1]],
                                                                            test.data    = data.cont.cont$test.data,
                                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                                            noise.var    = data.cont.cont$noise.var,
                                                                            corr.split   = data.cont.cont$corr.split,
                                                                            where.split  = data.cont.cont$where.split,
                                                                            dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("7")

#####################################################################################################################
##################### 8. dr: Noisy GLM Model in split, True propensity score model in split, Cv2 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                           est.used          = "DR",
                                                                           type.var          = "cont",
                                                                           propsc.mod.out    = F,
                                                                           propsc.mthd       = "GLM",
                                                                           propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                           propsc.mod.insplt = T,
                                                                           adj.mod.out       = F, 
                                                                           adj.mthd          = "GLM", 
                                                                           adj.form.true     = NULL, 
                                                                           adj.mod.insplt    = T, 
                                                                           num.truc.obs      = 30,
                                                                           min.node          = 20)

final.tree.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                          tree.list         = seq.created.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2$tree.list,
                                                                          type.var          = "cont",
                                                                          seed              = a[job.number], 
                                                                          n.cv              = 5,
                                                                          propsc.mod.out    = F,
                                                                          propsc.mthd       = "GLM",
                                                                          propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                          min.obs.propsc    = 5,
                                                                          adj.mod.out       = F,
                                                                          adj.mthd          = "GLM",
                                                                          adj.form.true     = NULL, 
                                                                          min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2[[1]],
                                                                            test.data    = data.cont.cont$test.data,
                                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                                            noise.var    = data.cont.cont$noise.var,
                                                                            corr.split   = data.cont.cont$corr.split,
                                                                            where.split  = data.cont.cont$where.split,
                                                                            dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("8")

#####################################################################################################################
#################### 9. dr: Noisy GLM Model in split, Noisy propensity score model in split, Cv2 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                              est.used          = "DR",
                                                                              type.var          = "cont",
                                                                              propsc.mod.out    = F,
                                                                              propsc.mthd       = "GLM",
                                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                              propsc.mod.insplt = T,
                                                                              adj.mod.out       = F, 
                                                                              adj.mthd          = "GLM", 
                                                                              adj.form.true     = NULL, 
                                                                              adj.mod.insplt    = T, 
                                                                              num.truc.obs      = 30,
                                                                              min.node          = 20)

final.tree.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                             tree.list         = seq.created.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2$tree.list,
                                                                             type.var          = "cont",
                                                                             seed              = a[job.number], 
                                                                             n.cv              = 5,
                                                                             propsc.mod.out    = F,
                                                                             propsc.mthd       = "GLM",
                                                                             propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                             min.obs.propsc    = 10,
                                                                             adj.mod.out       = F,
                                                                             adj.mthd          = "GLM",
                                                                             adj.form.true     = NULL, 
                                                                             min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2[[1]],
                                                                               test.data    = data.cont.cont$test.data,
                                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                                               noise.var    = data.cont.cont$noise.var,
                                                                               corr.split   = data.cont.cont$corr.split,
                                                                               where.split  = data.cont.cont$where.split,
                                                                               dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("9")

#####################################################################################################################
############# 10. dr: Misspecified GLM Model in split, Misspecified propensity score model in split, Cv2 ##############
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                        est.used          = "DR",
                                                                        type.var          = "cont",
                                                                        propsc.mod.out    = F,
                                                                        propsc.mthd       = "GLM",
                                                                        propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                        propsc.mod.insplt = T,
                                                                        adj.mod.out       = F, 
                                                                        adj.mthd          = "GLM", 
                                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                        adj.mod.insplt    = T, 
                                                                        num.truc.obs      = 30,
                                                                        min.node          = 20)

final.tree.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                       tree.list         = seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2$tree.list,
                                                                       type.var          = "cont",
                                                                       seed              = a[job.number], 
                                                                       n.cv              = 5,
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                       min.obs.propsc    = 10,
                                                                       adj.mod.out       = F,
                                                                       adj.mthd          = "GLM",
                                                                       adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                       min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2[[1]],
                                                                         test.data    = data.cont.cont$test.data,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("10")

performance.hetero.drInsplt <- list(adjTGlmInsplt.propscTGlmInsplt.cv1       = eval.final.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1,
                                    adjTGlmInsplt.propscNoisGlmInsplt.cv1    = eval.final.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1,
                                    adjNoisGlmInsplt.propscTGlmInsplt.cv1    = eval.final.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1,
                                    adjNoisGlmInsplt.propscNoisGlmInsplt.cv1 = eval.final.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1,
                                    adjFGlmInsplt.propscFGlmInsplt.cv1       = eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1,
                                    adjTGlmInsplt.propscTGlmInsplt.cv2       = eval.final.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2,
                                    adjTGlmInsplt.propscNoisGlmInsplt.cv2    = eval.final.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2,
                                    adjNoisGlmInsplt.propscTGlmInsplt.cv2    = eval.final.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2,
                                    adjNoisGlmInsplt.propscNoisGlmInsplt.cv2 = eval.final.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2,
                                    adjFGlmInsplt.propscFGlmInsplt.cv2       = eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2)



#####################################################################################################################
################################################### Homogeneous #####################################################
#####################################################################################################################
data.cont.cont            <- makeData.cont.noeff.cont(1000, 1000, coeff.prop.sc = 0.6)
data.used.full.cont.cont  <- data.cont.cont$data.used
data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
# val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ]  

#####################################################################################################################
######################### 7. ipw: GLM Model, inside split, True propensity score model, cv1 #########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.true.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                est.used          = "IPW",
                                                                type.var          = "cont",
                                                                propsc.mod.out    = F,
                                                                propsc.mthd       = "GLM",
                                                                propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                w                 = NULL,
                                                                propsc.mod.insplt = T,
                                                                num.truc.obs      = 30,
                                                                min.node          = 20)

final.tree.estipw.glm.propscinsplt.true.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont,
                                                                tree.list         = seq.created.estipw.glm.propscinsplt.true.cv1$tree.list,
                                                                lambda.used       = qchisq(0.95, 1),
                                                                val.sample        = data.validation.cont.cont,
                                                                type.var          = "cont",
                                                                propsc.mod.out    = F,
                                                                propsc.mthd       = "GLM",
                                                                propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                val.w             = NULL,
                                                                propsc.mod.insplt = T,
                                                                min.obs.mod       = 5)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.true.cv1[[1]],
                                                                 test.data    = data.cont.cont$test.data,
                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                 noise.var    = data.cont.cont$noise.var,
                                                                 corr.split   = data.cont.cont$corr.split,
                                                                 where.split  = data.cont.cont$where.split,
                                                                 dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinsplt.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("7")

#####################################################################################################################
######################## 8. ipw: GLM Model, inside split, Mis Func propensity score model, cv1 ######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.nois.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                est.used          = "IPW",
                                                                type.var          = "cont",
                                                                propsc.mod.out    = F,
                                                                propsc.mthd       = "GLM",
                                                                propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                w                 = NULL,
                                                                propsc.mod.insplt = T,
                                                                num.truc.obs      = 30,
                                                                min.node          = 20)

final.tree.estipw.glm.propscinsplt.nois.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont,
                                                                tree.list         = seq.created.estipw.glm.propscinsplt.nois.cv1$tree.list,
                                                                lambda.used       = qchisq(0.95, 1),
                                                                val.sample        = data.validation.cont.cont,
                                                                type.var          = "cont",
                                                                propsc.mod.out    = F,
                                                                propsc.mthd       = "GLM",
                                                                propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                val.w             = NULL,
                                                                propsc.mod.insplt = T,
                                                                min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.nois.cv1[[1]],
                                                                 test.data    = data.cont.cont$test.data,
                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                 noise.var    = data.cont.cont$noise.var,
                                                                 corr.split   = data.cont.cont$corr.split,
                                                                 where.split  = data.cont.cont$where.split,
                                                                 dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinsplt.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("8")

#####################################################################################################################
######################### 9. ipw: GLM Model, inside split, Mis Cov propensity score model, cv1 ######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                               est.used          = "IPW",
                                                               type.var          = "cont",
                                                               propsc.mod.out    = F,
                                                               propsc.mthd       = "GLM",
                                                               propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                               w                 = NULL,
                                                               propsc.mod.insplt = T,
                                                               num.truc.obs      = 30,
                                                               min.node          = 20)

final.tree.estipw.glm.propscinsplt.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont,
                                                               tree.list         = seq.created.estipw.glm.propscinsplt.mis.cv1$tree.list,
                                                               lambda.used       = qchisq(0.95, 1),
                                                               val.sample        = data.validation.cont.cont,
                                                               type.var          = "cont",
                                                               propsc.mod.out    = F,
                                                               propsc.mthd       = "GLM",
                                                               propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                               val.w             = NULL,
                                                               propsc.mod.insplt = T,
                                                               min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.mis.cv1[[1]],
                                                                test.data    = data.cont.cont$test.data,
                                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                                noise.var    = data.cont.cont$noise.var,
                                                                corr.split   = data.cont.cont$corr.split,
                                                                where.split  = data.cont.cont$where.split,
                                                                dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinsplt.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("9")

#####################################################################################################################
######################## 10. ipw: GLM Model, inside split, True propensity score model, Cv2 #########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                est.used          = "IPW",
                                                                type.var          = "cont",
                                                                propsc.mod.out    = F,
                                                                propsc.mthd       = "GLM",
                                                                propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                w                 = NULL,
                                                                propsc.mod.insplt = T,
                                                                num.truc.obs      = 30,
                                                                min.node          = 20)

final.tree.estipw.glm.propscinsplt.true.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont,
                                                                tree.list        = seq.created.estipw.glm.propscinsplt.true.cv2$tree.list,
                                                                type.var         = "cont",
                                                                seed             = a[job.number],
                                                                n.cv             = 5,
                                                                propsc.mod.out   = F, 
                                                                propsc.mthd      = "GLM", 
                                                                propsc.form.true = "A ~ X1 + X2 + X3",
                                                                min.obs.mod      = 5)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.true.cv2[[1]],
                                                                 test.data    = data.cont.cont$test.data,
                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                 noise.var    = data.cont.cont$noise.var,
                                                                 corr.split   = data.cont.cont$corr.split,
                                                                 where.split  = data.cont.cont$where.split,
                                                                 dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinsplt.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("10")

#####################################################################################################################
####################### 11. ipw: GLM Model, inside node, Mis Func propensity score model, Cv2 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                est.used          = "IPW",
                                                                type.var          = "cont",
                                                                propsc.mod.out    = F,
                                                                propsc.mthd       = "GLM",
                                                                propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                w                 = NULL,
                                                                propsc.mod.insplt = T,
                                                                num.truc.obs      = 30,
                                                                min.node          = 20)

final.tree.estipw.glm.propscinsplt.nois.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont,
                                                                tree.list        = seq.created.estipw.glm.propscinsplt.nois.cv2$tree.list,
                                                                type.var         = "cont",
                                                                seed             = a[job.number],
                                                                n.cv             = 5,
                                                                propsc.mod.out   = F, 
                                                                propsc.mthd      = "GLM", 
                                                                propsc.form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                min.obs.mod      = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.nois.cv2[[1]],
                                                                 test.data    = data.cont.cont$test.data,
                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                 noise.var    = data.cont.cont$noise.var,
                                                                 corr.split   = data.cont.cont$corr.split,
                                                                 where.split  = data.cont.cont$where.split,
                                                                 dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinsplt.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("11")

#####################################################################################################################
######################## 12. ipw: GLM Model, inside node, Mis Cov propensity score model, Cv2 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                               est.used          = "IPW",
                                                               type.var          = "cont",
                                                               propsc.mod.out    = F,
                                                               propsc.mthd       = "GLM",
                                                               propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                               w                 = NULL,
                                                               propsc.mod.insplt = T,
                                                               num.truc.obs      = 30,
                                                               min.node          = 20)

final.tree.estipw.glm.propscinsplt.mis.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont,
                                                               tree.list        = seq.created.estipw.glm.propscinsplt.mis.cv2$tree.list,
                                                               type.var         = "cont",
                                                               seed             = a[job.number],
                                                               n.cv             = 5,
                                                               propsc.mod.out   = F, 
                                                               propsc.mthd      = "GLM", 
                                                               propsc.form.true = "A ~ X1 + X3 + X4 + X5 + X6",
                                                               min.obs.mod      = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.mis.cv2[[1]],
                                                                test.data    = data.cont.cont$test.data,
                                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                                noise.var    = data.cont.cont$noise.var,
                                                                corr.split   = data.cont.cont$corr.split,
                                                                where.split  = data.cont.cont$where.split,
                                                                dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinsplt.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("12")

performance.homo.ipw <- list(glm.propscinsplt.true.cv1 = eval.final.estipw.glm.propscinsplt.true.cv1,
                             glm.propscinsplt.nois.cv1 = eval.final.estipw.glm.propscinsplt.nois.cv1,
                             glm.propscinsplt.mis.cv1  = eval.final.estipw.glm.propscinsplt.mis.cv1,
                             glm.propscinsplt.true.cv2 = eval.final.estipw.glm.propscinsplt.true.cv2,
                             glm.propscinsplt.nois.cv2 = eval.final.estipw.glm.propscinsplt.nois.cv2,
                             glm.propscinsplt.mis.cv2  = eval.final.estipw.glm.propscinsplt.mis.cv2)

#####################################################################################################################
############################## 7. g: GLM Model, inside split, True adjustment model, cv1 #############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.true.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                           est.used          = "G",
                                                           type.var          = "cont",
                                                           adj.mod.out       = F,
                                                           adj.mthd          = "GLM",
                                                           adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                           w                 = NULL,
                                                           adj.mod.insplt    = T,
                                                           num.truc.obs      = 10,
                                                           min.node          = 10)

final.tree.estg.glm.modinsplt.true.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                         tree.list         = seq.created.estg.glm.modinsplt.true.cv1$tree.list,
                                                         lambda.used       = qchisq(0.95, 1),
                                                         val.sample        = data.validation.cont.cont,              
                                                         type.var          = "cont",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                         val.w             = NULL,
                                                         adj.mod.insplt    = T,
                                                         min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estg.glm.modinsplt.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinsplt.true.cv1[[1]],
                                                            test.data    = data.cont.cont$test.data,
                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                            noise.var    = data.cont.cont$noise.var,
                                                            corr.split   = data.cont.cont$corr.split,
                                                            where.split  = data.cont.cont$where.split,
                                                            dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modinsplt.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("7")

#####################################################################################################################
############################## 8. g: GLM Model, inside split, Noisy adjustment model, cv1 ############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.nois.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                           est.used          = "G",
                                                           type.var          = "cont",
                                                           adj.mod.out       = F,
                                                           adj.mthd          = "GLM",
                                                           adj.form.true     = NULL,
                                                           w                 = NULL,
                                                           adj.mod.insplt    = T,
                                                           num.truc.obs      = 15,
                                                           min.node          = 15)

final.tree.estg.glm.modinsplt.nois.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                         tree.list         = seq.created.estg.glm.modinsplt.nois.cv1$tree.list,
                                                         lambda.used       = qchisq(0.95, 1),
                                                         val.sample        = data.validation.cont.cont,              
                                                         type.var          = "cont",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = NULL,
                                                         val.w             = NULL,
                                                         adj.mod.insplt    = T,
                                                         min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinsplt.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinsplt.nois.cv1[[1]],
                                                            test.data    = data.cont.cont$test.data,
                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                            noise.var    = data.cont.cont$noise.var,
                                                            corr.split   = data.cont.cont$corr.split,
                                                            where.split  = data.cont.cont$where.split,
                                                            dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modinsplt.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("8")

#####################################################################################################################
########################## 9. g: GLM Model, inside split, Misspecified adjustment model, cv1 #########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                          est.used          = "G",
                                                          type.var          = "cont",
                                                          adj.mod.out       = F,
                                                          adj.mthd          = "GLM",
                                                          adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                          w                 = NULL,
                                                          adj.mod.insplt    = T,
                                                          num.truc.obs      = 15,
                                                          min.node          = 15)

final.tree.estg.glm.modinsplt.mis.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                        tree.list         = seq.created.estg.glm.modinsplt.mis.cv1$tree.list,
                                                        lambda.used       = qchisq(0.95, 1),
                                                        val.sample        = data.validation.cont.cont,              
                                                        type.var          = "cont",
                                                        adj.mod.out       = F,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                        val.w             = NULL,
                                                        adj.mod.insplt    = T,
                                                        min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinsplt.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinsplt.mis.cv1[[1]],
                                                           test.data    = data.cont.cont$test.data,
                                                           true.trt.eff = data.cont.cont$true.trt.eff,
                                                           noise.var    = data.cont.cont$noise.var,
                                                           corr.split   = data.cont.cont$corr.split,
                                                           where.split  = data.cont.cont$where.split,
                                                           dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modinsplt.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("9")

#####################################################################################################################
############################# 10. g: GLM Model, inside split, True adjustment model, Cv2 #############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                           est.used          = "G",
                                                           type.var          = "cont",
                                                           adj.mod.out       = F,
                                                           adj.mthd          = "GLM",
                                                           adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                           w                 = NULL,
                                                           adj.mod.insplt    = T,
                                                           num.truc.obs      = 10,
                                                           min.node          = 10)

final.tree.estg.glm.modinsplt.true.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                         tree.list         = seq.created.estg.glm.modinsplt.true.cv2$tree.list,             
                                                         type.var          = "cont",
                                                         seed              = a[job.number],
                                                         n.cv              = 5,
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                         min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estg.glm.modinsplt.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinsplt.true.cv2[[1]],
                                                            test.data    = data.cont.cont$test.data,
                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                            noise.var    = data.cont.cont$noise.var,
                                                            corr.split   = data.cont.cont$corr.split,
                                                            where.split  = data.cont.cont$where.split,
                                                            dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modinsplt.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("10")

#####################################################################################################################
############################# 11. g: GLM Model, inside split, Noisy adjustment model, Cv2 ############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                           est.used          = "G",
                                                           type.var          = "cont",
                                                           adj.mod.out       = F,
                                                           adj.mthd          = "GLM",
                                                           adj.form.true     = NULL,
                                                           w                 = NULL,
                                                           adj.mod.insplt    = T,
                                                           num.truc.obs      = 15,
                                                           min.node          = 15)

final.tree.estg.glm.modinsplt.nois.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                         tree.list         = seq.created.estg.glm.modinsplt.nois.cv2$tree.list,             
                                                         type.var          = "cont",
                                                         seed              = a[job.number],
                                                         n.cv              = 5,
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = NULL,
                                                         min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinsplt.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinsplt.nois.cv2[[1]],
                                                            test.data    = data.cont.cont$test.data,
                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                            noise.var    = data.cont.cont$noise.var,
                                                            corr.split   = data.cont.cont$corr.split,
                                                            where.split  = data.cont.cont$where.split,
                                                            dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modinsplt.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("11")

#####################################################################################################################
########################## 12. g: GLM Model, inside split, Misspecified adjustment model, Cv2 ########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                          est.used          = "G",
                                                          type.var          = "cont",
                                                          adj.mod.out       = F,
                                                          adj.mthd          = "GLM",
                                                          adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                          w                 = NULL,
                                                          adj.mod.insplt    = T,
                                                          num.truc.obs      = 15,
                                                          min.node          = 15)

final.tree.estg.glm.modinsplt.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                        tree.list         = seq.created.estg.glm.modinsplt.mis.cv2$tree.list,             
                                                        type.var          = "cont",
                                                        seed              = a[job.number],
                                                        n.cv              = 5,
                                                        adj.mod.out       = F,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                        min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinsplt.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinsplt.mis.cv2[[1]],
                                                           test.data    = data.cont.cont$test.data,
                                                           true.trt.eff = data.cont.cont$true.trt.eff,
                                                           noise.var    = data.cont.cont$noise.var,
                                                           corr.split   = data.cont.cont$corr.split,
                                                           where.split  = data.cont.cont$where.split,
                                                           dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modinsplt.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("12")

performance.homo.g <- list(glm.modinsplt.true.cv1 = eval.final.estg.glm.modinsplt.true.cv1,
                           glm.modinsplt.nois.cv1 = eval.final.estg.glm.modinsplt.nois.cv1,
                           glm.modinsplt.mis.cv1  = eval.final.estg.glm.modinsplt.mis.cv1,
                           glm.modinsplt.true.cv2 = eval.final.estg.glm.modinsplt.true.cv2,
                           glm.modinsplt.nois.cv2 = eval.final.estg.glm.modinsplt.nois.cv2,
                           glm.modinsplt.mis.cv2  = eval.final.estg.glm.modinsplt.mis.cv2)

#####################################################################################################################
#################### 11. dr: True GLM Model in split, True propensity score model in split, Cv1 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                        est.used          = "DR",
                                                                        type.var          = "cont",
                                                                        propsc.mod.out    = F,
                                                                        propsc.mthd       = "GLM",
                                                                        propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                        propsc.mod.insplt = T,
                                                                        adj.mod.out       = F, 
                                                                        adj.mthd          = "GLM", 
                                                                        adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                        adj.mod.insplt    = T, 
                                                                        num.truc.obs      = 30,
                                                                        min.node          = 20)

final.tree.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                       tree.list         = seq.created.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1$tree.list, 
                                                                       lambda.used       = qchisq(0.95, 1),
                                                                       val.sample        = data.validation.cont.cont,
                                                                       type.var          = "cont",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                       propsc.mod.insplt = T,
                                                                       adj.mod.out       = F,
                                                                       adj.mthd          = "GLM",
                                                                       adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                       adj.mod.insplt    = T,
                                                                       min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1[[1]],
                                                                         test.data    = data.cont.cont$test.data,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("11")

#####################################################################################################################
##################### 12. dr: True GLM Model in split, Noisy propensity score model in split, Cv1 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                           est.used          = "DR",
                                                                           type.var          = "cont",
                                                                           propsc.mod.out    = F,
                                                                           propsc.mthd       = "GLM",
                                                                           propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                           propsc.mod.insplt = T,
                                                                           adj.mod.out       = F, 
                                                                           adj.mthd          = "GLM", 
                                                                           adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                           adj.mod.insplt    = T, 
                                                                           num.truc.obs      = 30,
                                                                           min.node          = 20)

final.tree.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                          tree.list         = seq.created.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1$tree.list, 
                                                                          lambda.used       = qchisq(0.95, 1),
                                                                          val.sample        = data.validation.cont.cont,
                                                                          type.var          = "cont",
                                                                          propsc.mod.out    = F,
                                                                          propsc.mthd       = "GLM",
                                                                          propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                          propsc.mod.insplt = T,
                                                                          adj.mod.out       = F,
                                                                          adj.mthd          = "GLM",
                                                                          adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                          adj.mod.insplt    = T,
                                                                          min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1[[1]],
                                                                            test.data    = data.cont.cont$test.data,
                                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                                            noise.var    = data.cont.cont$noise.var,
                                                                            corr.split   = data.cont.cont$corr.split,
                                                                            where.split  = data.cont.cont$where.split,
                                                                            dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("12")

#####################################################################################################################
#################### 13. dr: Noisy GLM Model in split, True propensity score model in split, Cv1 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                           est.used          = "DR",
                                                                           type.var          = "cont",
                                                                           propsc.mod.out    = F,
                                                                           propsc.mthd       = "GLM",
                                                                           propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                           propsc.mod.insplt = T,
                                                                           adj.mod.out       = F, 
                                                                           adj.mthd          = "GLM", 
                                                                           adj.form.true     = NULL, 
                                                                           adj.mod.insplt    = T, 
                                                                           num.truc.obs      = 30,
                                                                           min.node          = 20)

final.tree.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                          tree.list         = seq.created.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1$tree.list, 
                                                                          lambda.used       = qchisq(0.95, 1),
                                                                          val.sample        = data.validation.cont.cont,
                                                                          type.var          = "cont",
                                                                          propsc.mod.out    = F,
                                                                          propsc.mthd       = "GLM",
                                                                          propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                          propsc.mod.insplt = T,
                                                                          adj.mod.out       = F,
                                                                          adj.mthd          = "GLM",
                                                                          adj.form.true     = NULL, 
                                                                          adj.mod.insplt    = T,
                                                                          min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1[[1]],
                                                                            test.data    = data.cont.cont$test.data,
                                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                                            noise.var    = data.cont.cont$noise.var,
                                                                            corr.split   = data.cont.cont$corr.split,
                                                                            where.split  = data.cont.cont$where.split,
                                                                            dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("13")

#####################################################################################################################
#################### 14. dr: Noisy GLM Model in split, Noisy propensity score model in split, Cv1 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                              est.used          = "DR",
                                                                              type.var          = "cont",
                                                                              propsc.mod.out    = F,
                                                                              propsc.mthd       = "GLM",
                                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                              propsc.mod.insplt = T,
                                                                              adj.mod.out       = F, 
                                                                              adj.mthd          = "GLM", 
                                                                              adj.form.true     = NULL, 
                                                                              adj.mod.insplt    = T, 
                                                                              num.truc.obs      = 30,
                                                                              min.node          = 20)

final.tree.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                             tree.list         = seq.created.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1$tree.list, 
                                                                             lambda.used       = qchisq(0.95, 1),
                                                                             val.sample        = data.validation.cont.cont,
                                                                             type.var          = "cont",
                                                                             propsc.mod.out    = F,
                                                                             propsc.mthd       = "GLM",
                                                                             propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                             propsc.mod.insplt = T,
                                                                             adj.mod.out       = F,
                                                                             adj.mthd          = "GLM",
                                                                             adj.form.true     = NULL, 
                                                                             adj.mod.insplt    = T,
                                                                             min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1[[1]],
                                                                               test.data    = data.cont.cont$test.data,
                                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                                               noise.var    = data.cont.cont$noise.var,
                                                                               corr.split   = data.cont.cont$corr.split,
                                                                               where.split  = data.cont.cont$where.split,
                                                                               dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("14")

#####################################################################################################################
############ 15. dr: Misspecified GLM Model in split, Misspecified propensity score model in split, Cv1 #############
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                        est.used          = "DR",
                                                                        type.var          = "cont",
                                                                        propsc.mod.out    = F,
                                                                        propsc.mthd       = "GLM",
                                                                        propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                        propsc.mod.insplt = T,
                                                                        adj.mod.out       = F, 
                                                                        adj.mthd          = "GLM", 
                                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                        adj.mod.insplt    = T, 
                                                                        num.truc.obs      = 30,
                                                                        min.node          = 20)

final.tree.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                       tree.list         = seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1$tree.list, 
                                                                       lambda.used       = qchisq(0.95, 1),
                                                                       val.sample        = data.validation.cont.cont,
                                                                       type.var          = "cont",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                       propsc.mod.insplt = T,
                                                                       adj.mod.out       = F,
                                                                       adj.mthd          = "GLM",
                                                                       adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                       adj.mod.insplt    = T,
                                                                       min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1[[1]],
                                                                         test.data    = data.cont.cont$test.data,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("15")

#####################################################################################################################
##################### 16. dr: True GLM Model in split, True propensity score model in split, Cv2 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                        est.used          = "DR",
                                                                        type.var          = "cont",
                                                                        propsc.mod.out    = F,
                                                                        propsc.mthd       = "GLM",
                                                                        propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                        propsc.mod.insplt = T,
                                                                        adj.mod.out       = F, 
                                                                        adj.mthd          = "GLM", 
                                                                        adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                        adj.mod.insplt    = T, 
                                                                        num.truc.obs      = 30,
                                                                        min.node          = 20)

final.tree.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                       tree.list         = seq.created.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2$tree.list,
                                                                       type.var          = "cont",
                                                                       seed              = a[job.number], 
                                                                       n.cv              = 5,
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                       min.obs.propsc    = 5,
                                                                       adj.mod.out       = F,
                                                                       adj.mthd          = "GLM",
                                                                       adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                       min.obs.adj       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2[[1]],
                                                                         test.data    = data.cont.cont$test.data,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("16")

#####################################################################################################################
#################### 17. dr: True GLM Model in split, Noisy propensity score model in split, Cv2 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                           est.used          = "DR",
                                                                           type.var          = "cont",
                                                                           propsc.mod.out    = F,
                                                                           propsc.mthd       = "GLM",
                                                                           propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                           propsc.mod.insplt = T,
                                                                           adj.mod.out       = F, 
                                                                           adj.mthd          = "GLM", 
                                                                           adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                           adj.mod.insplt    = T, 
                                                                           num.truc.obs      = 30,
                                                                           min.node          = 20)

final.tree.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                          tree.list         = seq.created.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2$tree.list,
                                                                          type.var          = "cont",
                                                                          seed              = a[job.number], 
                                                                          n.cv              = 5,
                                                                          propsc.mod.out    = F,
                                                                          propsc.mthd       = "GLM",
                                                                          propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                          min.obs.propsc    = 10,
                                                                          adj.mod.out       = F,
                                                                          adj.mthd          = "GLM",
                                                                          adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                          min.obs.adj       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2[[1]],
                                                                            test.data    = data.cont.cont$test.data,
                                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                                            noise.var    = data.cont.cont$noise.var,
                                                                            corr.split   = data.cont.cont$corr.split,
                                                                            where.split  = data.cont.cont$where.split,
                                                                            dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("17")

#####################################################################################################################
##################### 18. dr: Noisy GLM Model in split, True propensity score model in split, Cv2 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                           est.used          = "DR",
                                                                           type.var          = "cont",
                                                                           propsc.mod.out    = F,
                                                                           propsc.mthd       = "GLM",
                                                                           propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                           propsc.mod.insplt = T,
                                                                           adj.mod.out       = F, 
                                                                           adj.mthd          = "GLM", 
                                                                           adj.form.true     = NULL, 
                                                                           adj.mod.insplt    = T, 
                                                                           num.truc.obs      = 30,
                                                                           min.node          = 20)

final.tree.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                          tree.list         = seq.created.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2$tree.list,
                                                                          type.var          = "cont",
                                                                          seed              = a[job.number], 
                                                                          n.cv              = 5,
                                                                          propsc.mod.out    = F,
                                                                          propsc.mthd       = "GLM",
                                                                          propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                          min.obs.propsc    = 5,
                                                                          adj.mod.out       = F,
                                                                          adj.mthd          = "GLM",
                                                                          adj.form.true     = NULL, 
                                                                          min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2[[1]],
                                                                            test.data    = data.cont.cont$test.data,
                                                                            true.trt.eff = data.cont.cont$true.trt.eff,
                                                                            noise.var    = data.cont.cont$noise.var,
                                                                            corr.split   = data.cont.cont$corr.split,
                                                                            where.split  = data.cont.cont$where.split,
                                                                            dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("18")

#####################################################################################################################
#################### 19. dr: Noisy GLM Model in split, Noisy propensity score model in split, Cv2 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                              est.used          = "DR",
                                                                              type.var          = "cont",
                                                                              propsc.mod.out    = F,
                                                                              propsc.mthd       = "GLM",
                                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                              propsc.mod.insplt = T,
                                                                              adj.mod.out       = F, 
                                                                              adj.mthd          = "GLM", 
                                                                              adj.form.true     = NULL, 
                                                                              adj.mod.insplt    = T, 
                                                                              num.truc.obs      = 30,
                                                                              min.node          = 20)

final.tree.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                             tree.list         = seq.created.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2$tree.list,
                                                                             type.var          = "cont",
                                                                             seed              = a[job.number], 
                                                                             n.cv              = 5,
                                                                             propsc.mod.out    = F,
                                                                             propsc.mthd       = "GLM",
                                                                             propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                             min.obs.propsc    = 10,
                                                                             adj.mod.out       = F,
                                                                             adj.mthd          = "GLM",
                                                                             adj.form.true     = NULL, 
                                                                             min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2[[1]],
                                                                               test.data    = data.cont.cont$test.data,
                                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                                               noise.var    = data.cont.cont$noise.var,
                                                                               corr.split   = data.cont.cont$corr.split,
                                                                               where.split  = data.cont.cont$where.split,
                                                                               dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("19")

#####################################################################################################################
############# 20. dr: Misspecified GLM Model in split, Misspecified propensity score model in split, Cv2 ############
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                        est.used          = "DR",
                                                                        type.var          = "cont",
                                                                        propsc.mod.out    = F,
                                                                        propsc.mthd       = "GLM",
                                                                        propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                        propsc.mod.insplt = T,
                                                                        adj.mod.out       = F, 
                                                                        adj.mthd          = "GLM", 
                                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                        adj.mod.insplt    = T, 
                                                                        num.truc.obs      = 30,
                                                                        min.node          = 20)

final.tree.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                       tree.list         = seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2$tree.list,
                                                                       type.var          = "cont",
                                                                       seed              = a[job.number], 
                                                                       n.cv              = 5,
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                       min.obs.propsc    = 10,
                                                                       adj.mod.out       = F,
                                                                       adj.mthd          = "GLM",
                                                                       adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                       min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2[[1]],
                                                                         test.data    = data.cont.cont$test.data,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("20")

performance.homo.drInsplt <- list(adjTGlmInsplt.propscTGlmInsplt.cv1       = eval.final.estdr.adjTGlmInsplt.propscTGlmInsplt.cv1,
                                  adjTGlmInsplt.propscNoisGlmInsplt.cv1    = eval.final.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv1,
                                  adjNoisGlmInsplt.propscTGlmInsplt.cv1    = eval.final.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv1,
                                  adjNoisGlmInsplt.propscNoisGlmInsplt.cv1 = eval.final.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv1,
                                  adjFGlmInsplt.propscFGlmInsplt.cv1       = eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1,
                                  adjTGlmInsplt.propscTGlmInsplt.cv2       = eval.final.estdr.adjTGlmInsplt.propscTGlmInsplt.cv2,
                                  adjTGlmInsplt.propscNoisGlmInsplt.cv2    = eval.final.estdr.adjTGlmInsplt.propscNoisGlmInsplt.cv2,
                                  adjNoisGlmInsplt.propscTGlmInsplt.cv2    = eval.final.estdr.adjNoisGlmInsplt.propscTGlmInsplt.cv2,
                                  adjNoisGlmInsplt.propscNoisGlmInsplt.cv2 = eval.final.estdr.adjNoisGlmInsplt.propscNoisGlmInsplt.cv2,
                                  adjFGlmInsplt.propscFGlmInsplt.cv2       = eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2)

