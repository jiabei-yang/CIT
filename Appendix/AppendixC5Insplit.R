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

# Pretend X2 is unmeasured for unmeasured cov
data.used.full.cont.cont.mis <- data.used.full.cont.cont %>%
  select(-X2)
data.used.cont.cont.mis <- data.used.cont.cont %>%
  select(-X2)
data.validation.cont.cont.mis <- data.validation.cont.cont %>%
  select(-X2)
test.data.mis <- data.cont.cont$test.data %>%
  select(-X2)

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
#################### 3. ipw: GLM Model, inside split, Unmeasured Cov propensity score model, cv1 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                               est.used          = "IPW",
                                                               type.var          = "cont",
                                                               propsc.mod.out    = F,
                                                               propsc.mthd       = "GLM",
                                                               propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                               w                 = NULL,
                                                               propsc.mod.insplt = T,
                                                               num.truc.obs      = 30,
                                                               min.node          = 20)

final.tree.estipw.glm.propscinsplt.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                               tree.list         = seq.created.estipw.glm.propscinsplt.mis.cv1$tree.list,
                                                               lambda.used       = qchisq(0.95, 1),
                                                               val.sample        = data.validation.cont.cont.mis,
                                                               type.var          = "cont",
                                                               propsc.mod.out    = F,
                                                               propsc.mthd       = "GLM",
                                                               propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                               val.w             = NULL,
                                                               propsc.mod.insplt = T,
                                                               min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.mis.cv1[[1]],
                                                                test.data    = test.data.mis,
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
#################### 6. ipw: GLM Model, inside split, Unmeasured cov propensity score model, Cv2 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                               est.used          = "IPW",
                                                               type.var          = "cont",
                                                               propsc.mod.out    = F,
                                                               propsc.mthd       = "GLM",
                                                               propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                               w                 = NULL,
                                                               propsc.mod.insplt = T,
                                                               num.truc.obs      = 30,
                                                               min.node          = 20)

final.tree.estipw.glm.propscinsplt.mis.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont.mis,
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
                                                                test.data    = test.data.mis,
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
############################ 7. g: GLM Model, inside split, True adjustment model, Cv1 ##############################
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
print("7")

#####################################################################################################################
########################### 8. g: GLM Model, inside split, Mis func adjustment model, cv1 ###########################
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
print("8")

#####################################################################################################################
######################### 9. g: GLM Model, inside split, Unmeasured cov adjustment model, cv1 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                          est.used          = "G",
                                                          type.var          = "cont",
                                                          adj.mod.out       = F,
                                                          adj.mthd          = "GLM",
                                                          adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                          w                 = NULL,
                                                          adj.mod.insplt    = T,
                                                          num.truc.obs      = 15,
                                                          min.node          = 15)

final.tree.estg.glm.modinsplt.mis.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                        tree.list         = seq.created.estg.glm.modinsplt.mis.cv1$tree.list,
                                                        lambda.used       = qchisq(0.95, 1),
                                                        val.sample        = data.validation.cont.cont.mis,              
                                                        type.var          = "cont",
                                                        adj.mod.out       = F,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                        val.w             = NULL,
                                                        adj.mod.insplt    = T,
                                                        min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinsplt.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinsplt.mis.cv1[[1]],
                                                           test.data    = test.data.mis,
                                                           true.trt.eff = data.cont.cont$true.trt.eff,
                                                           noise.var    = data.cont.cont$noise.var,
                                                           corr.split   = data.cont.cont$corr.split,
                                                           where.split  = data.cont.cont$where.split,
                                                           dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modinsplt.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinsplt.mis.cv1$corr.frst.splt <- as.character(seq.created.estg.glm.modinsplt.mis.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("9")

#####################################################################################################################
########################### 10. g: GLM Model, inside split, True adjustment model, Cv2 ##############################
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
print("10")

#####################################################################################################################
########################### 11. g: GLM Model, inside split, Mis func adjustment model, Cv2 ##########################
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
print("11")

#####################################################################################################################
######################## 12. g: GLM Model, inside split, Unmeasured cov adjustment model, Cv2 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                          est.used          = "G",
                                                          type.var          = "cont",
                                                          adj.mod.out       = F,
                                                          adj.mthd          = "GLM",
                                                          adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                          w                 = NULL,
                                                          adj.mod.insplt    = T,
                                                          num.truc.obs      = 15,
                                                          min.node          = 15)

final.tree.estg.glm.modinsplt.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont.mis,
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
                                                           test.data    = test.data.mis,
                                                           true.trt.eff = data.cont.cont$true.trt.eff,
                                                           noise.var    = data.cont.cont$noise.var,
                                                           corr.split   = data.cont.cont$corr.split,
                                                           where.split  = data.cont.cont$where.split,
                                                           dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modinsplt.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinsplt.mis.cv2$corr.frst.splt <- as.character(seq.created.estg.glm.modinsplt.mis.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("12")

performance.hetero.g <- list(glm.modinsplt.true.cv1 = eval.final.estg.glm.modinsplt.true.cv1,
                             glm.modinsplt.nois.cv1 = eval.final.estg.glm.modinsplt.nois.cv1,
                             glm.modinsplt.mis.cv1  = eval.final.estg.glm.modinsplt.mis.cv1,
                             glm.modinsplt.true.cv2 = eval.final.estg.glm.modinsplt.true.cv2,
                             glm.modinsplt.nois.cv2 = eval.final.estg.glm.modinsplt.nois.cv2,
                             glm.modinsplt.mis.cv2  = eval.final.estg.glm.modinsplt.mis.cv2)

#####################################################################################################################
#################### 13. dr: True GLM Model in split, True propensity score model in split, Cv1 #####################
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
print("13")

#####################################################################################################################
#################### 14. dr: True GLM Model in split, Mis func propensity score model in split, Cv1 #################
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
print("14")

#####################################################################################################################
################## 15. dr: Mis func GLM Model in split, True propensity score model in split, Cv1 ###################
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
print("15")

#####################################################################################################################
################# 16. dr: Mis func GLM Model in split, Mis func propensity score model in split, Cv1 ################
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
print("16")

#####################################################################################################################
######### 17. dr: Unmeasured cov GLM Model in split, Unmeasured cov propensity score model in split, Cv1 ############
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
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

final.tree.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                                       tree.list         = seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1$tree.list, 
                                                                       lambda.used       = qchisq(0.95, 1),
                                                                       val.sample        = data.validation.cont.cont.mis,
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
                                                                         test.data    = test.data.mis,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("17")

#####################################################################################################################
###################### 18. dr: True GLM Model in split, True propensity score model in split, Cv2 ###################
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
print("18")

#####################################################################################################################
#################### 19. dr: True GLM Model in split, Mis func propensity score model in split, Cv2 #################
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
print("19")

#####################################################################################################################
################### 20. dr: Mis func GLM Model in split, True propensity score model in split, Cv2 ##################
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
print("20")

#####################################################################################################################
################ 21. dr: Mis func GLM Model in split, Mis func propensity score model in split, Cv2 #################
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
print("21")

#####################################################################################################################
########## 22. dr: Unmeasured cov GLM Model in split, Unmeasured cov propensity score model in split, Cv2 ###########
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
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

final.tree.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont.mis,
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
                                                                         test.data    = test.data.mis,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("22")

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

# Pretend X2 is unmeasured for unmeasured cov
data.used.full.cont.cont.mis <- data.used.full.cont.cont %>%
  select(-X2)
data.used.cont.cont.mis <- data.used.cont.cont %>%
  select(-X2)
data.validation.cont.cont.mis <- data.validation.cont.cont %>%
  select(-X2)
test.data.mis <- data.cont.cont$test.data %>%
  select(-X2)

#####################################################################################################################
######################### 23. ipw: GLM Model, inside split, True propensity score model, cv1 ########################
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
print("23")

#####################################################################################################################
####################### 24. ipw: GLM Model, inside split, Mis Func propensity score model, cv1 ######################
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
print("24")

#####################################################################################################################
#################### 25. ipw: GLM Model, inside split, Unmeasured cov propensity score model, cv1 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                               est.used          = "IPW",
                                                               type.var          = "cont",
                                                               propsc.mod.out    = F,
                                                               propsc.mthd       = "GLM",
                                                               propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                               w                 = NULL,
                                                               propsc.mod.insplt = T,
                                                               num.truc.obs      = 30,
                                                               min.node          = 20)

final.tree.estipw.glm.propscinsplt.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                               tree.list         = seq.created.estipw.glm.propscinsplt.mis.cv1$tree.list,
                                                               lambda.used       = qchisq(0.95, 1),
                                                               val.sample        = data.validation.cont.cont.mis,
                                                               type.var          = "cont",
                                                               propsc.mod.out    = F,
                                                               propsc.mthd       = "GLM",
                                                               propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                               val.w             = NULL,
                                                               propsc.mod.insplt = T,
                                                               min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.mis.cv1[[1]],
                                                                test.data    = test.data.mis,
                                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                                noise.var    = data.cont.cont$noise.var,
                                                                corr.split   = data.cont.cont$corr.split,
                                                                where.split  = data.cont.cont$where.split,
                                                                dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinsplt.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("25")

#####################################################################################################################
######################## 26. ipw: GLM Model, inside split, True propensity score model, Cv2 #########################
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
print("26")

#####################################################################################################################
####################### 27. ipw: GLM Model, inside node, Mis Func propensity score model, Cv2 #######################
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
print("27")

#####################################################################################################################
##################### 28. ipw: GLM Model, inside node, Unmeasured cov propensity score model, Cv2 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                               est.used          = "IPW",
                                                               type.var          = "cont",
                                                               propsc.mod.out    = F,
                                                               propsc.mthd       = "GLM",
                                                               propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                               w                 = NULL,
                                                               propsc.mod.insplt = T,
                                                               num.truc.obs      = 30,
                                                               min.node          = 20)

final.tree.estipw.glm.propscinsplt.mis.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont.mis,
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
                                                                test.data    = test.data.mis,
                                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                                noise.var    = data.cont.cont$noise.var,
                                                                corr.split   = data.cont.cont$corr.split,
                                                                where.split  = data.cont.cont$where.split,
                                                                dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinsplt.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("28")

performance.homo.ipw <- list(glm.propscinsplt.true.cv1 = eval.final.estipw.glm.propscinsplt.true.cv1,
                             glm.propscinsplt.nois.cv1 = eval.final.estipw.glm.propscinsplt.nois.cv1,
                             glm.propscinsplt.mis.cv1  = eval.final.estipw.glm.propscinsplt.mis.cv1,
                             glm.propscinsplt.true.cv2 = eval.final.estipw.glm.propscinsplt.true.cv2,
                             glm.propscinsplt.nois.cv2 = eval.final.estipw.glm.propscinsplt.nois.cv2,
                             glm.propscinsplt.mis.cv2  = eval.final.estipw.glm.propscinsplt.mis.cv2)

#####################################################################################################################
############################ 29. g: GLM Model, inside split, True adjustment model, cv1 #############################
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
print("29")

#####################################################################################################################
########################### 30. g: GLM Model, inside split, Mis func adjustment model, cv1 ##########################
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
print("30")

#####################################################################################################################
######################## 31. g: GLM Model, inside split, Unmeasured cov adjustment model, cv1 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                          est.used          = "G",
                                                          type.var          = "cont",
                                                          adj.mod.out       = F,
                                                          adj.mthd          = "GLM",
                                                          adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                          w                 = NULL,
                                                          adj.mod.insplt    = T,
                                                          num.truc.obs      = 15,
                                                          min.node          = 15)

final.tree.estg.glm.modinsplt.mis.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                        tree.list         = seq.created.estg.glm.modinsplt.mis.cv1$tree.list,
                                                        lambda.used       = qchisq(0.95, 1),
                                                        val.sample        = data.validation.cont.cont.mis,              
                                                        type.var          = "cont",
                                                        adj.mod.out       = F,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                        val.w             = NULL,
                                                        adj.mod.insplt    = T,
                                                        min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinsplt.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinsplt.mis.cv1[[1]],
                                                           test.data    = test.data.mis,
                                                           true.trt.eff = data.cont.cont$true.trt.eff,
                                                           noise.var    = data.cont.cont$noise.var,
                                                           corr.split   = data.cont.cont$corr.split,
                                                           where.split  = data.cont.cont$where.split,
                                                           dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modinsplt.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("31")

#####################################################################################################################
############################# 32. g: GLM Model, inside split, True adjustment model, Cv2 ############################
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
print("32")

#####################################################################################################################
############################ 33. g: GLM Model, inside split, Mis func adjustment model, Cv2 #########################
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
print("33")

#####################################################################################################################
######################## 34. g: GLM Model, inside split, Unmeasured cov adjustment model, Cv2 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinsplt.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                          est.used          = "G",
                                                          type.var          = "cont",
                                                          adj.mod.out       = F,
                                                          adj.mthd          = "GLM",
                                                          adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                          w                 = NULL,
                                                          adj.mod.insplt    = T,
                                                          num.truc.obs      = 15,
                                                          min.node          = 15)

final.tree.estg.glm.modinsplt.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont.mis,
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
                                                           test.data    = test.data.mis,
                                                           true.trt.eff = data.cont.cont$true.trt.eff,
                                                           noise.var    = data.cont.cont$noise.var,
                                                           corr.split   = data.cont.cont$corr.split,
                                                           where.split  = data.cont.cont$where.split,
                                                           dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modinsplt.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("34")

performance.homo.g <- list(glm.modinsplt.true.cv1 = eval.final.estg.glm.modinsplt.true.cv1,
                           glm.modinsplt.nois.cv1 = eval.final.estg.glm.modinsplt.nois.cv1,
                           glm.modinsplt.mis.cv1  = eval.final.estg.glm.modinsplt.mis.cv1,
                           glm.modinsplt.true.cv2 = eval.final.estg.glm.modinsplt.true.cv2,
                           glm.modinsplt.nois.cv2 = eval.final.estg.glm.modinsplt.nois.cv2,
                           glm.modinsplt.mis.cv2  = eval.final.estg.glm.modinsplt.mis.cv2)

#####################################################################################################################
#################### 35. dr: True GLM Model in split, True propensity score model in split, Cv1 #####################
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
print("dr-adjT-propscT-cv1")

#####################################################################################################################
#################### 36. dr: True GLM Model in split, Mis func propensity score model in split, Cv1 #################
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
print("dr-adjT-propscNois-cv1")

#####################################################################################################################
################## 37. dr: Mis func GLM Model in split, True propensity score model in split, Cv1 ###################
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
print("dr-adjNois-propscT-cv1")

#####################################################################################################################
################# 38. dr: Mis func GLM Model in split, Mis func propensity score model in split, Cv1 ################
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
print("dr-adjNois-propscNois-cv1")

#####################################################################################################################
########## 39. dr: Unmeasured cov GLM Model in split, Unmeasured cov propensity score model in split, Cv1 ###########
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
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

final.tree.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                                       tree.list         = seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1$tree.list, 
                                                                       lambda.used       = qchisq(0.95, 1),
                                                                       val.sample        = data.validation.cont.cont.mis,
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
                                                                         test.data    = test.data.mis,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("dr-adjF-propscF-cv1")

#####################################################################################################################
##################### 40. dr: True GLM Model in split, True propensity score model in split, Cv2 ####################
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
print("dr-adjT-propscT-cv2")

#####################################################################################################################
################### 41. dr: True GLM Model in split, Mis func propensity score model in split, Cv2 ##################
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
print("dr-adjT-propscNois-cv2")

#####################################################################################################################
################### 42. dr: Mis func GLM Model in split, True propensity score model in split, Cv2 ##################
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
print("dr-adjNois-propscT-cv2")

#####################################################################################################################
################# 43. dr: Mis func GLM Model in split, Mis func propensity score model in split, Cv2 ################
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
print("dr-adjNois-propscNois-cv2")

#####################################################################################################################
########### 44. dr: Unmeasured cov GLM Model in split, Unmeasured cov propensity score model in split, Cv2 ##########
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
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

final.tree.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont.mis,
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
                                                                         test.data    = test.data.mis,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmInsplt.propscFGlmInsplt.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("dr-adjF-propscF-cv2")

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

file.name = paste("../Data/AppendixC5Insplit/", toString(job.number), ".RData", sep = "")
save(performance.hetero.ipw, performance.homo.ipw,
     performance.hetero.g, performance.homo.g,
     performance.hetero.drInsplt, performance.homo.drInsplt,
     file = file.name)
