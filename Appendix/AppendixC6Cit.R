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
data.used.bin.mixed.mis <- data.used.bin.mixed %>%
  select(-X2)
data.validation.bin.mixed.mis <- data.validation.bin.mixed %>%
  select(-X2)
test.data.mis <- data.bin.mixed$test.data %>%
  select(-X2)

#####################################################################################################################
######################### 1. ipw: GLM Model, inside node, True propensity score model, cv1 ##########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.true.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                              est.used          = "IPW",
                                                              type.var          = "bin",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                              w                 = NULL,
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.true.cv1 <- EstIpw.CvMethod1(data.used         = data.used.bin.mixed,
                                                              tree.list         = seq.created.estipw.glm.propscinnd.true.cv1$tree.list,
                                                              lambda.used       = qchisq(0.95, 1),
                                                              val.sample        = data.validation.bin.mixed,
                                                              type.var          = "bin",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                              val.w             = NULL,
                                                              propsc.mod.insplt = F,
                                                              min.obs.mod       = 5)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.true.cv1[[1]],
                                                               test.data    = data.bin.mixed$test.data,
                                                               true.trt.eff = data.bin.mixed$true.trt.eff,
                                                               noise.var    = data.bin.mixed$noise.var,
                                                               corr.split   = data.bin.mixed$corr.split,
                                                               where.split  = data.bin.mixed$where.split,
                                                               dir.split    = data.bin.mixed$dir.split,
                                                               split.cate   = data.bin.mixed$split.cate)
eval.final.estipw.glm.propscinnd.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinnd.true.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estipw.glm.propscinnd.true.cv1$tree.list[[1]],
                                                                                     corr.split = data.bin.mixed$corr.split,
                                                                                     split.cate = data.bin.mixed$split.cate)
print("1")

#####################################################################################################################
######################### 2. ipw: GLM Model, inside node, Mis func propensity score model, cv1 ######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.nois.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                              est.used          = "IPW",
                                                              type.var          = "bin",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                              w                 = NULL,
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.nois.cv1 <- EstIpw.CvMethod1(data.used         = data.used.bin.mixed,
                                                              tree.list         = seq.created.estipw.glm.propscinnd.nois.cv1$tree.list,
                                                              lambda.used       = qchisq(0.95, 1),
                                                              val.sample        = data.validation.bin.mixed,
                                                              type.var          = "bin",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                              val.w             = NULL,
                                                              propsc.mod.insplt = F,
                                                              min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.nois.cv1[[1]],
                                                               test.data    = data.bin.mixed$test.data,
                                                               true.trt.eff = data.bin.mixed$true.trt.eff,
                                                               noise.var    = data.bin.mixed$noise.var,
                                                               corr.split   = data.bin.mixed$corr.split,
                                                               where.split  = data.bin.mixed$where.split,
                                                               dir.split    = data.bin.mixed$dir.split,
                                                               split.cate   = data.bin.mixed$split.cate)
eval.final.estipw.glm.propscinnd.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinnd.nois.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estipw.glm.propscinnd.nois.cv1$tree.list[[1]],
                                                                                     corr.split = data.bin.mixed$corr.split,
                                                                                     split.cate = data.bin.mixed$split.cate)
print("2")

#####################################################################################################################
##################### 3. ipw: GLM Model, inside node, Unmeasured cov propensity score model, cv1 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.mis.cv1 <- create.sequence(data.used         = data.used.bin.mixed.mis,
                                                             est.used          = "IPW",
                                                             type.var          = "bin",
                                                             propsc.mod.out    = F,
                                                             propsc.mthd       = "GLM",
                                                             propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                             w                 = NULL,
                                                             propsc.mod.insplt = F,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscinnd.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.bin.mixed.mis,
                                                             tree.list         = seq.created.estipw.glm.propscinnd.mis.cv1$tree.list,
                                                             lambda.used       = qchisq(0.95, 1),
                                                             val.sample        = data.validation.bin.mixed.mis,
                                                             type.var          = "bin",
                                                             propsc.mod.out    = F,
                                                             propsc.mthd       = "GLM",
                                                             propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                             val.w             = NULL,
                                                             propsc.mod.insplt = F,
                                                             min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.mis.cv1[[1]],
                                                              test.data    = test.data.mis,
                                                              true.trt.eff = data.bin.mixed$true.trt.eff,
                                                              noise.var    = data.bin.mixed$noise.var,
                                                              corr.split   = data.bin.mixed$corr.split,
                                                              where.split  = data.bin.mixed$where.split,
                                                              dir.split    = data.bin.mixed$dir.split,
                                                              split.cate   = data.bin.mixed$split.cate)
eval.final.estipw.glm.propscinnd.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinnd.mis.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estipw.glm.propscinnd.mis.cv1$tree.list[[1]],
                                                                                    corr.split = data.bin.mixed$corr.split,
                                                                                    split.cate = data.bin.mixed$split.cate)
print("3")

#####################################################################################################################
########################## 4. ipw: GLM Model, inside node, True propensity score model, Cv2 #########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.true.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                              est.used          = "IPW",
                                                              type.var          = "bin",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                              w                 = NULL,
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.true.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.bin.mixed,
                                                              tree.list        = seq.created.estipw.glm.propscinnd.true.cv2$tree.list,
                                                              type.var         = "bin",
                                                              seed             = a[job.number],
                                                              n.cv             = 5,
                                                              propsc.mod.out   = F,
                                                              propsc.mthd      = "GLM",
                                                              propsc.form.true = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                              min.obs.mod      = 5)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.true.cv2[[1]],
                                                               test.data    = data.bin.mixed$test.data,
                                                               true.trt.eff = data.bin.mixed$true.trt.eff,
                                                               noise.var    = data.bin.mixed$noise.var,
                                                               corr.split   = data.bin.mixed$corr.split,
                                                               where.split  = data.bin.mixed$where.split,
                                                               dir.split    = data.bin.mixed$dir.split,
                                                               split.cate   = data.bin.mixed$split.cate)
eval.final.estipw.glm.propscinnd.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinnd.true.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estipw.glm.propscinnd.true.cv2$tree.list[[1]],
                                                                                     corr.split = data.bin.mixed$corr.split,
                                                                                     split.cate = data.bin.mixed$split.cate)
print("4")

#####################################################################################################################
######################## 5. ipw: GLM Model, inside node, Mis func propensity score model, Cv2 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.nois.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                              est.used          = "IPW",
                                                              type.var          = "bin",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                              w                 = NULL,
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.nois.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.bin.mixed,
                                                              tree.list        = seq.created.estipw.glm.propscinnd.nois.cv2$tree.list,
                                                              type.var         = "bin",
                                                              seed             = a[job.number],
                                                              n.cv             = 5,
                                                              propsc.mod.out   = F,
                                                              propsc.mthd      = "GLM",
                                                              propsc.form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                              min.obs.mod      = 15)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.nois.cv2[[1]],
                                                               test.data    = data.bin.mixed$test.data,
                                                               true.trt.eff = data.bin.mixed$true.trt.eff,
                                                               noise.var    = data.bin.mixed$noise.var,
                                                               corr.split   = data.bin.mixed$corr.split,
                                                               where.split  = data.bin.mixed$where.split,
                                                               dir.split    = data.bin.mixed$dir.split,
                                                               split.cate   = data.bin.mixed$split.cate)
eval.final.estipw.glm.propscinnd.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinnd.nois.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estipw.glm.propscinnd.nois.cv2$tree.list[[1]],
                                                                                     corr.split = data.bin.mixed$corr.split,
                                                                                     split.cate = data.bin.mixed$split.cate)
print("5")

#####################################################################################################################
##################### 6. ipw: GLM Model, inside node, Unmeasured cov propensity score model, Cv2 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.mis.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed.mis,
                                                             est.used          = "IPW",
                                                             type.var          = "bin",
                                                             propsc.mod.out    = F,
                                                             propsc.mthd       = "GLM",
                                                             propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                             w                 = NULL,
                                                             propsc.mod.insplt = F,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscinnd.mis.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.bin.mixed.mis,
                                                             tree.list        = seq.created.estipw.glm.propscinnd.mis.cv2$tree.list,
                                                             type.var         = "bin",
                                                             seed             = a[job.number],
                                                             n.cv             = 5,
                                                             propsc.mod.out   = F,
                                                             propsc.mthd      = "GLM",
                                                             propsc.form.true = "A ~ X1 + X3 + X4 + X5 + X6",
                                                             min.obs.mod      = 15)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.mis.cv2[[1]],
                                                              test.data    = test.data.mis,
                                                              true.trt.eff = data.bin.mixed$true.trt.eff,
                                                              noise.var    = data.bin.mixed$noise.var,
                                                              corr.split   = data.bin.mixed$corr.split,
                                                              where.split  = data.bin.mixed$where.split,
                                                              dir.split    = data.bin.mixed$dir.split,
                                                              split.cate   = data.bin.mixed$split.cate)
eval.final.estipw.glm.propscinnd.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinnd.mis.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estipw.glm.propscinnd.mis.cv2$tree.list[[1]],
                                                                                    corr.split = data.bin.mixed$corr.split,
                                                                                    split.cate = data.bin.mixed$split.cate)
print("6")

performance.hetero.ipw <- list(glm.propscinnd.true.cv1 = eval.final.estipw.glm.propscinnd.true.cv1,
                               glm.propscinnd.nois.cv1 = eval.final.estipw.glm.propscinnd.nois.cv1,
                               glm.propscinnd.mis.cv1  = eval.final.estipw.glm.propscinnd.mis.cv1,
                               glm.propscinnd.true.cv2 = eval.final.estipw.glm.propscinnd.true.cv2,
                               glm.propscinnd.nois.cv2 = eval.final.estipw.glm.propscinnd.nois.cv2,
                               glm.propscinnd.mis.cv2  = eval.final.estipw.glm.propscinnd.mis.cv2)

#####################################################################################################################
############################# 7. g: GLM Model, inside node, True adjustment model, cv1 ##############################
#####################################################################################################################
t0 <- Sys.time()
# The difference between true adjustment model and the misspecified model is very dependent on the true adjustment model
seq.created.estg.glm.modinnd.true.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                         est.used          = "G",
                                                         type.var          = "bin",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))",
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 10,
                                                         min.node          = 10)

final.tree.estg.glm.modinnd.true.cv1 <- EstG.CvMethod1(data.used         = data.used.bin.mixed,
                                                       tree.list         = seq.created.estg.glm.modinnd.true.cv1$tree.list,
                                                       lambda.used       = qchisq(0.95, 1),
                                                       val.sample        = data.validation.bin.mixed,              
                                                       type.var          = "bin",
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))",
                                                       val.w             = NULL,
                                                       adj.mod.insplt    = F,
                                                       min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.true.cv1[[1]],
                                                          test.data    = data.bin.mixed$test.data,
                                                          true.trt.eff = data.bin.mixed$true.trt.eff,
                                                          noise.var    = data.bin.mixed$noise.var,
                                                          corr.split   = data.bin.mixed$corr.split,
                                                          where.split  = data.bin.mixed$where.split,
                                                          dir.split    = data.bin.mixed$dir.split,
                                                          split.cate   = data.bin.mixed$split.cate)
eval.final.estg.glm.modinnd.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinnd.true.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estg.glm.modinnd.true.cv1$tree.list[[1]],
                                                                                corr.split = data.bin.mixed$corr.split,
                                                                                split.cate = data.bin.mixed$split.cate)
print("7")

#####################################################################################################################
########################### 8. g: GLM Model, inside node, Mis func adjustment model, cv1 ############################
#####################################################################################################################
t0 <- Sys.time()
# The difference between true adjustment model and the misspecified model is very dependent on the true adjustment model
seq.created.estg.glm.modinnd.nois.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                         est.used          = "G",
                                                         type.var          = "bin",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = NULL,
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 15,
                                                         min.node          = 15)

final.tree.estg.glm.modinnd.nois.cv1 <- EstG.CvMethod1(data.used         = data.used.bin.mixed,
                                                       tree.list         = seq.created.estg.glm.modinnd.nois.cv1$tree.list,
                                                       lambda.used       = qchisq(0.95, 1),
                                                       val.sample        = data.validation.bin.mixed,              
                                                       type.var          = "bin",
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = NULL,
                                                       val.w             = NULL,
                                                       adj.mod.insplt    = F,
                                                       min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.nois.cv1[[1]],
                                                          test.data    = data.bin.mixed$test.data,
                                                          true.trt.eff = data.bin.mixed$true.trt.eff,
                                                          noise.var    = data.bin.mixed$noise.var,
                                                          corr.split   = data.bin.mixed$corr.split,
                                                          where.split  = data.bin.mixed$where.split,
                                                          dir.split    = data.bin.mixed$dir.split,
                                                          split.cate   = data.bin.mixed$split.cate)
eval.final.estg.glm.modinnd.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinnd.nois.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estg.glm.modinnd.nois.cv1$tree.list[[1]],
                                                                                corr.split = data.bin.mixed$corr.split,
                                                                                split.cate = data.bin.mixed$split.cate)
print("8")

#####################################################################################################################
########################## 9. g: GLM Model, inside node, Unmeasured cov adjustment model, cv1 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.mis.cv1 <- create.sequence(data.used         = data.used.bin.mixed.mis,
                                                        est.used          = "G",
                                                        type.var          = "bin",
                                                        adj.mod.out       = F,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                        w                 = NULL,
                                                        adj.mod.insplt    = F,
                                                        num.truc.obs      = 15,
                                                        min.node          = 15)

final.tree.estg.glm.modinnd.mis.cv1 <- EstG.CvMethod1(data.used         = data.used.bin.mixed.mis,
                                                      tree.list         = seq.created.estg.glm.modinnd.mis.cv1$tree.list,
                                                      lambda.used       = qchisq(0.95, 1),
                                                      val.sample        = data.validation.bin.mixed.mis,              
                                                      type.var          = "bin",
                                                      adj.mod.out       = F,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                      val.w             = NULL,
                                                      adj.mod.insplt    = F,
                                                      min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.mis.cv1[[1]],
                                                         test.data    = test.data.mis,
                                                         true.trt.eff = data.bin.mixed$true.trt.eff,
                                                         noise.var    = data.bin.mixed$noise.var,
                                                         corr.split   = data.bin.mixed$corr.split,
                                                         where.split  = data.bin.mixed$where.split,
                                                         dir.split    = data.bin.mixed$dir.split,
                                                         split.cate   = data.bin.mixed$split.cate)
eval.final.estg.glm.modinnd.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinnd.mis.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estg.glm.modinnd.mis.cv1$tree.list[[1]],
                                                                               corr.split = data.bin.mixed$corr.split,
                                                                               split.cate = data.bin.mixed$split.cate)
print("9")

#####################################################################################################################
############################# 10. g: GLM Model, inside node, True adjustment model, Cv2 #############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.true.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                         est.used          = "G",
                                                         type.var          = "bin",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))",
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 10,
                                                         min.node          = 10)

final.tree.estg.glm.modinnd.true.cv2 <- EstG.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                       tree.list         = seq.created.estg.glm.modinnd.true.cv2$tree.list,             
                                                       type.var          = "bin",
                                                       seed              = a[job.number],
                                                       n.cv              = 5,
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))",
                                                       min.obs.mod       = 10)

t1 <- Sys.time()

eval.final.estg.glm.modinnd.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.true.cv2[[1]],
                                                          test.data    = data.bin.mixed$test.data,
                                                          true.trt.eff = data.bin.mixed$true.trt.eff,
                                                          noise.var    = data.bin.mixed$noise.var,
                                                          corr.split   = data.bin.mixed$corr.split,
                                                          where.split  = data.bin.mixed$where.split,
                                                          dir.split    = data.bin.mixed$dir.split,
                                                          split.cate   = data.bin.mixed$split.cate)
eval.final.estg.glm.modinnd.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinnd.true.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estg.glm.modinnd.true.cv2$tree.list[[1]],
                                                                                corr.split = data.bin.mixed$corr.split,
                                                                                split.cate = data.bin.mixed$split.cate)
print("10")

#####################################################################################################################
############################ 11. g: GLM Model, inside node, Mis func adjustment model, Cv2 ##########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.nois.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                         est.used          = "G",
                                                         type.var          = "bin",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = NULL,
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 15,
                                                         min.node          = 15)

final.tree.estg.glm.modinnd.nois.cv2 <- EstG.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                       tree.list         = seq.created.estg.glm.modinnd.nois.cv2$tree.list,             
                                                       type.var          = "bin",
                                                       seed              = a[job.number],
                                                       n.cv              = 5,
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = NULL,
                                                       min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.nois.cv2[[1]],
                                                          test.data    = data.bin.mixed$test.data,
                                                          true.trt.eff = data.bin.mixed$true.trt.eff,
                                                          noise.var    = data.bin.mixed$noise.var,
                                                          corr.split   = data.bin.mixed$corr.split,
                                                          where.split  = data.bin.mixed$where.split,
                                                          dir.split    = data.bin.mixed$dir.split,
                                                          split.cate   = data.bin.mixed$split.cate)
eval.final.estg.glm.modinnd.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinnd.nois.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estg.glm.modinnd.nois.cv2$tree.list[[1]],
                                                                                corr.split = data.bin.mixed$corr.split,
                                                                                split.cate = data.bin.mixed$split.cate)
print("11")

#####################################################################################################################
######################### 12. g: GLM Model, inside node, Unmeasured cov adjustment model, Cv2 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.mis.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed.mis,
                                                        est.used          = "G",
                                                        type.var          = "bin",
                                                        adj.mod.out       = F,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                        w                 = NULL,
                                                        adj.mod.insplt    = F,
                                                        num.truc.obs      = 15,
                                                        min.node          = 15)

final.tree.estg.glm.modinnd.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.bin.mixed.mis,
                                                      tree.list         = seq.created.estg.glm.modinnd.mis.cv2$tree.list,             
                                                      type.var          = "bin",
                                                      seed              = a[job.number],
                                                      n.cv              = 5,
                                                      adj.mod.out       = F,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                      min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.mis.cv2[[1]],
                                                         test.data    = test.data.mis,
                                                         true.trt.eff = data.bin.mixed$true.trt.eff,
                                                         noise.var    = data.bin.mixed$noise.var,
                                                         corr.split   = data.bin.mixed$corr.split,
                                                         where.split  = data.bin.mixed$where.split,
                                                         dir.split    = data.bin.mixed$dir.split,
                                                         split.cate   = data.bin.mixed$split.cate)
eval.final.estg.glm.modinnd.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinnd.mis.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estg.glm.modinnd.mis.cv2$tree.list[[1]],
                                                                               corr.split = data.bin.mixed$corr.split,
                                                                               split.cate = data.bin.mixed$split.cate)
print("12")

performance.hetero.g <- list(glm.modinnd.true.cv1 = eval.final.estg.glm.modinnd.true.cv1,
                             glm.modinnd.nois.cv1 = eval.final.estg.glm.modinnd.nois.cv1,
                             glm.modinnd.mis.cv1  = eval.final.estg.glm.modinnd.mis.cv1,
                             glm.modinnd.true.cv2 = eval.final.estg.glm.modinnd.true.cv2,
                             glm.modinnd.nois.cv2 = eval.final.estg.glm.modinnd.nois.cv2,
                             glm.modinnd.mis.cv2  = eval.final.estg.glm.modinnd.mis.cv2)

#####################################################################################################################
##################### 13. dr: True GLM Model in node, True propensity score model in node, Cv1 ######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                    est.used          = "DR",
                                                                    type.var          = "bin",
                                                                    propsc.mod.out    = F,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                    propsc.mod.insplt = F,
                                                                    adj.mod.out       = F, 
                                                                    adj.mthd          = "GLM", 
                                                                    adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))", 
                                                                    adj.mod.insplt    = F, 
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)

final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                   tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                   lambda.used       = qchisq(0.95, 1),
                                                                   val.sample        = data.validation.bin.mixed,
                                                                   type.var          = "bin",
                                                                   propsc.mod.out    = F,
                                                                   propsc.mthd       = "GLM",
                                                                   propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                   propsc.mod.insplt = F,
                                                                   adj.mod.out       = F,
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))", 
                                                                   adj.mod.insplt    = F,
                                                                   min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1[[1]],
                                                                     test.data    = data.bin.mixed$test.data,
                                                                     true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                     noise.var    = data.bin.mixed$noise.var,
                                                                     corr.split   = data.bin.mixed$corr.split,
                                                                     where.split  = data.bin.mixed$where.split,
                                                                     dir.split    = data.bin.mixed$dir.split,
                                                                     split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1$tree.list[[1]],
                                                                                           corr.split = data.bin.mixed$corr.split,
                                                                                           split.cate = data.bin.mixed$split.cate)
print("13")

#####################################################################################################################
##################### 14. dr: True GLM Model in node, Mis func propensity score model in node, Cv1 ##################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                       est.used          = "DR",
                                                                       type.var          = "bin",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                       propsc.mod.insplt = F,
                                                                       adj.mod.out       = F, 
                                                                       adj.mthd          = "GLM", 
                                                                       adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))", 
                                                                       adj.mod.insplt    = F, 
                                                                       num.truc.obs      = 30,
                                                                       min.node          = 20)

final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                      tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                      lambda.used       = qchisq(0.95, 1),
                                                                      val.sample        = data.validation.bin.mixed,
                                                                      type.var          = "bin",
                                                                      propsc.mod.out    = F,
                                                                      propsc.mthd       = "GLM",
                                                                      propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                      propsc.mod.insplt = F,
                                                                      adj.mod.out       = F,
                                                                      adj.mthd          = "GLM",
                                                                      adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))", 
                                                                      adj.mod.insplt    = F,
                                                                      min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1[[1]],
                                                                        test.data    = data.bin.mixed$test.data,
                                                                        true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                        noise.var    = data.bin.mixed$noise.var,
                                                                        corr.split   = data.bin.mixed$corr.split,
                                                                        where.split  = data.bin.mixed$where.split,
                                                                        dir.split    = data.bin.mixed$dir.split,
                                                                        split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$tree.list[[1]],
                                                                                              corr.split = data.bin.mixed$corr.split,
                                                                                              split.cate = data.bin.mixed$split.cate)
print("14")

#####################################################################################################################
#################### 15. dr: Mis func GLM Model in node, True propensity score model in node, Cv1 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                       est.used          = "DR",
                                                                       type.var          = "bin",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                       propsc.mod.insplt = F,
                                                                       adj.mod.out       = F, 
                                                                       adj.mthd          = "GLM", 
                                                                       adj.form.true     = NULL, 
                                                                       adj.mod.insplt    = F, 
                                                                       num.truc.obs      = 30,
                                                                       min.node          = 20)

final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                      tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                      lambda.used       = qchisq(0.95, 1),
                                                                      val.sample        = data.validation.bin.mixed,
                                                                      type.var          = "bin",
                                                                      propsc.mod.out    = F,
                                                                      propsc.mthd       = "GLM",
                                                                      propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                      propsc.mod.insplt = F,
                                                                      adj.mod.out       = F,
                                                                      adj.mthd          = "GLM",
                                                                      adj.form.true     = NULL, 
                                                                      adj.mod.insplt    = F,
                                                                      min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1[[1]],
                                                                        test.data    = data.bin.mixed$test.data,
                                                                        true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                        noise.var    = data.bin.mixed$noise.var,
                                                                        corr.split   = data.bin.mixed$corr.split,
                                                                        where.split  = data.bin.mixed$where.split,
                                                                        dir.split    = data.bin.mixed$dir.split,
                                                                        split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$tree.list[[1]],
                                                                                              corr.split = data.bin.mixed$corr.split,
                                                                                              split.cate = data.bin.mixed$split.cate)
print("15")

#####################################################################################################################
################## 16. dr: Mis func GLM Model in node, Mis func propensity score model in node, Cv1 #################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                          est.used          = "DR",
                                                                          type.var          = "bin",
                                                                          propsc.mod.out    = F,
                                                                          propsc.mthd       = "GLM",
                                                                          propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                          propsc.mod.insplt = F,
                                                                          adj.mod.out       = F, 
                                                                          adj.mthd          = "GLM", 
                                                                          adj.form.true     = NULL, 
                                                                          adj.mod.insplt    = F, 
                                                                          num.truc.obs      = 30,
                                                                          min.node          = 20)

final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                         tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                         lambda.used       = qchisq(0.95, 1),
                                                                         val.sample        = data.validation.bin.mixed,
                                                                         type.var          = "bin",
                                                                         propsc.mod.out    = F,
                                                                         propsc.mthd       = "GLM",
                                                                         propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                         propsc.mod.insplt = F,
                                                                         adj.mod.out       = F,
                                                                         adj.mthd          = "GLM",
                                                                         adj.form.true     = NULL, 
                                                                         adj.mod.insplt    = F,
                                                                         min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1[[1]],
                                                                           test.data    = data.bin.mixed$test.data,
                                                                           true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                           noise.var    = data.bin.mixed$noise.var,
                                                                           corr.split   = data.bin.mixed$corr.split,
                                                                           where.split  = data.bin.mixed$where.split,
                                                                           dir.split    = data.bin.mixed$dir.split,
                                                                           split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$tree.list[[1]],
                                                                                                 corr.split = data.bin.mixed$corr.split,
                                                                                                 split.cate = data.bin.mixed$split.cate)
print("16")

#####################################################################################################################
########### 17. dr: Unmeasured cov GLM Model in node, Unmeasured cov propensity score model in node, Cv1 ############
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed.mis,
                                                                    est.used          = "DR",
                                                                    type.var          = "bin",
                                                                    propsc.mod.out    = F,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                    propsc.mod.insplt = F,
                                                                    adj.mod.out       = F, 
                                                                    adj.mthd          = "GLM", 
                                                                    adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                    adj.mod.insplt    = F, 
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)

final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed.mis,
                                                                   tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1$tree.list, 
                                                                   lambda.used       = qchisq(0.95, 1),
                                                                   val.sample        = data.validation.bin.mixed.mis,
                                                                   type.var          = "bin",
                                                                   propsc.mod.out    = F,
                                                                   propsc.mthd       = "GLM",
                                                                   propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                   propsc.mod.insplt = F,
                                                                   adj.mod.out       = F,
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                   adj.mod.insplt    = F,
                                                                   min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1[[1]],
                                                                     test.data    = test.data.mis,
                                                                     true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                     noise.var    = data.bin.mixed$noise.var,
                                                                     corr.split   = data.bin.mixed$corr.split,
                                                                     where.split  = data.bin.mixed$where.split,
                                                                     dir.split    = data.bin.mixed$dir.split,
                                                                     split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1$tree.list[[1]],
                                                                                           corr.split = data.bin.mixed$corr.split,
                                                                                           split.cate = data.bin.mixed$split.cate)
print("17")

#####################################################################################################################
###################### 18. dr: True GLM Model in node, True propensity score model in node, Cv2 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                    est.used          = "DR",
                                                                    type.var          = "bin",
                                                                    propsc.mod.out    = F,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                    propsc.mod.insplt = F,
                                                                    adj.mod.out       = F, 
                                                                    adj.mthd          = "GLM", 
                                                                    adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))", 
                                                                    adj.mod.insplt    = F, 
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)

final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                   tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2$tree.list,
                                                                   type.var          = "bin",
                                                                   seed              = a[job.number], 
                                                                   n.cv              = 5,
                                                                   propsc.mod.out    = F,
                                                                   propsc.mthd       = "GLM",
                                                                   propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                   min.obs.propsc    = 5,
                                                                   adj.mod.out       = F,
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))", 
                                                                   min.obs.adj       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2[[1]],
                                                                     test.data    = data.bin.mixed$test.data,
                                                                     true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                     noise.var    = data.bin.mixed$noise.var,
                                                                     corr.split   = data.bin.mixed$corr.split,
                                                                     where.split  = data.bin.mixed$where.split,
                                                                     dir.split    = data.bin.mixed$dir.split,
                                                                     split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2$tree.list[[1]],
                                                                                           corr.split = data.bin.mixed$corr.split,
                                                                                           split.cate = data.bin.mixed$split.cate)
print("18")

#####################################################################################################################
#################### 19. dr: True GLM Model in node, Mis func propensity score model in node, Cv2 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                       est.used          = "DR",
                                                                       type.var          = "bin",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                       propsc.mod.insplt = F,
                                                                       adj.mod.out       = F, 
                                                                       adj.mthd          = "GLM", 
                                                                       adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))", 
                                                                       adj.mod.insplt    = F, 
                                                                       num.truc.obs      = 30,
                                                                       min.node          = 20)

final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                      tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$tree.list,
                                                                      type.var          = "bin",
                                                                      seed              = a[job.number], 
                                                                      n.cv              = 5,
                                                                      propsc.mod.out    = F,
                                                                      propsc.mthd       = "GLM",
                                                                      propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                      min.obs.propsc    = 10,
                                                                      adj.mod.out       = F,
                                                                      adj.mthd          = "GLM",
                                                                      adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))", 
                                                                      min.obs.adj       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2[[1]],
                                                                        test.data    = data.bin.mixed$test.data,
                                                                        true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                        noise.var    = data.bin.mixed$noise.var,
                                                                        corr.split   = data.bin.mixed$corr.split,
                                                                        where.split  = data.bin.mixed$where.split,
                                                                        dir.split    = data.bin.mixed$dir.split,
                                                                        split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$tree.list[[1]],
                                                                                              corr.split = data.bin.mixed$corr.split,
                                                                                              split.cate = data.bin.mixed$split.cate)
print("19")

#####################################################################################################################
#################### 20. dr: Mis func GLM Model in node, True propensity score model in node, Cv2 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                       est.used          = "DR",
                                                                       type.var          = "bin",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                       propsc.mod.insplt = F,
                                                                       adj.mod.out       = F, 
                                                                       adj.mthd          = "GLM", 
                                                                       adj.form.true     = NULL, 
                                                                       adj.mod.insplt    = F, 
                                                                       num.truc.obs      = 30,
                                                                       min.node          = 20)

final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                      tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$tree.list,
                                                                      type.var          = "bin",
                                                                      seed              = a[job.number], 
                                                                      n.cv              = 5,
                                                                      propsc.mod.out    = F,
                                                                      propsc.mthd       = "GLM",
                                                                      propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                      min.obs.propsc    = 5,
                                                                      adj.mod.out       = F,
                                                                      adj.mthd          = "GLM",
                                                                      adj.form.true     = NULL, 
                                                                      min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2[[1]],
                                                                        test.data    = data.bin.mixed$test.data,
                                                                        true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                        noise.var    = data.bin.mixed$noise.var,
                                                                        corr.split   = data.bin.mixed$corr.split,
                                                                        where.split  = data.bin.mixed$where.split,
                                                                        dir.split    = data.bin.mixed$dir.split,
                                                                        split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$tree.list[[1]],
                                                                                              corr.split = data.bin.mixed$corr.split,
                                                                                              split.cate = data.bin.mixed$split.cate)
print("20")

#####################################################################################################################
################# 21. dr: Mis func GLM Model in node, Mis func propensity score model in node, Cv2 ##################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                          est.used          = "DR",
                                                                          type.var          = "bin",
                                                                          propsc.mod.out    = F,
                                                                          propsc.mthd       = "GLM",
                                                                          propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                          propsc.mod.insplt = F,
                                                                          adj.mod.out       = F, 
                                                                          adj.mthd          = "GLM", 
                                                                          adj.form.true     = NULL, 
                                                                          adj.mod.insplt    = F, 
                                                                          num.truc.obs      = 30,
                                                                          min.node          = 20)

final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                         tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$tree.list,
                                                                         type.var          = "bin",
                                                                         seed              = a[job.number], 
                                                                         n.cv              = 5,
                                                                         propsc.mod.out    = F,
                                                                         propsc.mthd       = "GLM",
                                                                         propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                         min.obs.propsc    = 10,
                                                                         adj.mod.out       = F,
                                                                         adj.mthd          = "GLM",
                                                                         adj.form.true     = NULL, 
                                                                         min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2[[1]],
                                                                           test.data    = data.bin.mixed$test.data,
                                                                           true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                           noise.var    = data.bin.mixed$noise.var,
                                                                           corr.split   = data.bin.mixed$corr.split,
                                                                           where.split  = data.bin.mixed$where.split,
                                                                           dir.split    = data.bin.mixed$dir.split,
                                                                           split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$tree.list[[1]],
                                                                                                 corr.split = data.bin.mixed$corr.split,
                                                                                                 split.cate = data.bin.mixed$split.cate)
print("21")

#####################################################################################################################
############ 22. dr: Unmeasured cov GLM Model in node, Unmeasured cov propensity score model in node, Cv2 ###########
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed.mis,
                                                                    est.used          = "DR",
                                                                    type.var          = "bin",
                                                                    propsc.mod.out    = F,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                    propsc.mod.insplt = F,
                                                                    adj.mod.out       = F, 
                                                                    adj.mthd          = "GLM", 
                                                                    adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                    adj.mod.insplt    = F, 
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)

final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed.mis,
                                                                   tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2$tree.list,
                                                                   type.var          = "bin",
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

eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv2[[1]],
                                                                     test.data    = test.data.mis,
                                                                     true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                     noise.var    = data.bin.mixed$noise.var,
                                                                     corr.split   = data.bin.mixed$corr.split,
                                                                     where.split  = data.bin.mixed$where.split,
                                                                     dir.split    = data.bin.mixed$dir.split,
                                                                     split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2$tree.list[[1]],
                                                                                           corr.split = data.bin.mixed$corr.split,
                                                                                           split.cate = data.bin.mixed$split.cate)
print("22")

performance.hetero.drInnd <- list(adjTGlmInnd.propscTGlmInnd.cv1       = eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1,
                                  adjTGlmInnd.propscNoisGlmInnd.cv1    = eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1,
                                  adjNoisGlmInnd.propscTGlmInnd.cv1    = eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1,
                                  adjNoisGlmInnd.propscNoisGlmInnd.cv1 = eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1,
                                  adjFGlmInnd.propscFGlmInnd.cv1       = eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1,
                                  adjTGlmInnd.propscTGlmInnd.cv2       = eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2,
                                  adjTGlmInnd.propscNoisGlmInnd.cv2    = eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2,
                                  adjNoisGlmInnd.propscTGlmInnd.cv2    = eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2,
                                  adjNoisGlmInnd.propscNoisGlmInnd.cv2 = eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2,
                                  adjFGlmInnd.propscFGlmInnd.cv2       = eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)

#####################################################################################################################
################################################### Homogeneous #####################################################
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
data.used.bin.mixed.mis <- data.used.bin.mixed %>%
  select(-X2)
data.validation.bin.mixed.mis <- data.validation.bin.mixed %>%
  select(-X2)
test.data.mis <- data.bin.mixed$test.data %>%
  select(-X2)

#####################################################################################################################
######################### 23. ipw: GLM Model, inside node, True propensity score model, cv1 #########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.true.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                              est.used          = "IPW",
                                                              type.var          = "bin",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                              w                 = NULL,
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.true.cv1 <- EstIpw.CvMethod1(data.used         = data.used.bin.mixed,
                                                              tree.list         = seq.created.estipw.glm.propscinnd.true.cv1$tree.list,
                                                              lambda.used       = qchisq(0.95, 1),
                                                              val.sample        = data.validation.bin.mixed,
                                                              type.var          = "bin",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                              val.w             = NULL,
                                                              propsc.mod.insplt = F,
                                                              min.obs.mod       = 5)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.true.cv1[[1]],
                                                               test.data    = data.bin.mixed$test.data,
                                                               true.trt.eff = data.bin.mixed$true.trt.eff,
                                                               noise.var    = data.bin.mixed$noise.var,
                                                               corr.split   = data.bin.mixed$corr.split,
                                                               where.split  = data.bin.mixed$where.split,
                                                               dir.split    = data.bin.mixed$dir.split,
                                                               split.cate   = data.bin.mixed$split.cate)
eval.final.estipw.glm.propscinnd.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("23")

#####################################################################################################################
####################### 24. ipw: GLM Model, inside node, Mis func propensity score model, cv1 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.nois.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                              est.used          = "IPW",
                                                              type.var          = "bin",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                              w                 = NULL,
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.nois.cv1 <- EstIpw.CvMethod1(data.used         = data.used.bin.mixed,
                                                              tree.list         = seq.created.estipw.glm.propscinnd.nois.cv1$tree.list,
                                                              lambda.used       = qchisq(0.95, 1),
                                                              val.sample        = data.validation.bin.mixed,
                                                              type.var          = "bin",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                              val.w             = NULL,
                                                              propsc.mod.insplt = F,
                                                              min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.nois.cv1[[1]],
                                                               test.data    = data.bin.mixed$test.data,
                                                               true.trt.eff = data.bin.mixed$true.trt.eff,
                                                               noise.var    = data.bin.mixed$noise.var,
                                                               corr.split   = data.bin.mixed$corr.split,
                                                               where.split  = data.bin.mixed$where.split,
                                                               dir.split    = data.bin.mixed$dir.split,
                                                               split.cate   = data.bin.mixed$split.cate)
eval.final.estipw.glm.propscinnd.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("24")

#####################################################################################################################
#################### 25. ipw: GLM Model, inside node, Unmeasured cov propensity score model, cv1 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.mis.cv1 <- create.sequence(data.used         = data.used.bin.mixed.mis,
                                                             est.used          = "IPW",
                                                             type.var          = "bin",
                                                             propsc.mod.out    = F,
                                                             propsc.mthd       = "GLM",
                                                             propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                             w                 = NULL,
                                                             propsc.mod.insplt = F,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscinnd.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.bin.mixed.mis,
                                                             tree.list         = seq.created.estipw.glm.propscinnd.mis.cv1$tree.list,
                                                             lambda.used       = qchisq(0.95, 1),
                                                             val.sample        = data.validation.bin.mixed.mis,
                                                             type.var          = "bin",
                                                             propsc.mod.out    = F,
                                                             propsc.mthd       = "GLM",
                                                             propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                             val.w             = NULL,
                                                             propsc.mod.insplt = F,
                                                             min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.mis.cv1[[1]],
                                                              test.data    = test.data.mis,
                                                              true.trt.eff = data.bin.mixed$true.trt.eff,
                                                              noise.var    = data.bin.mixed$noise.var,
                                                              corr.split   = data.bin.mixed$corr.split,
                                                              where.split  = data.bin.mixed$where.split,
                                                              dir.split    = data.bin.mixed$dir.split,
                                                              split.cate   = data.bin.mixed$split.cate)
eval.final.estipw.glm.propscinnd.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("25")

#####################################################################################################################
######################### 26. ipw: GLM Model, inside node, True propensity score model, Cv2 #########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.true.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                              est.used          = "IPW",
                                                              type.var          = "bin",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                              w                 = NULL,
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.true.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.bin.mixed,
                                                              tree.list        = seq.created.estipw.glm.propscinnd.true.cv2$tree.list,
                                                              type.var         = "bin",
                                                              seed             = a[job.number],
                                                              n.cv             = 5,
                                                              propsc.mod.out   = F,
                                                              propsc.mthd      = "GLM",
                                                              propsc.form.true = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                              min.obs.mod      = 5)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.true.cv2[[1]],
                                                               test.data    = data.bin.mixed$test.data,
                                                               true.trt.eff = data.bin.mixed$true.trt.eff,
                                                               noise.var    = data.bin.mixed$noise.var,
                                                               corr.split   = data.bin.mixed$corr.split,
                                                               where.split  = data.bin.mixed$where.split,
                                                               dir.split    = data.bin.mixed$dir.split,
                                                               split.cate   = data.bin.mixed$split.cate)
eval.final.estipw.glm.propscinnd.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("26")

#####################################################################################################################
####################### 27. ipw: GLM Model, inside node, Mis func propensity score model, Cv2 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.nois.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                              est.used          = "IPW",
                                                              type.var          = "bin",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                              w                 = NULL,
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.nois.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.bin.mixed,
                                                              tree.list        = seq.created.estipw.glm.propscinnd.nois.cv2$tree.list,
                                                              type.var         = "bin",
                                                              seed             = a[job.number],
                                                              n.cv             = 5,
                                                              propsc.mod.out   = F,
                                                              propsc.mthd      = "GLM",
                                                              propsc.form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                              min.obs.mod      = 15)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.nois.cv2[[1]],
                                                               test.data    = data.bin.mixed$test.data,
                                                               true.trt.eff = data.bin.mixed$true.trt.eff,
                                                               noise.var    = data.bin.mixed$noise.var,
                                                               corr.split   = data.bin.mixed$corr.split,
                                                               where.split  = data.bin.mixed$where.split,
                                                               dir.split    = data.bin.mixed$dir.split,
                                                               split.cate   = data.bin.mixed$split.cate)
eval.final.estipw.glm.propscinnd.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("27")

#####################################################################################################################
###################### 28. ipw: GLM Model, inside node, Unmeasured cov propensity score model, Cv2 ##################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.mis.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed.mis,
                                                             est.used          = "IPW",
                                                             type.var          = "bin",
                                                             propsc.mod.out    = F,
                                                             propsc.mthd       = "GLM",
                                                             propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                             w                 = NULL,
                                                             propsc.mod.insplt = F,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscinnd.mis.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.bin.mixed.mis,
                                                             tree.list        = seq.created.estipw.glm.propscinnd.mis.cv2$tree.list,
                                                             type.var         = "bin",
                                                             seed             = a[job.number],
                                                             n.cv             = 5,
                                                             propsc.mod.out   = F,
                                                             propsc.mthd      = "GLM",
                                                             propsc.form.true = "A ~ X1 + X3 + X4 + X5 + X6",
                                                             min.obs.mod      = 15)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.mis.cv2[[1]],
                                                              test.data    = test.data.mis,
                                                              true.trt.eff = data.bin.mixed$true.trt.eff,
                                                              noise.var    = data.bin.mixed$noise.var,
                                                              corr.split   = data.bin.mixed$corr.split,
                                                              where.split  = data.bin.mixed$where.split,
                                                              dir.split    = data.bin.mixed$dir.split,
                                                              split.cate   = data.bin.mixed$split.cate)
eval.final.estipw.glm.propscinnd.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("28")

performance.homo.ipw <- list(glm.propscinnd.true.cv1 = eval.final.estipw.glm.propscinnd.true.cv1,
                             glm.propscinnd.nois.cv1 = eval.final.estipw.glm.propscinnd.nois.cv1,
                             glm.propscinnd.mis.cv1  = eval.final.estipw.glm.propscinnd.mis.cv1,
                             glm.propscinnd.true.cv2 = eval.final.estipw.glm.propscinnd.true.cv2,
                             glm.propscinnd.nois.cv2 = eval.final.estipw.glm.propscinnd.nois.cv2,
                             glm.propscinnd.mis.cv2  = eval.final.estipw.glm.propscinnd.mis.cv2)

#####################################################################################################################
############################ 29. g: GLM Model, inside node, True adjustment model, cv1 ##############################
#####################################################################################################################
t0 <- Sys.time()
# The difference between true adjustment model and the misspecified model is very dependent on the true adjustment model
seq.created.estg.glm.modinnd.true.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                         est.used          = "G",
                                                         type.var          = "bin",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))",
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 10,
                                                         min.node          = 10)

final.tree.estg.glm.modinnd.true.cv1 <- EstG.CvMethod1(data.used         = data.used.bin.mixed,
                                                       tree.list         = seq.created.estg.glm.modinnd.true.cv1$tree.list,
                                                       lambda.used       = qchisq(0.95, 1),
                                                       val.sample        = data.validation.bin.mixed,              
                                                       type.var          = "bin",
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))",
                                                       val.w             = NULL,
                                                       adj.mod.insplt    = F,
                                                       min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.true.cv1[[1]],
                                                          test.data    = data.bin.mixed$test.data,
                                                          true.trt.eff = data.bin.mixed$true.trt.eff,
                                                          noise.var    = data.bin.mixed$noise.var,
                                                          corr.split   = data.bin.mixed$corr.split,
                                                          where.split  = data.bin.mixed$where.split,
                                                          dir.split    = data.bin.mixed$dir.split,
                                                          split.cate   = data.bin.mixed$split.cate)
eval.final.estg.glm.modinnd.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("29")

#####################################################################################################################
########################### 30. g: GLM Model, inside node, Mis func adjustment model, cv1 ###########################
#####################################################################################################################
t0 <- Sys.time()
# The difference between true adjustment model and the misspecified model is very dependent on the true adjustment model
seq.created.estg.glm.modinnd.nois.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                         est.used          = "G",
                                                         type.var          = "bin",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = NULL,
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 15,
                                                         min.node          = 15)

final.tree.estg.glm.modinnd.nois.cv1 <- EstG.CvMethod1(data.used         = data.used.bin.mixed,
                                                       tree.list         = seq.created.estg.glm.modinnd.nois.cv1$tree.list,
                                                       lambda.used       = qchisq(0.95, 1),
                                                       val.sample        = data.validation.bin.mixed,              
                                                       type.var          = "bin",
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = NULL,
                                                       val.w             = NULL,
                                                       adj.mod.insplt    = F,
                                                       min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.nois.cv1[[1]],
                                                          test.data    = data.bin.mixed$test.data,
                                                          true.trt.eff = data.bin.mixed$true.trt.eff,
                                                          noise.var    = data.bin.mixed$noise.var,
                                                          corr.split   = data.bin.mixed$corr.split,
                                                          where.split  = data.bin.mixed$where.split,
                                                          dir.split    = data.bin.mixed$dir.split,
                                                          split.cate   = data.bin.mixed$split.cate)
eval.final.estg.glm.modinnd.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("30")

#####################################################################################################################
######################### 31. g: GLM Model, inside node, Unmeasured cov adjustment model, cv1 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.mis.cv1 <- create.sequence(data.used         = data.used.bin.mixed.mis,
                                                        est.used          = "G",
                                                        type.var          = "bin",
                                                        adj.mod.out       = F,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                        w                 = NULL,
                                                        adj.mod.insplt    = F,
                                                        num.truc.obs      = 15,
                                                        min.node          = 15)

final.tree.estg.glm.modinnd.mis.cv1 <- EstG.CvMethod1(data.used         = data.used.bin.mixed.mis,
                                                      tree.list         = seq.created.estg.glm.modinnd.mis.cv1$tree.list,
                                                      lambda.used       = qchisq(0.95, 1),
                                                      val.sample        = data.validation.bin.mixed.mis,              
                                                      type.var          = "bin",
                                                      adj.mod.out       = F,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                      val.w             = NULL,
                                                      adj.mod.insplt    = F,
                                                      min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.mis.cv1[[1]],
                                                         test.data    = test.data.mis,
                                                         true.trt.eff = data.bin.mixed$true.trt.eff,
                                                         noise.var    = data.bin.mixed$noise.var,
                                                         corr.split   = data.bin.mixed$corr.split,
                                                         where.split  = data.bin.mixed$where.split,
                                                         dir.split    = data.bin.mixed$dir.split,
                                                         split.cate   = data.bin.mixed$split.cate)
eval.final.estg.glm.modinnd.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("31")

#####################################################################################################################
############################# 32. g: GLM Model, inside node, True adjustment model, Cv2 #############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.true.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                         est.used          = "G",
                                                         type.var          = "bin",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))",
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 10,
                                                         min.node          = 10)

final.tree.estg.glm.modinnd.true.cv2 <- EstG.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                       tree.list         = seq.created.estg.glm.modinnd.true.cv2$tree.list,             
                                                       type.var          = "bin",
                                                       seed              = a[job.number],
                                                       n.cv              = 5,
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))",
                                                       min.obs.mod       = 10)

t1 <- Sys.time()

eval.final.estg.glm.modinnd.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.true.cv2[[1]],
                                                          test.data    = data.bin.mixed$test.data,
                                                          true.trt.eff = data.bin.mixed$true.trt.eff,
                                                          noise.var    = data.bin.mixed$noise.var,
                                                          corr.split   = data.bin.mixed$corr.split,
                                                          where.split  = data.bin.mixed$where.split,
                                                          dir.split    = data.bin.mixed$dir.split,
                                                          split.cate   = data.bin.mixed$split.cate)
eval.final.estg.glm.modinnd.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("32")

#####################################################################################################################
############################ 33. g: GLM Model, inside node, Mis func adjustment model, Cv2 ##########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.nois.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                         est.used          = "G",
                                                         type.var          = "bin",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = NULL,
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 15,
                                                         min.node          = 15)

final.tree.estg.glm.modinnd.nois.cv2 <- EstG.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                       tree.list         = seq.created.estg.glm.modinnd.nois.cv2$tree.list,             
                                                       type.var          = "bin",
                                                       seed              = a[job.number],
                                                       n.cv              = 5,
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = NULL,
                                                       min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.nois.cv2[[1]],
                                                          test.data    = data.bin.mixed$test.data,
                                                          true.trt.eff = data.bin.mixed$true.trt.eff,
                                                          noise.var    = data.bin.mixed$noise.var,
                                                          corr.split   = data.bin.mixed$corr.split,
                                                          where.split  = data.bin.mixed$where.split,
                                                          dir.split    = data.bin.mixed$dir.split,
                                                          split.cate   = data.bin.mixed$split.cate)
eval.final.estg.glm.modinnd.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("33")

#####################################################################################################################
######################### 34. g: GLM Model, inside node, Unmeasured cov adjustment model, Cv2 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.mis.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed.mis,
                                                        est.used          = "G",
                                                        type.var          = "bin",
                                                        adj.mod.out       = F,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                        w                 = NULL,
                                                        adj.mod.insplt    = F,
                                                        num.truc.obs      = 15,
                                                        min.node          = 15)

final.tree.estg.glm.modinnd.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.bin.mixed.mis,
                                                      tree.list         = seq.created.estg.glm.modinnd.mis.cv2$tree.list,             
                                                      type.var          = "bin",
                                                      seed              = a[job.number],
                                                      n.cv              = 5,
                                                      adj.mod.out       = F,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                      min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.mis.cv2[[1]],
                                                         test.data    = test.data.mis,
                                                         true.trt.eff = data.bin.mixed$true.trt.eff,
                                                         noise.var    = data.bin.mixed$noise.var,
                                                         corr.split   = data.bin.mixed$corr.split,
                                                         where.split  = data.bin.mixed$where.split,
                                                         dir.split    = data.bin.mixed$dir.split,
                                                         split.cate   = data.bin.mixed$split.cate)
eval.final.estg.glm.modinnd.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("34")

performance.homo.g <- list(glm.modinnd.true.cv1 = eval.final.estg.glm.modinnd.true.cv1,
                           glm.modinnd.nois.cv1 = eval.final.estg.glm.modinnd.nois.cv1,
                           glm.modinnd.mis.cv1  = eval.final.estg.glm.modinnd.mis.cv1,
                           glm.modinnd.true.cv2 = eval.final.estg.glm.modinnd.true.cv2,
                           glm.modinnd.nois.cv2 = eval.final.estg.glm.modinnd.nois.cv2,
                           glm.modinnd.mis.cv2  = eval.final.estg.glm.modinnd.mis.cv2)

#####################################################################################################################
###################### 35. dr: True GLM Model in node, True propensity score model in node, Cv1 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                    est.used          = "DR",
                                                                    type.var          = "bin",
                                                                    propsc.mod.out    = F,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                    propsc.mod.insplt = F,
                                                                    adj.mod.out       = F, 
                                                                    adj.mthd          = "GLM", 
                                                                    adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))", 
                                                                    adj.mod.insplt    = F, 
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)

final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                   tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                   lambda.used       = qchisq(0.95, 1),
                                                                   val.sample        = data.validation.bin.mixed,
                                                                   type.var          = "bin",
                                                                   propsc.mod.out    = F,
                                                                   propsc.mthd       = "GLM",
                                                                   propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                   propsc.mod.insplt = F,
                                                                   adj.mod.out       = F,
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))", 
                                                                   adj.mod.insplt    = F,
                                                                   min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1[[1]],
                                                                     test.data    = data.bin.mixed$test.data,
                                                                     true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                     noise.var    = data.bin.mixed$noise.var,
                                                                     corr.split   = data.bin.mixed$corr.split,
                                                                     where.split  = data.bin.mixed$where.split,
                                                                     dir.split    = data.bin.mixed$dir.split,
                                                                     split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("35")

#####################################################################################################################
################### 36. dr: True GLM Model in node, Mis func propensity score model in node, Cv1 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                       est.used          = "DR",
                                                                       type.var          = "bin",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                       propsc.mod.insplt = F,
                                                                       adj.mod.out       = F, 
                                                                       adj.mthd          = "GLM", 
                                                                       adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))", 
                                                                       adj.mod.insplt    = F, 
                                                                       num.truc.obs      = 30,
                                                                       min.node          = 20)

final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                      tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                      lambda.used       = qchisq(0.95, 1),
                                                                      val.sample        = data.validation.bin.mixed,
                                                                      type.var          = "bin",
                                                                      propsc.mod.out    = F,
                                                                      propsc.mthd       = "GLM",
                                                                      propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                      propsc.mod.insplt = F,
                                                                      adj.mod.out       = F,
                                                                      adj.mthd          = "GLM",
                                                                      adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))", 
                                                                      adj.mod.insplt    = F,
                                                                      min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1[[1]],
                                                                        test.data    = data.bin.mixed$test.data,
                                                                        true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                        noise.var    = data.bin.mixed$noise.var,
                                                                        corr.split   = data.bin.mixed$corr.split,
                                                                        where.split  = data.bin.mixed$where.split,
                                                                        dir.split    = data.bin.mixed$dir.split,
                                                                        split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("36")

#####################################################################################################################
##################### 37. dr: Mis func GLM Model in node, True propensity score model in node, Cv1 ##################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                       est.used          = "DR",
                                                                       type.var          = "bin",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                       propsc.mod.insplt = F,
                                                                       adj.mod.out       = F, 
                                                                       adj.mthd          = "GLM", 
                                                                       adj.form.true     = NULL, 
                                                                       adj.mod.insplt    = F, 
                                                                       num.truc.obs      = 30,
                                                                       min.node          = 20)

final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                      tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                      lambda.used       = qchisq(0.95, 1),
                                                                      val.sample        = data.validation.bin.mixed,
                                                                      type.var          = "bin",
                                                                      propsc.mod.out    = F,
                                                                      propsc.mthd       = "GLM",
                                                                      propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                      propsc.mod.insplt = F,
                                                                      adj.mod.out       = F,
                                                                      adj.mthd          = "GLM",
                                                                      adj.form.true     = NULL, 
                                                                      adj.mod.insplt    = F,
                                                                      min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1[[1]],
                                                                        test.data    = data.bin.mixed$test.data,
                                                                        true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                        noise.var    = data.bin.mixed$noise.var,
                                                                        corr.split   = data.bin.mixed$corr.split,
                                                                        where.split  = data.bin.mixed$where.split,
                                                                        dir.split    = data.bin.mixed$dir.split,
                                                                        split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("37")

#####################################################################################################################
################## 38. dr: Mis func GLM Model in node, Mis func propensity score model in node, Cv1 #################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                          est.used          = "DR",
                                                                          type.var          = "bin",
                                                                          propsc.mod.out    = F,
                                                                          propsc.mthd       = "GLM",
                                                                          propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                          propsc.mod.insplt = F,
                                                                          adj.mod.out       = F, 
                                                                          adj.mthd          = "GLM", 
                                                                          adj.form.true     = NULL, 
                                                                          adj.mod.insplt    = F, 
                                                                          num.truc.obs      = 30,
                                                                          min.node          = 20)

final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                         tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                         lambda.used       = qchisq(0.95, 1),
                                                                         val.sample        = data.validation.bin.mixed,
                                                                         type.var          = "bin",
                                                                         propsc.mod.out    = F,
                                                                         propsc.mthd       = "GLM",
                                                                         propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                         propsc.mod.insplt = F,
                                                                         adj.mod.out       = F,
                                                                         adj.mthd          = "GLM",
                                                                         adj.form.true     = NULL, 
                                                                         adj.mod.insplt    = F,
                                                                         min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1[[1]],
                                                                           test.data    = data.bin.mixed$test.data,
                                                                           true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                           noise.var    = data.bin.mixed$noise.var,
                                                                           corr.split   = data.bin.mixed$corr.split,
                                                                           where.split  = data.bin.mixed$where.split,
                                                                           dir.split    = data.bin.mixed$dir.split,
                                                                           split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("38")

#####################################################################################################################
############# 39. dr: Unmeasured cov GLM Model in node, Unmeasured cov propensity score model in node, Cv1 ##########
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed.mis,
                                                                    est.used          = "DR",
                                                                    type.var          = "bin",
                                                                    propsc.mod.out    = F,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                    propsc.mod.insplt = F,
                                                                    adj.mod.out       = F, 
                                                                    adj.mthd          = "GLM", 
                                                                    adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                    adj.mod.insplt    = F, 
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)

final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed.mis,
                                                                   tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1$tree.list, 
                                                                   lambda.used       = qchisq(0.95, 1),
                                                                   val.sample        = data.validation.bin.mixed.mis,
                                                                   type.var          = "bin",
                                                                   propsc.mod.out    = F,
                                                                   propsc.mthd       = "GLM",
                                                                   propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                   propsc.mod.insplt = F,
                                                                   adj.mod.out       = F,
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                   adj.mod.insplt    = F,
                                                                   min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1[[1]],
                                                                     test.data    = test.data.mis,
                                                                     true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                     noise.var    = data.bin.mixed$noise.var,
                                                                     corr.split   = data.bin.mixed$corr.split,
                                                                     where.split  = data.bin.mixed$where.split,
                                                                     dir.split    = data.bin.mixed$dir.split,
                                                                     split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("39")

#####################################################################################################################
###################### 40. dr: True GLM Model in node, True propensity score model in node, Cv2 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                    est.used          = "DR",
                                                                    type.var          = "bin",
                                                                    propsc.mod.out    = F,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                    propsc.mod.insplt = F,
                                                                    adj.mod.out       = F, 
                                                                    adj.mthd          = "GLM", 
                                                                    adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))", 
                                                                    adj.mod.insplt    = F, 
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)

final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                   tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2$tree.list,
                                                                   type.var          = "bin",
                                                                   seed              = a[job.number], 
                                                                   n.cv              = 5,
                                                                   propsc.mod.out    = F,
                                                                   propsc.mthd       = "GLM",
                                                                   propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                   min.obs.propsc    = 5,
                                                                   adj.mod.out       = F,
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))", 
                                                                   min.obs.adj       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2[[1]],
                                                                     test.data    = data.bin.mixed$test.data,
                                                                     true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                     noise.var    = data.bin.mixed$noise.var,
                                                                     corr.split   = data.bin.mixed$corr.split,
                                                                     where.split  = data.bin.mixed$where.split,
                                                                     dir.split    = data.bin.mixed$dir.split,
                                                                     split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("40")

#####################################################################################################################
################### 41. dr: True GLM Model in node, Mis func propensity score model in node, Cv2 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                       est.used          = "DR",
                                                                       type.var          = "bin",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                       propsc.mod.insplt = F,
                                                                       adj.mod.out       = F, 
                                                                       adj.mthd          = "GLM", 
                                                                       adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))", 
                                                                       adj.mod.insplt    = F, 
                                                                       num.truc.obs      = 30,
                                                                       min.node          = 20)

final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                      tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$tree.list,
                                                                      type.var          = "bin",
                                                                      seed              = a[job.number], 
                                                                      n.cv              = 5,
                                                                      propsc.mod.out    = F,
                                                                      propsc.mthd       = "GLM",
                                                                      propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                      min.obs.propsc    = 10,
                                                                      adj.mod.out       = F,
                                                                      adj.mthd          = "GLM",
                                                                      adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))", 
                                                                      min.obs.adj       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2[[1]],
                                                                        test.data    = data.bin.mixed$test.data,
                                                                        true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                        noise.var    = data.bin.mixed$noise.var,
                                                                        corr.split   = data.bin.mixed$corr.split,
                                                                        where.split  = data.bin.mixed$where.split,
                                                                        dir.split    = data.bin.mixed$dir.split,
                                                                        split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("41")

#####################################################################################################################
#################### 42. dr: Mis func GLM Model in node, True propensity score model in node, Cv2 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                       est.used          = "DR",
                                                                       type.var          = "bin",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                       propsc.mod.insplt = F,
                                                                       adj.mod.out       = F, 
                                                                       adj.mthd          = "GLM", 
                                                                       adj.form.true     = NULL, 
                                                                       adj.mod.insplt    = F, 
                                                                       num.truc.obs      = 30,
                                                                       min.node          = 20)

final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                      tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$tree.list,
                                                                      type.var          = "bin",
                                                                      seed              = a[job.number], 
                                                                      n.cv              = 5,
                                                                      propsc.mod.out    = F,
                                                                      propsc.mthd       = "GLM",
                                                                      propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                      min.obs.propsc    = 5,
                                                                      adj.mod.out       = F,
                                                                      adj.mthd          = "GLM",
                                                                      adj.form.true     = NULL, 
                                                                      min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2[[1]],
                                                                        test.data    = data.bin.mixed$test.data,
                                                                        true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                        noise.var    = data.bin.mixed$noise.var,
                                                                        corr.split   = data.bin.mixed$corr.split,
                                                                        where.split  = data.bin.mixed$where.split,
                                                                        dir.split    = data.bin.mixed$dir.split,
                                                                        split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("42")

#####################################################################################################################
################## 43. dr: Mis func GLM Model in node, Mis func propensity score model in node, Cv2 #################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                          est.used          = "DR",
                                                                          type.var          = "bin",
                                                                          propsc.mod.out    = F,
                                                                          propsc.mthd       = "GLM",
                                                                          propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                          propsc.mod.insplt = F,
                                                                          adj.mod.out       = F, 
                                                                          adj.mthd          = "GLM", 
                                                                          adj.form.true     = NULL, 
                                                                          adj.mod.insplt    = F, 
                                                                          num.truc.obs      = 30,
                                                                          min.node          = 20)

final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                         tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$tree.list,
                                                                         type.var          = "bin",
                                                                         seed              = a[job.number], 
                                                                         n.cv              = 5,
                                                                         propsc.mod.out    = F,
                                                                         propsc.mthd       = "GLM",
                                                                         propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                         min.obs.propsc    = 10,
                                                                         adj.mod.out       = F,
                                                                         adj.mthd          = "GLM",
                                                                         adj.form.true     = NULL, 
                                                                         min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2[[1]],
                                                                           test.data    = data.bin.mixed$test.data,
                                                                           true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                           noise.var    = data.bin.mixed$noise.var,
                                                                           corr.split   = data.bin.mixed$corr.split,
                                                                           where.split  = data.bin.mixed$where.split,
                                                                           dir.split    = data.bin.mixed$dir.split,
                                                                           split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("43")

#####################################################################################################################
############ 44. dr: Unmeasured cov GLM Model in node, Unmeasured cov propensity score model in node, Cv2 ###########
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed.mis,
                                                                    est.used          = "DR",
                                                                    type.var          = "bin",
                                                                    propsc.mod.out    = F,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                    propsc.mod.insplt = F,
                                                                    adj.mod.out       = F, 
                                                                    adj.mthd          = "GLM", 
                                                                    adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                    adj.mod.insplt    = F, 
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)

final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed.mis,
                                                                   tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2$tree.list,
                                                                   type.var          = "bin",
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

eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv2[[1]],
                                                                     test.data    = test.data.mis,
                                                                     true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                     noise.var    = data.bin.mixed$noise.var,
                                                                     corr.split   = data.bin.mixed$corr.split,
                                                                     where.split  = data.bin.mixed$where.split,
                                                                     dir.split    = data.bin.mixed$dir.split,
                                                                     split.cate   = data.bin.mixed$split.cate)
eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("44")

performance.homo.drInnd <- list(adjTGlmInnd.propscTGlmInnd.cv1       = eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1,
                                adjTGlmInnd.propscNoisGlmInnd.cv1    = eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1,
                                adjNoisGlmInnd.propscTGlmInnd.cv1    = eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1,
                                adjNoisGlmInnd.propscNoisGlmInnd.cv1 = eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1,
                                adjFGlmInnd.propscFGlmInnd.cv1       = eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1,
                                adjTGlmInnd.propscTGlmInnd.cv2       = eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2,
                                adjTGlmInnd.propscNoisGlmInnd.cv2    = eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2,
                                adjNoisGlmInnd.propscTGlmInnd.cv2    = eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2,
                                adjNoisGlmInnd.propscNoisGlmInnd.cv2 = eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2,
                                adjFGlmInnd.propscFGlmInnd.cv2       = eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)

