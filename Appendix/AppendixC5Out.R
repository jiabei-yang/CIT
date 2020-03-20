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
######################### 1. ipw: GLM Model, outside node, True propensity score model, cv1 #########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscout.true.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                             est.used          = "IPW",
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T,
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = "A ~ X1 + X2 + X3",
                                                             w                 = NULL,
                                                             propsc.mod.insplt = NULL,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscout.true.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                             tree.list         = seq.created.estipw.glm.propscout.true.cv1$tree.list, 
                                                             lambda.used       = qchisq(0.95, 1), 
                                                             val.sample        = data.validation.cont.cont, 
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T, 
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = "A ~ X1 + X2 + X3",
                                                             val.w             = NULL,
                                                             propsc.mod.insplt = NULL, 
                                                             min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estipw.glm.propscout.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscout.true.cv1[[1]], 
                                                              test.data    = data.cont.cont$test.data,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split,
                                                              where.split  = data.cont.cont$where.split,
                                                              dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscout.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscout.true.cv1$corr.frst.splt <- as.character(seq.created.estipw.glm.propscout.true.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("1")

#####################################################################################################################
######################## 2. ipw: GLM Model, outside node, Mis func propensity score model, cv1 ######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscout.nois.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                             est.used          = "IPW",
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T,
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)", 
                                                             w                 = NULL, 
                                                             propsc.mod.insplt = NULL,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscout.nois.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                             tree.list         = seq.created.estipw.glm.propscout.nois.cv1$tree.list, 
                                                             lambda.used       = qchisq(0.95, 1), 
                                                             val.sample        = data.validation.cont.cont, 
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T, 
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)", 
                                                             val.w             = NULL,
                                                             propsc.mod.insplt = NULL, 
                                                             min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estipw.glm.propscout.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscout.nois.cv1[[1]], 
                                                              test.data    = data.cont.cont$test.data,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split,
                                                              where.split  = data.cont.cont$where.split,
                                                              dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscout.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscout.nois.cv1$corr.frst.splt <- as.character(seq.created.estipw.glm.propscout.nois.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("2")

#####################################################################################################################
##################### 3. ipw: GLM Model, outside node, Unmeasured cov propensity score model, cv1 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscout.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                            est.used          = "IPW",
                                                            type.var          = "cont",
                                                            propsc.mod.out    = T,
                                                            propsc.mthd       = "GLM", 
                                                            propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6", 
                                                            w                 = NULL, 
                                                            propsc.mod.insplt = NULL,
                                                            num.truc.obs      = 30,
                                                            min.node          = 20)

final.tree.estipw.glm.propscout.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont.mis, 
                                                            tree.list         = seq.created.estipw.glm.propscout.mis.cv1$tree.list, 
                                                            lambda.used       = qchisq(0.95, 1), 
                                                            val.sample        = data.validation.cont.cont.mis, 
                                                            type.var          = "cont",
                                                            propsc.mod.out    = T, 
                                                            propsc.mthd       = "GLM", 
                                                            propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6", 
                                                            val.w             = NULL,
                                                            propsc.mod.insplt = NULL, 
                                                            min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estipw.glm.propscout.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscout.mis.cv1[[1]], 
                                                             test.data    = test.data.mis,
                                                             true.trt.eff = data.cont.cont$true.trt.eff,
                                                             noise.var    = data.cont.cont$noise.var,
                                                             corr.split   = data.cont.cont$corr.split,
                                                             where.split  = data.cont.cont$where.split,
                                                             dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscout.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscout.mis.cv1$corr.frst.splt <- as.character(seq.created.estipw.glm.propscout.mis.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("3")

#####################################################################################################################
######################### 4. ipw: GLM Model, outside node, True propensity score model, Cv2 #########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscout.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                             est.used          = "IPW",
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T,
                                                             propsc.mthd       = "GLM",
                                                             propsc.form.true  = "A ~ X1 + X2 + X3",
                                                             w                 = NULL,
                                                             propsc.mod.insplt = NULL,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscout.true.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont,
                                                             tree.list        = seq.created.estipw.glm.propscout.true.cv2$tree.list,
                                                             type.var         = "cont",
                                                             seed             = a[job.number],
                                                             n.cv             = 5,
                                                             propsc.mod.out   = T, 
                                                             propsc.mthd      = "GLM", 
                                                             propsc.form.true = "A ~ X1 + X2 + X3",
                                                             min.obs.mod      = NULL)
t1 <- Sys.time()

eval.final.estipw.glm.propscout.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscout.true.cv2[[1]],
                                                              test.data    = data.cont.cont$test.data,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split,
                                                              where.split  = data.cont.cont$where.split,
                                                              dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscout.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscout.true.cv2$corr.frst.splt <- as.character(seq.created.estipw.glm.propscout.true.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("4")

#####################################################################################################################
####################### 5. ipw: GLM Model, outside node, Mis func propensity score model, Cv2 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscout.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                             est.used          = "IPW",
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T,
                                                             propsc.mthd       = "GLM",
                                                             propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                             w                 = NULL,
                                                             propsc.mod.insplt = NULL,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscout.nois.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont,
                                                             tree.list        = seq.created.estipw.glm.propscout.nois.cv2$tree.list,
                                                             type.var         = "cont",
                                                             seed             = a[job.number],
                                                             n.cv             = 5,
                                                             propsc.mod.out   = T, 
                                                             propsc.mthd      = "GLM", 
                                                             propsc.form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                             min.obs.mod      = NULL)
t1 <- Sys.time()

eval.final.estipw.glm.propscout.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscout.nois.cv2[[1]],
                                                              test.data    = data.cont.cont$test.data,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split,
                                                              where.split  = data.cont.cont$where.split,
                                                              dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscout.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscout.nois.cv2$corr.frst.splt <- as.character(seq.created.estipw.glm.propscout.nois.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("5")

#####################################################################################################################
###################### 6. ipw: GLM Model, outside node, Unmeasured cov propensity score model, Cv2 ##################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscout.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                            est.used          = "IPW",
                                                            type.var          = "cont",
                                                            propsc.mod.out    = T,
                                                            propsc.mthd       = "GLM",
                                                            propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                            w                 = NULL,
                                                            propsc.mod.insplt = NULL,
                                                            num.truc.obs      = 30,
                                                            min.node          = 20)

final.tree.estipw.glm.propscout.mis.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont.mis,
                                                            tree.list        = seq.created.estipw.glm.propscout.mis.cv2$tree.list,
                                                            type.var         = "cont",
                                                            seed             = a[job.number],
                                                            n.cv             = 5,
                                                            propsc.mod.out   = T, 
                                                            propsc.mthd      = "GLM", 
                                                            propsc.form.true = "A ~ X1 + X3 + X4 + X5 + X6",
                                                            min.obs.mod      = NULL)
t1 <- Sys.time()

eval.final.estipw.glm.propscout.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscout.mis.cv2[[1]],
                                                             test.data    = test.data.mis,
                                                             true.trt.eff = data.cont.cont$true.trt.eff,
                                                             noise.var    = data.cont.cont$noise.var,
                                                             corr.split   = data.cont.cont$corr.split,
                                                             where.split  = data.cont.cont$where.split,
                                                             dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscout.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscout.mis.cv2$corr.frst.splt <- as.character(seq.created.estipw.glm.propscout.mis.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("6")

performance.hetero.ipw <- list(glm.propscout.true.cv1 = eval.final.estipw.glm.propscout.true.cv1,
                               glm.propscout.nois.cv1 = eval.final.estipw.glm.propscout.nois.cv1,
                               glm.propscout.mis.cv1  = eval.final.estipw.glm.propscout.mis.cv1,
                               glm.propscout.true.cv2 = eval.final.estipw.glm.propscout.true.cv2,
                               glm.propscout.nois.cv2 = eval.final.estipw.glm.propscout.nois.cv2,
                               glm.propscout.mis.cv2  = eval.final.estipw.glm.propscout.mis.cv2)

#####################################################################################################################
############################# 7. g: GLM Model, outside node, True adjustment model, cv1 #############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modout.true.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                        est.used          = "G",
                                                        type.var          = "cont",
                                                        adj.mod.out       = T,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                        w                 = NULL,
                                                        adj.mod.insplt    = NULL,
                                                        num.truc.obs      = 10,
                                                        min.node          = 10)

final.tree.estg.glm.modout.true.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                      tree.list         = seq.created.estg.glm.modout.true.cv1$tree.list,
                                                      lambda.used       = qchisq(0.95, 1),
                                                      val.sample        = data.validation.cont.cont,              
                                                      type.var          = "cont",
                                                      adj.mod.out       = T,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                      val.w             = NULL,
                                                      adj.mod.insplt    = NULL,
                                                      min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estg.glm.modout.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.true.cv1[[1]],
                                                         test.data    = data.cont.cont$test.data,
                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                         noise.var    = data.cont.cont$noise.var,
                                                         corr.split   = data.cont.cont$corr.split,
                                                         where.split  = data.cont.cont$where.split,
                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modout.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modout.true.cv1$corr.frst.splt <- as.character(seq.created.estg.glm.modout.true.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("7")

#####################################################################################################################
########################### 8. g: GLM Model, outside node, Mis func adjustment model, cv1 ###########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modout.nois.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                        est.used          = "G",
                                                        type.var          = "cont",
                                                        adj.mod.out       = T,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = NULL,
                                                        w                 = NULL,
                                                        adj.mod.insplt    = NULL,
                                                        num.truc.obs      = 15,
                                                        min.node          = 15)

final.tree.estg.glm.modout.nois.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                      tree.list         = seq.created.estg.glm.modout.nois.cv1$tree.list,
                                                      lambda.used       = qchisq(0.95, 1),
                                                      val.sample        = data.validation.cont.cont,              
                                                      type.var          = "cont",
                                                      adj.mod.out       = T,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = NULL,
                                                      val.w             = NULL,
                                                      adj.mod.insplt    = NULL,
                                                      min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estg.glm.modout.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.nois.cv1[[1]],
                                                         test.data    = data.cont.cont$test.data,
                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                         noise.var    = data.cont.cont$noise.var,
                                                         corr.split   = data.cont.cont$corr.split,
                                                         where.split  = data.cont.cont$where.split,
                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modout.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modout.nois.cv1$corr.frst.splt <- as.character(seq.created.estg.glm.modout.nois.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("8")

#####################################################################################################################
######################## 9. g: GLM Model, outside node, Unmeasured cov adjustment model, cv1 ########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modout.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                       est.used          = "G",
                                                       type.var          = "cont",
                                                       adj.mod.out       = T,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                       w                 = NULL,
                                                       adj.mod.insplt    = NULL,
                                                       num.truc.obs      = 15,
                                                       min.node          = 15)

final.tree.estg.glm.modout.mis.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                     tree.list         = seq.created.estg.glm.modout.mis.cv1$tree.list,
                                                     lambda.used       = qchisq(0.95, 1),
                                                     val.sample        = data.validation.cont.cont.mis,              
                                                     type.var          = "cont",
                                                     adj.mod.out       = T,
                                                     adj.mthd          = "GLM",
                                                     adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                     val.w             = NULL,
                                                     adj.mod.insplt    = NULL,
                                                     min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estg.glm.modout.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.mis.cv1[[1]],
                                                        test.data    = test.data.mis,
                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                        noise.var    = data.cont.cont$noise.var,
                                                        corr.split   = data.cont.cont$corr.split,
                                                        where.split  = data.cont.cont$where.split,
                                                        dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modout.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modout.mis.cv1$corr.frst.splt <- as.character(seq.created.estg.glm.modout.mis.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("9")

#####################################################################################################################
########################### 10. g: GLM Model, outside node, True adjustment model, Cv2 ##############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modout.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                        est.used          = "G",
                                                        type.var          = "cont",
                                                        adj.mod.out       = T,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                        w                 = NULL,
                                                        adj.mod.insplt    = NULL,
                                                        num.truc.obs      = 10,
                                                        min.node          = 10)

final.tree.estg.glm.modout.true.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                      tree.list         = seq.created.estg.glm.modout.true.cv2$tree.list,             
                                                      type.var          = "cont",
                                                      seed              = a[job.number],
                                                      n.cv              = 5,
                                                      adj.mod.out       = T,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                      min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estg.glm.modout.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.true.cv2[[1]],
                                                         test.data    = data.cont.cont$test.data,
                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                         noise.var    = data.cont.cont$noise.var,
                                                         corr.split   = data.cont.cont$corr.split,
                                                         where.split  = data.cont.cont$where.split,
                                                         dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modout.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modout.true.cv2$corr.frst.splt <- as.character(seq.created.estg.glm.modout.true.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("10")

#####################################################################################################################
########################## 11. g: GLM Model, outside node, Mis func adjustment model, Cv2 ###########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modout.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                        est.used          = "G",
                                                        type.var          = "cont",
                                                        adj.mod.out       = T,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = NULL,
                                                        w                 = NULL,
                                                        adj.mod.insplt    = NULL,
                                                        num.truc.obs      = 15,
                                                        min.node          = 15)

final.tree.estg.glm.modout.nois.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                      tree.list         = seq.created.estg.glm.modout.nois.cv2$tree.list,             
                                                      type.var          = "cont",
                                                      seed              = a[job.number],
                                                      n.cv              = 5,
                                                      adj.mod.out       = T,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = NULL,
                                                      min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estg.glm.modout.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.nois.cv2[[1]],
                                                         test.data    = data.cont.cont$test.data,
                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                         noise.var    = data.cont.cont$noise.var,
                                                         corr.split   = data.cont.cont$corr.split,
                                                         where.split  = data.cont.cont$where.split,
                                                         dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modout.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modout.nois.cv2$corr.frst.splt <- as.character(seq.created.estg.glm.modout.nois.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("11")

#####################################################################################################################
####################### 12. g: GLM Model, outside node, Unmeasured cov adjustment model, Cv2 ########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modout.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                       est.used          = "G",
                                                       type.var          = "cont",
                                                       adj.mod.out       = T,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                       w                 = NULL,
                                                       adj.mod.insplt    = NULL,
                                                       num.truc.obs      = 15,
                                                       min.node          = 15)

final.tree.estg.glm.modout.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                     tree.list         = seq.created.estg.glm.modout.mis.cv2$tree.list,             
                                                     type.var          = "cont",
                                                     seed              = a[job.number],
                                                     n.cv              = 5,
                                                     adj.mod.out       = T,
                                                     adj.mthd          = "GLM",
                                                     adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                     min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estg.glm.modout.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.mis.cv2[[1]],
                                                        test.data    = test.data.mis,
                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                        noise.var    = data.cont.cont$noise.var,
                                                        corr.split   = data.cont.cont$corr.split,
                                                        where.split  = data.cont.cont$where.split,
                                                        dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modout.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modout.mis.cv2$corr.frst.splt <- as.character(seq.created.estg.glm.modout.mis.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("12")

performance.hetero.g <- list(glm.modout.true.cv1 = eval.final.estg.glm.modout.true.cv1,
                             glm.modout.nois.cv1 = eval.final.estg.glm.modout.nois.cv1,
                             glm.modout.mis.cv1  = eval.final.estg.glm.modout.mis.cv1,
                             glm.modout.true.cv2 = eval.final.estg.glm.modout.true.cv2,
                             glm.modout.nois.cv2 = eval.final.estg.glm.modout.nois.cv2,
                             glm.modout.mis.cv2  = eval.final.estg.glm.modout.mis.cv2)

#####################################################################################################################
###################### 13. dr: True GLM Model outside, True propensity score model outside, Cv1 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmOut.propscTGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                  est.used          = "DR",
                                                                  type.var          = "cont",
                                                                  propsc.mod.out    = T,
                                                                  propsc.mthd       = "GLM",
                                                                  propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                  propsc.mod.insplt = NULL,
                                                                  adj.mod.out       = T, 
                                                                  adj.mthd          = "GLM", 
                                                                  adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                  adj.mod.insplt    = NULL, 
                                                                  num.truc.obs      = 30,
                                                                  min.node          = 20)

final.tree.estdr.adjTGlmOut.propscTGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                 tree.list         = seq.created.estdr.adjTGlmOut.propscTGlmOut.cv1$tree.list, 
                                                                 lambda.used       = qchisq(0.95, 1),
                                                                 val.sample        = data.validation.cont.cont,
                                                                 type.var          = "cont",
                                                                 propsc.mod.out    = T,
                                                                 propsc.mthd       = "GLM",
                                                                 propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                 propsc.mod.insplt = NULL,
                                                                 adj.mod.out       = T,
                                                                 adj.mthd          = "GLM",
                                                                 adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                 adj.mod.insplt    = NULL,
                                                                 min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmOut.propscTGlmOut.cv1[[1]],
                                                                   test.data    = data.cont.cont$test.data,
                                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                                   noise.var    = data.cont.cont$noise.var,
                                                                   corr.split   = data.cont.cont$corr.split,
                                                                   where.split  = data.cont.cont$where.split,
                                                                   dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjTGlmOut.propscTGlmOut.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("13")

#####################################################################################################################
#################### 14. dr: True GLM Model outside, Mis func propensity score model outside, Cv1 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                     est.used          = "DR",
                                                                     type.var          = "cont",
                                                                     propsc.mod.out    = T,
                                                                     propsc.mthd       = "GLM",
                                                                     propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                     propsc.mod.insplt = NULL,
                                                                     adj.mod.out       = T, 
                                                                     adj.mthd          = "GLM", 
                                                                     adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                     adj.mod.insplt    = NULL, 
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)

final.tree.estdr.adjTGlmOut.propscNoisGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                    tree.list         = seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv1$tree.list, 
                                                                    lambda.used       = qchisq(0.95, 1),
                                                                    val.sample        = data.validation.cont.cont,
                                                                    type.var          = "cont",
                                                                    propsc.mod.out    = T,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                    propsc.mod.insplt = NULL,
                                                                    adj.mod.out       = T,
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                    adj.mod.insplt    = NULL,
                                                                    min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmOut.propscNoisGlmOut.cv1[[1]],
                                                                      test.data    = data.cont.cont$test.data,
                                                                      true.trt.eff = data.cont.cont$true.trt.eff,
                                                                      noise.var    = data.cont.cont$noise.var,
                                                                      corr.split   = data.cont.cont$corr.split,
                                                                      where.split  = data.cont.cont$where.split,
                                                                      dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("14")

#####################################################################################################################
#################### 15. dr: Mis func GLM Model outside, True propensity score model outside, Cv1 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                     est.used          = "DR",
                                                                     type.var          = "cont",
                                                                     propsc.mod.out    = T,
                                                                     propsc.mthd       = "GLM",
                                                                     propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                     propsc.mod.insplt = NULL,
                                                                     adj.mod.out       = T, 
                                                                     adj.mthd          = "GLM", 
                                                                     adj.form.true     = NULL, 
                                                                     adj.mod.insplt    = NULL, 
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)

final.tree.estdr.adjNoisGlmOut.propscTGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                    tree.list         = seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv1$tree.list, 
                                                                    lambda.used       = qchisq(0.95, 1),
                                                                    val.sample        = data.validation.cont.cont,
                                                                    type.var          = "cont",
                                                                    propsc.mod.out    = T,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                    propsc.mod.insplt = NULL,
                                                                    adj.mod.out       = T,
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = NULL, 
                                                                    adj.mod.insplt    = NULL,
                                                                    min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmOut.propscTGlmOut.cv1[[1]],
                                                                      test.data    = data.cont.cont$test.data,
                                                                      true.trt.eff = data.cont.cont$true.trt.eff,
                                                                      noise.var    = data.cont.cont$noise.var,
                                                                      corr.split   = data.cont.cont$corr.split,
                                                                      where.split  = data.cont.cont$where.split,
                                                                      dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("15")

#####################################################################################################################
################## 16. dr: Mis func GLM Model outside, Mis func propensity score model outside, Cv1 #################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                        est.used          = "DR",
                                                                        type.var          = "cont",
                                                                        propsc.mod.out    = T,
                                                                        propsc.mthd       = "GLM",
                                                                        propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                        propsc.mod.insplt = NULL,
                                                                        adj.mod.out       = T, 
                                                                        adj.mthd          = "GLM", 
                                                                        adj.form.true     = NULL, 
                                                                        adj.mod.insplt    = NULL, 
                                                                        num.truc.obs      = 30,
                                                                        min.node          = 20)

final.tree.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                       tree.list         = seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1$tree.list, 
                                                                       lambda.used       = qchisq(0.95, 1),
                                                                       val.sample        = data.validation.cont.cont,
                                                                       type.var          = "cont",
                                                                       propsc.mod.out    = T,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                       propsc.mod.insplt = NULL,
                                                                       adj.mod.out       = T,
                                                                       adj.mthd          = "GLM",
                                                                       adj.form.true     = NULL, 
                                                                       adj.mod.insplt    = NULL,
                                                                       min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1[[1]],
                                                                         test.data    = data.cont.cont$test.data,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("16")

#####################################################################################################################
############ 17. dr: Unmeasured cov GLM Model outside, Unmeasured cov propensity score model outside, Cv1 ###########
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmOut.propscFGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                                  est.used          = "DR",
                                                                  type.var          = "cont",
                                                                  propsc.mod.out    = T,
                                                                  propsc.mthd       = "GLM",
                                                                  propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                  propsc.mod.insplt = NULL,
                                                                  adj.mod.out       = T, 
                                                                  adj.mthd          = "GLM", 
                                                                  adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                  adj.mod.insplt    = NULL, 
                                                                  num.truc.obs      = 30,
                                                                  min.node          = 20)

final.tree.estdr.adjFGlmOut.propscFGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                                 tree.list         = seq.created.estdr.adjFGlmOut.propscFGlmOut.cv1$tree.list, 
                                                                 lambda.used       = qchisq(0.95, 1),
                                                                 val.sample        = data.validation.cont.cont.mis,
                                                                 type.var          = "cont",
                                                                 propsc.mod.out    = T,
                                                                 propsc.mthd       = "GLM",
                                                                 propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                 propsc.mod.insplt = NULL,
                                                                 adj.mod.out       = T,
                                                                 adj.mthd          = "GLM",
                                                                 adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                 adj.mod.insplt    = NULL,
                                                                 min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmOut.propscFGlmOut.cv1[[1]],
                                                                   test.data    = test.data.mis,
                                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                                   noise.var    = data.cont.cont$noise.var,
                                                                   corr.split   = data.cont.cont$corr.split,
                                                                   where.split  = data.cont.cont$where.split,
                                                                   dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjFGlmOut.propscFGlmOut.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("17")

#####################################################################################################################
###################### 18. dr: True GLM Model outside, True propensity score model outside, Cv2 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmOut.propscTGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                  est.used          = "DR",
                                                                  type.var          = "cont",
                                                                  propsc.mod.out    = T,
                                                                  propsc.mthd       = "GLM",
                                                                  propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                  propsc.mod.insplt = NULL,
                                                                  adj.mod.out       = T, 
                                                                  adj.mthd          = "GLM", 
                                                                  adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                  adj.mod.insplt    = NULL, 
                                                                  num.truc.obs      = 30,
                                                                  min.node          = 20)

final.tree.estdr.adjTGlmOut.propscTGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                 tree.list         = seq.created.estdr.adjTGlmOut.propscTGlmOut.cv2$tree.list,
                                                                 type.var          = "cont",
                                                                 seed              = a[job.number], 
                                                                 n.cv              = 5,
                                                                 propsc.mod.out    = T,
                                                                 propsc.mthd       = "GLM",
                                                                 propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                 min.obs.propsc    = 5,
                                                                 adj.mod.out       = T,
                                                                 adj.mthd          = "GLM",
                                                                 adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                 min.obs.adj       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmOut.propscTGlmOut.cv2[[1]],
                                                                   test.data    = data.cont.cont$test.data,
                                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                                   noise.var    = data.cont.cont$noise.var,
                                                                   corr.split   = data.cont.cont$corr.split,
                                                                   where.split  = data.cont.cont$where.split,
                                                                   dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjTGlmOut.propscTGlmOut.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("18")

#####################################################################################################################
#################### 19. dr: True GLM Model outside, Mis func propensity score model outside, Cv2 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                     est.used          = "DR",
                                                                     type.var          = "cont",
                                                                     propsc.mod.out    = T,
                                                                     propsc.mthd       = "GLM",
                                                                     propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                     propsc.mod.insplt = NULL,
                                                                     adj.mod.out       = T, 
                                                                     adj.mthd          = "GLM", 
                                                                     adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                     adj.mod.insplt    = NULL, 
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)

final.tree.estdr.adjTGlmOut.propscNoisGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                    tree.list         = seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv2$tree.list,
                                                                    type.var          = "cont",
                                                                    seed              = a[job.number], 
                                                                    n.cv              = 5,
                                                                    propsc.mod.out    = T,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                    min.obs.propsc    = 10,
                                                                    adj.mod.out       = T,
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                    min.obs.adj       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmOut.propscNoisGlmOut.cv2[[1]],
                                                                      test.data    = data.cont.cont$test.data,
                                                                      true.trt.eff = data.cont.cont$true.trt.eff,
                                                                      noise.var    = data.cont.cont$noise.var,
                                                                      corr.split   = data.cont.cont$corr.split,
                                                                      where.split  = data.cont.cont$where.split,
                                                                      dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("19")

#####################################################################################################################
#################### 20. dr: Mis func GLM Model outside, True propensity score model outside, Cv2 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                     est.used          = "DR",
                                                                     type.var          = "cont",
                                                                     propsc.mod.out    = T,
                                                                     propsc.mthd       = "GLM",
                                                                     propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                     propsc.mod.insplt = NULL,
                                                                     adj.mod.out       = T, 
                                                                     adj.mthd          = "GLM", 
                                                                     adj.form.true     = NULL, 
                                                                     adj.mod.insplt    = NULL, 
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)

final.tree.estdr.adjNoisGlmOut.propscTGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                    tree.list         = seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv2$tree.list,
                                                                    type.var          = "cont",
                                                                    seed              = a[job.number], 
                                                                    n.cv              = 5,
                                                                    propsc.mod.out    = T,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                    min.obs.propsc    = 5,
                                                                    adj.mod.out       = T,
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = NULL, 
                                                                    min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmOut.propscTGlmOut.cv2[[1]],
                                                                      test.data    = data.cont.cont$test.data,
                                                                      true.trt.eff = data.cont.cont$true.trt.eff,
                                                                      noise.var    = data.cont.cont$noise.var,
                                                                      corr.split   = data.cont.cont$corr.split,
                                                                      where.split  = data.cont.cont$where.split,
                                                                      dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("20")

#####################################################################################################################
################# 21. dr: Mis func GLM Model outside, Mis func propensity score model outside, Cv2 ##################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                        est.used          = "DR",
                                                                        type.var          = "cont",
                                                                        propsc.mod.out    = T,
                                                                        propsc.mthd       = "GLM",
                                                                        propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                        propsc.mod.insplt = NULL,
                                                                        adj.mod.out       = T, 
                                                                        adj.mthd          = "GLM", 
                                                                        adj.form.true     = NULL, 
                                                                        adj.mod.insplt    = NULL, 
                                                                        num.truc.obs      = 30,
                                                                        min.node          = 20)

final.tree.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                       tree.list         = seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2$tree.list,
                                                                       type.var          = "cont",
                                                                       seed              = a[job.number], 
                                                                       n.cv              = 5,
                                                                       propsc.mod.out    = T,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                       min.obs.propsc    = 10,
                                                                       adj.mod.out       = T,
                                                                       adj.mthd          = "GLM",
                                                                       adj.form.true     = NULL, 
                                                                       min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2[[1]],
                                                                         test.data    = data.cont.cont$test.data,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("21")

#####################################################################################################################
############ 22. dr: Unmeasured cov GLM Model outside, Unmeasured cov propensity score model outside, Cv2 ###########
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmOut.propscFGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                                  est.used          = "DR",
                                                                  type.var          = "cont",
                                                                  propsc.mod.out    = T,
                                                                  propsc.mthd       = "GLM",
                                                                  propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                  propsc.mod.insplt = NULL,
                                                                  adj.mod.out       = T, 
                                                                  adj.mthd          = "GLM", 
                                                                  adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                  adj.mod.insplt    = NULL, 
                                                                  num.truc.obs      = 30,
                                                                  min.node          = 20)

final.tree.estdr.adjFGlmOut.propscFGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                                 tree.list         = seq.created.estdr.adjFGlmOut.propscFGlmOut.cv2$tree.list,
                                                                 type.var          = "cont",
                                                                 seed              = a[job.number], 
                                                                 n.cv              = 5,
                                                                 propsc.mod.out    = T,
                                                                 propsc.mthd       = "GLM",
                                                                 propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                 min.obs.propsc    = 10,
                                                                 adj.mod.out       = T,
                                                                 adj.mthd          = "GLM",
                                                                 adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                 min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmOut.propscFGlmOut.cv2[[1]],
                                                                   test.data    = test.data.mis,
                                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                                   noise.var    = data.cont.cont$noise.var,
                                                                   corr.split   = data.cont.cont$corr.split,
                                                                   where.split  = data.cont.cont$where.split,
                                                                   dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjFGlmOut.propscFGlmOut.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("22")

performance.hetero.drOut <- list(adjTGlmOut.propscTGlmOut.cv1       = eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1,
                                 adjTGlmOut.propscNoisGlmOut.cv1    = eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1,
                                 adjNoisGlmOut.propscTGlmOut.cv1    = eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1,
                                 adjNoisGlmOut.propscNoisGlmOut.cv1 = eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1,
                                 adjFGlmOut.propscFGlmOut.cv1       = eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1,
                                 adjTGlmOut.propscTGlmOut.cv2       = eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2,
                                 adjTGlmOut.propscNoisGlmOut.cv2    = eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2,
                                 adjNoisGlmOut.propscTGlmOut.cv2    = eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2,
                                 adjNoisGlmOut.propscNoisGlmOut.cv2 = eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2,
                                 adjFGlmOut.propscFGlmOut.cv2       = eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2)



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
######################### 23. ipw: GLM Model, outside node, True propensity score model, cv1 ########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscout.true.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                             est.used          = "IPW",
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T,
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = "A ~ X1 + X2 + X3",
                                                             w                 = NULL,
                                                             propsc.mod.insplt = NULL,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscout.true.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                             tree.list         = seq.created.estipw.glm.propscout.true.cv1$tree.list, 
                                                             lambda.used       = qchisq(0.95, 1), 
                                                             val.sample        = data.validation.cont.cont, 
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T, 
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = "A ~ X1 + X2 + X3",
                                                             val.w             = NULL,
                                                             propsc.mod.insplt = NULL, 
                                                             min.obs.mod       = 5)
t1 <- Sys.time()

eval.final.estipw.glm.propscout.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscout.true.cv1[[1]], 
                                                              test.data    = data.cont.cont$test.data,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split,
                                                              where.split  = data.cont.cont$where.split,
                                                              dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscout.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("23")

#####################################################################################################################
####################### 24. ipw: GLM Model, outside node, Mis func propensity score model, cv1 ######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscout.nois.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                             est.used          = "IPW",
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T,
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)", 
                                                             w                 = NULL, 
                                                             propsc.mod.insplt = NULL,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscout.nois.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                             tree.list         = seq.created.estipw.glm.propscout.nois.cv1$tree.list, 
                                                             lambda.used       = qchisq(0.95, 1), 
                                                             val.sample        = data.validation.cont.cont, 
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T, 
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)", 
                                                             val.w             = NULL,
                                                             propsc.mod.insplt = NULL, 
                                                             min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estipw.glm.propscout.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscout.nois.cv1[[1]], 
                                                              test.data    = data.cont.cont$test.data,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split,
                                                              where.split  = data.cont.cont$where.split,
                                                              dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscout.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("24")

#####################################################################################################################
#################### 25. ipw: GLM Model, outside node, Unmeasured cov propensity score model, cv1 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscout.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                            est.used          = "IPW",
                                                            type.var          = "cont",
                                                            propsc.mod.out    = T,
                                                            propsc.mthd       = "GLM", 
                                                            propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6", 
                                                            w                 = NULL, 
                                                            propsc.mod.insplt = NULL,
                                                            num.truc.obs      = 30,
                                                            min.node          = 20)

final.tree.estipw.glm.propscout.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont.mis, 
                                                            tree.list         = seq.created.estipw.glm.propscout.mis.cv1$tree.list, 
                                                            lambda.used       = qchisq(0.95, 1), 
                                                            val.sample        = data.validation.cont.cont.mis, 
                                                            type.var          = "cont",
                                                            propsc.mod.out    = T, 
                                                            propsc.mthd       = "GLM", 
                                                            propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6", 
                                                            val.w             = NULL,
                                                            propsc.mod.insplt = NULL, 
                                                            min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estipw.glm.propscout.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscout.mis.cv1[[1]], 
                                                             test.data    = test.data.mis,
                                                             true.trt.eff = data.cont.cont$true.trt.eff,
                                                             noise.var    = data.cont.cont$noise.var,
                                                             corr.split   = data.cont.cont$corr.split,
                                                             where.split  = data.cont.cont$where.split,
                                                             dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscout.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("25")

#####################################################################################################################
########################## 26. ipw: GLM Model, outside node, True propensity score model, Cv2 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscout.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                             est.used          = "IPW",
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T,
                                                             propsc.mthd       = "GLM",
                                                             propsc.form.true  = "A ~ X1 + X2 + X3",
                                                             w                 = NULL,
                                                             propsc.mod.insplt = NULL,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscout.true.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont,
                                                             tree.list        = seq.created.estipw.glm.propscout.true.cv2$tree.list,
                                                             type.var         = "cont",
                                                             seed             = a[job.number],
                                                             n.cv             = 5,
                                                             propsc.mod.out   = T, 
                                                             propsc.mthd      = "GLM", 
                                                             propsc.form.true = "A ~ X1 + X2 + X3",
                                                             min.obs.mod      = NULL)
t1 <- Sys.time()

eval.final.estipw.glm.propscout.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscout.true.cv2[[1]],
                                                              test.data    = data.cont.cont$test.data,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split,
                                                              where.split  = data.cont.cont$where.split,
                                                              dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscout.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("26")

#####################################################################################################################
####################### 27. ipw: GLM Model, outside node, Mis func propensity score model, Cv2 ######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscout.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                             est.used          = "IPW",
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T,
                                                             propsc.mthd       = "GLM",
                                                             propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                             w                 = NULL,
                                                             propsc.mod.insplt = NULL,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscout.nois.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont,
                                                             tree.list        = seq.created.estipw.glm.propscout.nois.cv2$tree.list,
                                                             type.var         = "cont",
                                                             seed             = a[job.number],
                                                             n.cv             = 5,
                                                             propsc.mod.out   = T, 
                                                             propsc.mthd      = "GLM", 
                                                             propsc.form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                             min.obs.mod      = NULL)
t1 <- Sys.time()

eval.final.estipw.glm.propscout.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscout.nois.cv2[[1]],
                                                              test.data    = data.cont.cont$test.data,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split,
                                                              where.split  = data.cont.cont$where.split,
                                                              dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscout.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("27")

#####################################################################################################################
#################### 28. ipw: GLM Model, outside node, Unmeasured cov propensity score model, Cv2 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscout.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                            est.used          = "IPW",
                                                            type.var          = "cont",
                                                            propsc.mod.out    = T,
                                                            propsc.mthd       = "GLM",
                                                            propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                            w                 = NULL,
                                                            propsc.mod.insplt = NULL,
                                                            num.truc.obs      = 30,
                                                            min.node          = 20)

final.tree.estipw.glm.propscout.mis.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont.mis,
                                                            tree.list        = seq.created.estipw.glm.propscout.mis.cv2$tree.list,
                                                            type.var         = "cont",
                                                            seed             = a[job.number],
                                                            n.cv             = 5,
                                                            propsc.mod.out   = T, 
                                                            propsc.mthd      = "GLM", 
                                                            propsc.form.true = "A ~ X1 + X3 + X4 + X5 + X6",
                                                            min.obs.mod      = NULL)
t1 <- Sys.time()

eval.final.estipw.glm.propscout.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscout.mis.cv2[[1]],
                                                             test.data    = test.data.mis,
                                                             true.trt.eff = data.cont.cont$true.trt.eff,
                                                             noise.var    = data.cont.cont$noise.var,
                                                             corr.split   = data.cont.cont$corr.split,
                                                             where.split  = data.cont.cont$where.split,
                                                             dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscout.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("28")

performance.homo.ipw <- list(glm.propscout.true.cv1 = eval.final.estipw.glm.propscout.true.cv1,
                             glm.propscout.nois.cv1 = eval.final.estipw.glm.propscout.nois.cv1,
                             glm.propscout.mis.cv1  = eval.final.estipw.glm.propscout.mis.cv1,
                             glm.propscout.true.cv2 = eval.final.estipw.glm.propscout.true.cv2,
                             glm.propscout.nois.cv2 = eval.final.estipw.glm.propscout.nois.cv2,
                             glm.propscout.mis.cv2  = eval.final.estipw.glm.propscout.mis.cv2)

#####################################################################################################################
############################# 29. g: GLM Model, outside node, True adjustment model, cv1 ############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modout.true.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                        est.used          = "G",
                                                        type.var          = "cont",
                                                        adj.mod.out       = T,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                        w                 = NULL,
                                                        adj.mod.insplt    = NULL,
                                                        num.truc.obs      = 10,
                                                        min.node          = 10)

final.tree.estg.glm.modout.true.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                      tree.list         = seq.created.estg.glm.modout.true.cv1$tree.list,
                                                      lambda.used       = qchisq(0.95, 1),
                                                      val.sample        = data.validation.cont.cont,              
                                                      type.var          = "cont",
                                                      adj.mod.out       = T,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                      val.w             = NULL,
                                                      adj.mod.insplt    = NULL,
                                                      min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estg.glm.modout.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.true.cv1[[1]],
                                                         test.data    = data.cont.cont$test.data,
                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                         noise.var    = data.cont.cont$noise.var,
                                                         corr.split   = data.cont.cont$corr.split,
                                                         where.split  = data.cont.cont$where.split,
                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modout.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("29")

#####################################################################################################################
########################### 30. g: GLM Model, outside node, Mis func adjustment model, cv1 ##########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modout.nois.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                        est.used          = "G",
                                                        type.var          = "cont",
                                                        adj.mod.out       = T,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = NULL,
                                                        w                 = NULL,
                                                        adj.mod.insplt    = NULL,
                                                        num.truc.obs      = 15,
                                                        min.node          = 15)

final.tree.estg.glm.modout.nois.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                      tree.list         = seq.created.estg.glm.modout.nois.cv1$tree.list,
                                                      lambda.used       = qchisq(0.95, 1),
                                                      val.sample        = data.validation.cont.cont,              
                                                      type.var          = "cont",
                                                      adj.mod.out       = T,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = NULL,
                                                      val.w             = NULL,
                                                      adj.mod.insplt    = NULL,
                                                      min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estg.glm.modout.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.nois.cv1[[1]],
                                                         test.data    = data.cont.cont$test.data,
                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                         noise.var    = data.cont.cont$noise.var,
                                                         corr.split   = data.cont.cont$corr.split,
                                                         where.split  = data.cont.cont$where.split,
                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modout.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("30")

#####################################################################################################################
######################## 31. g: GLM Model, outside node, Unmeasured cov adjustment model, cv1 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modout.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                       est.used          = "G",
                                                       type.var          = "cont",
                                                       adj.mod.out       = T,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                       w                 = NULL,
                                                       adj.mod.insplt    = NULL,
                                                       num.truc.obs      = 15,
                                                       min.node          = 15)

final.tree.estg.glm.modout.mis.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                     tree.list         = seq.created.estg.glm.modout.mis.cv1$tree.list,
                                                     lambda.used       = qchisq(0.95, 1),
                                                     val.sample        = data.validation.cont.cont.mis,              
                                                     type.var          = "cont",
                                                     adj.mod.out       = T,
                                                     adj.mthd          = "GLM",
                                                     adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                     val.w             = NULL,
                                                     adj.mod.insplt    = NULL,
                                                     min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estg.glm.modout.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.mis.cv1[[1]],
                                                        test.data    = test.data.mis,
                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                        noise.var    = data.cont.cont$noise.var,
                                                        corr.split   = data.cont.cont$corr.split,
                                                        where.split  = data.cont.cont$where.split,
                                                        dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modout.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("31")

#####################################################################################################################
############################# 32. g: GLM Model, outside node, True adjustment model, Cv2 ############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modout.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                        est.used          = "G",
                                                        type.var          = "cont",
                                                        adj.mod.out       = T,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                        w                 = NULL,
                                                        adj.mod.insplt    = NULL,
                                                        num.truc.obs      = 10,
                                                        min.node          = 10)

final.tree.estg.glm.modout.true.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                      tree.list         = seq.created.estg.glm.modout.true.cv2$tree.list,             
                                                      type.var          = "cont",
                                                      seed              = a[job.number],
                                                      n.cv              = 5,
                                                      adj.mod.out       = T,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                      min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estg.glm.modout.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.true.cv2[[1]],
                                                         test.data    = data.cont.cont$test.data,
                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                         noise.var    = data.cont.cont$noise.var,
                                                         corr.split   = data.cont.cont$corr.split,
                                                         where.split  = data.cont.cont$where.split,
                                                         dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modout.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("32")

#####################################################################################################################
########################### 33. g: GLM Model, outside node, Mis func adjustment model, Cv2 ##########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modout.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                        est.used          = "G",
                                                        type.var          = "cont",
                                                        adj.mod.out       = T,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = NULL,
                                                        w                 = NULL,
                                                        adj.mod.insplt    = NULL,
                                                        num.truc.obs      = 15,
                                                        min.node          = 15)

final.tree.estg.glm.modout.nois.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                      tree.list         = seq.created.estg.glm.modout.nois.cv2$tree.list,             
                                                      type.var          = "cont",
                                                      seed              = a[job.number],
                                                      n.cv              = 5,
                                                      adj.mod.out       = T,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = NULL,
                                                      min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estg.glm.modout.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.nois.cv2[[1]],
                                                         test.data    = data.cont.cont$test.data,
                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                         noise.var    = data.cont.cont$noise.var,
                                                         corr.split   = data.cont.cont$corr.split,
                                                         where.split  = data.cont.cont$where.split,
                                                         dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modout.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("33")

#####################################################################################################################
######################### 34. g: GLM Model, outside node, Unmeasured cov adjustment model, Cv2 ######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modout.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                       est.used          = "G",
                                                       type.var          = "cont",
                                                       adj.mod.out       = T,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                       w                 = NULL,
                                                       adj.mod.insplt    = NULL,
                                                       num.truc.obs      = 15,
                                                       min.node          = 15)

final.tree.estg.glm.modout.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                     tree.list         = seq.created.estg.glm.modout.mis.cv2$tree.list,             
                                                     type.var          = "cont",
                                                     seed              = a[job.number],
                                                     n.cv              = 5,
                                                     adj.mod.out       = T,
                                                     adj.mthd          = "GLM",
                                                     adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                     min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estg.glm.modout.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.mis.cv2[[1]],
                                                        test.data    = test.data.mis,
                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                        noise.var    = data.cont.cont$noise.var,
                                                        corr.split   = data.cont.cont$corr.split,
                                                        where.split  = data.cont.cont$where.split,
                                                        dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modout.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("34")

performance.homo.g <- list(glm.modout.true.cv1 = eval.final.estg.glm.modout.true.cv1,
                           glm.modout.nois.cv1 = eval.final.estg.glm.modout.nois.cv1,
                           glm.modout.mis.cv1  = eval.final.estg.glm.modout.mis.cv1,
                           glm.modout.true.cv2 = eval.final.estg.glm.modout.true.cv2,
                           glm.modout.nois.cv2 = eval.final.estg.glm.modout.nois.cv2,
                           glm.modout.mis.cv2  = eval.final.estg.glm.modout.mis.cv2)

#####################################################################################################################
###################### 35. dr: True GLM Model outside, True propensity score model outside, Cv1 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmOut.propscTGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                  est.used          = "DR",
                                                                  type.var          = "cont",
                                                                  propsc.mod.out    = T,
                                                                  propsc.mthd       = "GLM",
                                                                  propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                  propsc.mod.insplt = NULL,
                                                                  adj.mod.out       = T, 
                                                                  adj.mthd          = "GLM", 
                                                                  adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                  adj.mod.insplt    = NULL, 
                                                                  num.truc.obs      = 30,
                                                                  min.node          = 20)

final.tree.estdr.adjTGlmOut.propscTGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                 tree.list         = seq.created.estdr.adjTGlmOut.propscTGlmOut.cv1$tree.list, 
                                                                 lambda.used       = qchisq(0.95, 1),
                                                                 val.sample        = data.validation.cont.cont,
                                                                 type.var          = "cont",
                                                                 propsc.mod.out    = T,
                                                                 propsc.mthd       = "GLM",
                                                                 propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                 propsc.mod.insplt = NULL,
                                                                 adj.mod.out       = T,
                                                                 adj.mthd          = "GLM",
                                                                 adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                 adj.mod.insplt    = NULL,
                                                                 min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmOut.propscTGlmOut.cv1[[1]],
                                                                   test.data    = data.cont.cont$test.data,
                                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                                   noise.var    = data.cont.cont$noise.var,
                                                                   corr.split   = data.cont.cont$corr.split,
                                                                   where.split  = data.cont.cont$where.split,
                                                                   dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("35")

#####################################################################################################################
#################### 36. dr: True GLM Model outside, Mis func propensity score model outside, Cv1 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                     est.used          = "DR",
                                                                     type.var          = "cont",
                                                                     propsc.mod.out    = T,
                                                                     propsc.mthd       = "GLM",
                                                                     propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                     propsc.mod.insplt = NULL,
                                                                     adj.mod.out       = T, 
                                                                     adj.mthd          = "GLM", 
                                                                     adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                     adj.mod.insplt    = NULL, 
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)

final.tree.estdr.adjTGlmOut.propscNoisGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                    tree.list         = seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv1$tree.list, 
                                                                    lambda.used       = qchisq(0.95, 1),
                                                                    val.sample        = data.validation.cont.cont,
                                                                    type.var          = "cont",
                                                                    propsc.mod.out    = T,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                    propsc.mod.insplt = NULL,
                                                                    adj.mod.out       = T,
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                    adj.mod.insplt    = NULL,
                                                                    min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmOut.propscNoisGlmOut.cv1[[1]],
                                                                      test.data    = data.cont.cont$test.data,
                                                                      true.trt.eff = data.cont.cont$true.trt.eff,
                                                                      noise.var    = data.cont.cont$noise.var,
                                                                      corr.split   = data.cont.cont$corr.split,
                                                                      where.split  = data.cont.cont$where.split,
                                                                      dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("36")

#####################################################################################################################
#################### 37. dr: Mis func GLM Model outside, True propensity score model outside, Cv1 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                     est.used          = "DR",
                                                                     type.var          = "cont",
                                                                     propsc.mod.out    = T,
                                                                     propsc.mthd       = "GLM",
                                                                     propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                     propsc.mod.insplt = NULL,
                                                                     adj.mod.out       = T, 
                                                                     adj.mthd          = "GLM", 
                                                                     adj.form.true     = NULL, 
                                                                     adj.mod.insplt    = NULL, 
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)

final.tree.estdr.adjNoisGlmOut.propscTGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                    tree.list         = seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv1$tree.list, 
                                                                    lambda.used       = qchisq(0.95, 1),
                                                                    val.sample        = data.validation.cont.cont,
                                                                    type.var          = "cont",
                                                                    propsc.mod.out    = T,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                    propsc.mod.insplt = NULL,
                                                                    adj.mod.out       = T,
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = NULL, 
                                                                    adj.mod.insplt    = NULL,
                                                                    min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmOut.propscTGlmOut.cv1[[1]],
                                                                      test.data    = data.cont.cont$test.data,
                                                                      true.trt.eff = data.cont.cont$true.trt.eff,
                                                                      noise.var    = data.cont.cont$noise.var,
                                                                      corr.split   = data.cont.cont$corr.split,
                                                                      where.split  = data.cont.cont$where.split,
                                                                      dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("37")

#####################################################################################################################
################## 38. dr: Mis func GLM Model outside, Mis func propensity score model outside, Cv1 #################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                        est.used          = "DR",
                                                                        type.var          = "cont",
                                                                        propsc.mod.out    = T,
                                                                        propsc.mthd       = "GLM",
                                                                        propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                        propsc.mod.insplt = NULL,
                                                                        adj.mod.out       = T, 
                                                                        adj.mthd          = "GLM", 
                                                                        adj.form.true     = NULL, 
                                                                        adj.mod.insplt    = NULL, 
                                                                        num.truc.obs      = 30,
                                                                        min.node          = 20)

final.tree.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                       tree.list         = seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1$tree.list, 
                                                                       lambda.used       = qchisq(0.95, 1),
                                                                       val.sample        = data.validation.cont.cont,
                                                                       type.var          = "cont",
                                                                       propsc.mod.out    = T,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                       propsc.mod.insplt = NULL,
                                                                       adj.mod.out       = T,
                                                                       adj.mthd          = "GLM",
                                                                       adj.form.true     = NULL, 
                                                                       adj.mod.insplt    = NULL,
                                                                       min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1[[1]],
                                                                         test.data    = data.cont.cont$test.data,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("38")

#####################################################################################################################
############ 39. dr: Unmeasured cov GLM Model outside, Unmeasured cov propensity score model outside, Cv1 ###########
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmOut.propscFGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                                  est.used          = "DR",
                                                                  type.var          = "cont",
                                                                  propsc.mod.out    = T,
                                                                  propsc.mthd       = "GLM",
                                                                  propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                  propsc.mod.insplt = NULL,
                                                                  adj.mod.out       = T, 
                                                                  adj.mthd          = "GLM", 
                                                                  adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                  adj.mod.insplt    = NULL, 
                                                                  num.truc.obs      = 30,
                                                                  min.node          = 20)

final.tree.estdr.adjFGlmOut.propscFGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                                 tree.list         = seq.created.estdr.adjFGlmOut.propscFGlmOut.cv1$tree.list, 
                                                                 lambda.used       = qchisq(0.95, 1),
                                                                 val.sample        = data.validation.cont.cont.mis,
                                                                 type.var          = "cont",
                                                                 propsc.mod.out    = T,
                                                                 propsc.mthd       = "GLM",
                                                                 propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                 propsc.mod.insplt = NULL,
                                                                 adj.mod.out       = T,
                                                                 adj.mthd          = "GLM",
                                                                 adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                 adj.mod.insplt    = NULL,
                                                                 min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmOut.propscFGlmOut.cv1[[1]],
                                                                   test.data    = test.data.mis,
                                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                                   noise.var    = data.cont.cont$noise.var,
                                                                   corr.split   = data.cont.cont$corr.split,
                                                                   where.split  = data.cont.cont$where.split,
                                                                   dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("39")

#####################################################################################################################
###################### 40. dr: True GLM Model outside, True propensity score model outside, Cv2 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmOut.propscTGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                  est.used          = "DR",
                                                                  type.var          = "cont",
                                                                  propsc.mod.out    = T,
                                                                  propsc.mthd       = "GLM",
                                                                  propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                  propsc.mod.insplt = NULL,
                                                                  adj.mod.out       = T, 
                                                                  adj.mthd          = "GLM", 
                                                                  adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                  adj.mod.insplt    = NULL, 
                                                                  num.truc.obs      = 30,
                                                                  min.node          = 20)

final.tree.estdr.adjTGlmOut.propscTGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                 tree.list         = seq.created.estdr.adjTGlmOut.propscTGlmOut.cv2$tree.list,
                                                                 type.var          = "cont",
                                                                 seed              = a[job.number], 
                                                                 n.cv              = 5,
                                                                 propsc.mod.out    = T,
                                                                 propsc.mthd       = "GLM",
                                                                 propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                 min.obs.propsc    = 5,
                                                                 adj.mod.out       = T,
                                                                 adj.mthd          = "GLM",
                                                                 adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                 min.obs.adj       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmOut.propscTGlmOut.cv2[[1]],
                                                                   test.data    = data.cont.cont$test.data,
                                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                                   noise.var    = data.cont.cont$noise.var,
                                                                   corr.split   = data.cont.cont$corr.split,
                                                                   where.split  = data.cont.cont$where.split,
                                                                   dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("40")

#####################################################################################################################
#################### 41. dr: True GLM Model outside, Mis func propensity score model outside, Cv2 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                     est.used          = "DR",
                                                                     type.var          = "cont",
                                                                     propsc.mod.out    = T,
                                                                     propsc.mthd       = "GLM",
                                                                     propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                     propsc.mod.insplt = NULL,
                                                                     adj.mod.out       = T, 
                                                                     adj.mthd          = "GLM", 
                                                                     adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                     adj.mod.insplt    = NULL, 
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)

final.tree.estdr.adjTGlmOut.propscNoisGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                    tree.list         = seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv2$tree.list,
                                                                    type.var          = "cont",
                                                                    seed              = a[job.number], 
                                                                    n.cv              = 5,
                                                                    propsc.mod.out    = T,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                    min.obs.propsc    = 10,
                                                                    adj.mod.out       = T,
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                    min.obs.adj       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmOut.propscNoisGlmOut.cv2[[1]],
                                                                      test.data    = data.cont.cont$test.data,
                                                                      true.trt.eff = data.cont.cont$true.trt.eff,
                                                                      noise.var    = data.cont.cont$noise.var,
                                                                      corr.split   = data.cont.cont$corr.split,
                                                                      where.split  = data.cont.cont$where.split,
                                                                      dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("41")

#####################################################################################################################
#################### 42. dr: Mis func GLM Model outside, True propensity score model outside, Cv2 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                     est.used          = "DR",
                                                                     type.var          = "cont",
                                                                     propsc.mod.out    = T,
                                                                     propsc.mthd       = "GLM",
                                                                     propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                     propsc.mod.insplt = NULL,
                                                                     adj.mod.out       = T, 
                                                                     adj.mthd          = "GLM", 
                                                                     adj.form.true     = NULL, 
                                                                     adj.mod.insplt    = NULL, 
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)

final.tree.estdr.adjNoisGlmOut.propscTGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                    tree.list         = seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv2$tree.list,
                                                                    type.var          = "cont",
                                                                    seed              = a[job.number], 
                                                                    n.cv              = 5,
                                                                    propsc.mod.out    = T,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                    min.obs.propsc    = 5,
                                                                    adj.mod.out       = T,
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = NULL, 
                                                                    min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmOut.propscTGlmOut.cv2[[1]],
                                                                      test.data    = data.cont.cont$test.data,
                                                                      true.trt.eff = data.cont.cont$true.trt.eff,
                                                                      noise.var    = data.cont.cont$noise.var,
                                                                      corr.split   = data.cont.cont$corr.split,
                                                                      where.split  = data.cont.cont$where.split,
                                                                      dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("42")

#####################################################################################################################
################# 43. dr: Mis func GLM Model outside, Mis func propensity score model outside, Cv2 ##################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                        est.used          = "DR",
                                                                        type.var          = "cont",
                                                                        propsc.mod.out    = T,
                                                                        propsc.mthd       = "GLM",
                                                                        propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                        propsc.mod.insplt = NULL,
                                                                        adj.mod.out       = T, 
                                                                        adj.mthd          = "GLM", 
                                                                        adj.form.true     = NULL, 
                                                                        adj.mod.insplt    = NULL, 
                                                                        num.truc.obs      = 30,
                                                                        min.node          = 20)

final.tree.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                       tree.list         = seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2$tree.list,
                                                                       type.var          = "cont",
                                                                       seed              = a[job.number], 
                                                                       n.cv              = 5,
                                                                       propsc.mod.out    = T,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                       min.obs.propsc    = 10,
                                                                       adj.mod.out       = T,
                                                                       adj.mthd          = "GLM",
                                                                       adj.form.true     = NULL, 
                                                                       min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2[[1]],
                                                                         test.data    = data.cont.cont$test.data,
                                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                                         noise.var    = data.cont.cont$noise.var,
                                                                         corr.split   = data.cont.cont$corr.split,
                                                                         where.split  = data.cont.cont$where.split,
                                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("43")

#####################################################################################################################
############ 44. dr: Unmeasured cov GLM Model outside, Unmeasured cov propensity score model outside, Cv2 ###########
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmOut.propscFGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                                  est.used          = "DR",
                                                                  type.var          = "cont",
                                                                  propsc.mod.out    = T,
                                                                  propsc.mthd       = "GLM",
                                                                  propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                  propsc.mod.insplt = NULL,
                                                                  adj.mod.out       = T, 
                                                                  adj.mthd          = "GLM", 
                                                                  adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                  adj.mod.insplt    = NULL, 
                                                                  num.truc.obs      = 30,
                                                                  min.node          = 20)

final.tree.estdr.adjFGlmOut.propscFGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                                 tree.list         = seq.created.estdr.adjFGlmOut.propscFGlmOut.cv2$tree.list,
                                                                 type.var          = "cont",
                                                                 seed              = a[job.number], 
                                                                 n.cv              = 5,
                                                                 propsc.mod.out    = T,
                                                                 propsc.mthd       = "GLM",
                                                                 propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                 min.obs.propsc    = 10,
                                                                 adj.mod.out       = T,
                                                                 adj.mthd          = "GLM",
                                                                 adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                 min.obs.adj       = 15)
t1 <- Sys.time()

eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmOut.propscFGlmOut.cv2[[1]],
                                                                   test.data    = test.data.mis,
                                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                                   noise.var    = data.cont.cont$noise.var,
                                                                   corr.split   = data.cont.cont$corr.split,
                                                                   where.split  = data.cont.cont$where.split,
                                                                   dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("44")

performance.homo.drOut <- list(adjTGlmOut.propscTGlmOut.cv1       = eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1,
                               adjTGlmOut.propscNoisGlmOut.cv1    = eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1,
                               adjNoisGlmOut.propscTGlmOut.cv1    = eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1,
                               adjNoisGlmOut.propscNoisGlmOut.cv1 = eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1,
                               adjFGlmOut.propscFGlmOut.cv1       = eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1,
                               adjTGlmOut.propscTGlmOut.cv2       = eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2,
                               adjTGlmOut.propscNoisGlmOut.cv2    = eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2,
                               adjNoisGlmOut.propscTGlmOut.cv2    = eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2,
                               adjNoisGlmOut.propscNoisGlmOut.cv2 = eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2,
                               adjFGlmOut.propscFGlmOut.cv2       = eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2)

file.name = paste("../Data/AppendixC5Out/", toString(job.number), ".RData", sep = "")
save(performance.hetero.ipw, performance.homo.ipw,
     performance.hetero.g, performance.homo.g,
     performance.hetero.drOut, performance.homo.drOut,
     file = file.name)
