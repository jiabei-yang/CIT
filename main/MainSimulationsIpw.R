





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

