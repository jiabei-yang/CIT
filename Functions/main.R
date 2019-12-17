functions <- list.files(getwd())
functions <- grep("main.R", functions, invert = T, value = T)
for (i in functions){
  source(i)
}

data.cont.cont            <- makeData.cont.eff.cont(1000, 1000)
data.used.full.cont.cont  <- data.cont.cont$data.used
data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ]

# ipw: GLM Model, outside node, True propensity score model
t0 <- Sys.time()
seq.created.estipw.glm.propscout.true <- create.sequence(data.used         = data.used.cont.cont,
                                                         est.used          = "IPW",
                                                         type.var          = "cont",
                                                         propsc.mod.out    = T,
                                                         propsc.mthd       = "GLM", 
                                                         propsc.form.true  = "A ~ X1 + X2 + X3", 
                                                         w                 = cbind(1, data.used.cont.cont[, 3:5]), 
                                                         propsc.mod.insplt = NULL,
                                                         num.truc.obs      = 30,
                                                         min.node          = 20)

# numtrees.ipw.glm.propscout.true <- length(seq.created.estipw.glm.propscout.true$tree.list)

final.tree.estipw.glm.propscout.true.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                             tree.list         = seq.created.estipw.glm.propscout.true$tree.list, 
                                                             lambda.used       = 4, 
                                                             val.sample        = data.validation.cont.cont, 
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T, 
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = "A ~ X1 + X2 + X3", 
                                                             val.w             = cbind(1, data.validation.cont.cont[, 3:5]),
                                                             propsc.mod.insplt = NULL, 
                                                             min.obs.mod       = NULL)
t1 <- Sys.time()

final.tree.estipw.glm.propscout.true.cv1

data.used         = data.used.cont.cont
tree.list         = seq.created.estipw.glm.propscout.true$tree.list
lambda.used       = 4
val.sample        = data.validation.cont.cont
type.var          = "cont"
propsc.mod.out    = T
propsc.mthd       = "GLM"
propsc.form.true  = "A ~ X1 + X2 + X3"
val.w             = cbind(1, data.validation.cont.cont[, 3:5])
propsc.mod.insplt = NULL
min.obs.mod       = NULL

eval.final.estipw.glm.propscout.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscout.true.cv1[[1]], 
                                                              test.data    = data.cont.cont$test.data,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split)
eval.final.estipw.glm.propscout.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))

# final.tree   = seq.created.estipw.glm.propscout.true$tree.list[[14]]
# final.tree   = final.tree.estipw.glm.propscout.true.cv1[[1]]
# test.data    = data.cont.cont$test.data
# true.trt.eff = data.cont.cont$true.trt.eff
# noise.var    = data.cont.cont$noise.var
# corr.split   = data.cont.cont$corr.split

# ipw: GLM Model, inside node, True propensity score model
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.true <- create.sequence(data.used         = data.used.cont.cont,
                                                          est.used          = "IPW",
                                                          type.var          = "cont",
                                                          propsc.mod.out    = F,
                                                          propsc.mthd       = "GLM", 
                                                          propsc.form.true  = "A ~ X1 + X2 + X3", 
                                                          w                 = cbind(1, data.used.cont.cont[, 3:5]), 
                                                          propsc.mod.insplt = F,
                                                          num.truc.obs      = 30,
                                                          min.node          = 20)
# numtrees.ipw.glm.propscinnd.true <- length(seq.created.estipw.glm.propscinnd.true$tree.list)

final.tree.estipw.glm.propscinnd.true.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                              tree.list         = seq.created.estipw.glm.propscinnd.true$tree.list, 
                                                              lambda.used       = 4, 
                                                              val.sample        = data.validation.cont.cont, 
                                                              type.var          = "cont",
                                                              propsc.mod.out    = F, 
                                                              propsc.mthd       = "GLM", 
                                                              propsc.form.true  = "A ~ X1 + X2 + X3", 
                                                              val.w             = cbind(1, data.validation.cont.cont[, 3:5]),
                                                              propsc.mod.insplt = F, 
                                                              min.obs.mod       = 5)
t1 <- Sys.time()
final.tree.estipw.glm.propscinnd.true.cv1

data.used         = data.used.cont.cont
tree.list         = seq.created.estipw.glm.propscinnd.true$tree.list
lambda.used       = 4
val.sample        = data.validation.cont.cont
type.var          = "cont"
propsc.mod.out    = F
propsc.mthd       = "GLM"
propsc.form.true  = "A ~ X1 + X2 + X3"
val.w             = cbind(1, data.validation.cont.cont[, 3:5])
propsc.mod.insplt = F
num.truc.obs      = NULL
min.node          = NULL

eval.final.estipw.glm.propscinnd.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.true.cv1[[1]], 
                                                               test.data    = data.cont.cont$test.data,
                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                               noise.var    = data.cont.cont$noise.var,
                                                               corr.split   = data.cont.cont$corr.split)
eval.final.estipw.glm.propscinnd.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))

# ipw: GLM Model, inside split, True propensity score model
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.true <- create.sequence(data.used         = data.used.cont.cont,
                                                            est.used          = "IPW",
                                                            type.var          = "cont",
                                                            propsc.mod.out    = F,
                                                            propsc.mthd       = "GLM", 
                                                            propsc.form.true  = "A ~ X1 + X2 + X3", 
                                                            w                 = cbind(1, data.used.cont.cont[, 3:5]), 
                                                            propsc.mod.insplt = T,
                                                            num.truc.obs      = 30,
                                                            min.node          = 20)
# numtrees.ipw.glm.propscinsplt.true <- length(seq.created.estipw.glm.propscinsplt.true$tree.list)

final.tree.estipw.glm.propscinsplt.true.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                                tree.list         = seq.created.estipw.glm.propscinsplt.true$tree.list, 
                                                                lambda.used       = 4, 
                                                                val.sample        = data.validation.cont.cont, 
                                                                type.var          = "cont",
                                                                propsc.mod.out    = F, 
                                                                propsc.mthd       = "GLM", 
                                                                propsc.form.true  = "A ~ X1 + X2 + X3", 
                                                                val.w             = cbind(1, data.validation.cont.cont[, 3:5]),
                                                                propsc.mod.insplt = T, 
                                                                min.obs.mod       = 5)
t1 <- Sys.time()
final.tree.estipw.glm.propscinsplt.true.cv1

data.used         = data.used.cont.cont
tree.list         = seq.created.estipw.glm.propscinsplt.true$tree.list
lambda.used       = 4
val.sample        = data.validation.cont.cont
type.var          = "cont"
propsc.mod.out    = F
propsc.mthd       = "GLM"
propsc.form.true  = "A ~ X1 + X2 + X3"
val.w             = cbind(1, data.validation.cont.cont[, 3:5])
propsc.mod.insplt = T
num.truc.obs      = NULL
min.node          = NULL

eval.final.estipw.glm.propscinsplt.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.true.cv1[[1]], 
                                                                 test.data    = data.cont.cont$test.data,
                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                 noise.var    = data.cont.cont$noise.var,
                                                                 corr.split   = data.cont.cont$corr.split)
eval.final.estipw.glm.propscinsplt.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))

# ipw: GLM Model, outside node, Misspecified propensity score model
t0 <- Sys.time()
seq.created.estipw.glm.propscout.mis <- create.sequence(data.used         = data.used.cont.cont,
                                                        est.used          = "IPW",
                                                        type.var          = "cont",
                                                        propsc.mod.out    = T,
                                                        propsc.mthd       = "GLM", 
                                                        propsc.form.true  = NULL, 
                                                        w                 = NULL, 
                                                        propsc.mod.insplt = NULL,
                                                        num.truc.obs      = 30,
                                                        min.node          = 20)
# numtrees.ipw.glm.propscout.mis <- length(seq.created.estipw.glm.propscout.mis$tree.list)

# val.sample: the order of the columns must be A, Y, X
final.tree.estipw.glm.propscout.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                             tree.list         = seq.created.estipw.glm.propscout.mis$tree.list, 
                                                             lambda.used       = 4, 
                                                             val.sample        = data.validation.cont.cont, 
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T, 
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = NULL, 
                                                             val.w             = NULL,
                                                             propsc.mod.insplt = NULL, 
                                                             min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estipw.glm.propscout.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscout.mis.cv1[[1]], 
                                                             test.data    = data.cont.cont$test.data,
                                                             true.trt.eff = data.cont.cont$true.trt.eff,
                                                             noise.var    = data.cont.cont$noise.var,
                                                             corr.split   = data.cont.cont$corr.split)
eval.final.estipw.glm.propscout.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))


# ipw: GLM Model, inside node, Misspecified propensity score model
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.mis <- create.sequence(data.used         = data.used.cont.cont,
                                                         est.used          = "IPW",
                                                         type.var          = "cont",
                                                         propsc.mod.out    = F,
                                                         propsc.mthd       = "GLM", 
                                                         propsc.form.true  = NULL, 
                                                         w                 = NULL, 
                                                         propsc.mod.insplt = F,
                                                         num.truc.obs      = 30,
                                                         min.node          = 20)
# numtrees.ipw.glm.propscinnd.mis <- length(seq.created.estipw.glm.propscinnd.mis$tree.list)

final.tree.estipw.glm.propscinnd.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                             tree.list         = seq.created.estipw.glm.propscinnd.mis$tree.list, 
                                                             lambda.used       = 4, 
                                                             val.sample        = data.validation.cont.cont, 
                                                             type.var          = "cont",
                                                             propsc.mod.out    = F, 
                                                             propsc.mthd       = "GLM", 
                                                             propsc.form.true  = NULL, 
                                                             val.w             = NULL,
                                                             propsc.mod.insplt = F, 
                                                             min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.mis.cv1[[1]], 
                                                              test.data    = data.cont.cont$test.data,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split)
eval.final.estipw.glm.propscinnd.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))


# ipw: GLM Model, inside split, Misspecified propensity score model
t0 <- Sys.time()
seq.created.estipw.glm.propscinsplt.mis <- create.sequence(data.used         = data.used.cont.cont,
                                                           est.used          = "IPW",
                                                           type.var          = "cont",
                                                           propsc.mod.out    = F,
                                                           propsc.mthd       = "GLM", 
                                                           propsc.form.true  = NULL, 
                                                           w                 = NULL, 
                                                           propsc.mod.insplt = T,
                                                           num.truc.obs      = 30,
                                                           min.node          = 20)
# numtrees.ipw.glm.propscinsplt.mis <- length(seq.created.estipw.glm.propscinsplt.mis$tree.list)

final.tree.estipw.glm.propscinsplt.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                               tree.list         = seq.created.estipw.glm.propscinsplt.mis$tree.list, 
                                                               lambda.used       = 4, 
                                                               val.sample        = data.validation.cont.cont, 
                                                               type.var          = "cont",
                                                               propsc.mod.out    = F, 
                                                               propsc.mthd       = "GLM", 
                                                               propsc.form.true  = NULL, 
                                                               val.w             = NULL,
                                                               propsc.mod.insplt = T, 
                                                               min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinsplt.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinsplt.mis.cv1[[1]], 
                                                                test.data    = data.cont.cont$test.data,
                                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                                noise.var    = data.cont.cont$noise.var,
                                                                corr.split   = data.cont.cont$corr.split)
eval.final.estipw.glm.propscinsplt.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))


# ipw: GAM Model, outside node, True propensity score model
t0 <- Sys.time()
seq.created.estipw.gam.propscout.true <- create.sequence(data.used         = data.used.cont.cont,
                                                         est.used          = "IPW",
                                                         type.var          = "cont",
                                                         propsc.mod.out    = T,
                                                         propsc.mthd       = "GAM",
                                                         propsc.form.true  = "A ~ s(X1) + s(X2) + s(X3)",
                                                         w                 = cbind(1, data.used.cont.cont[, 3:5]),
                                                         propsc.mod.insplt = NULL,
                                                         num.truc.obs      = 30,
                                                         min.node          = 20)
# numtrees.ipw.gam.propscout.true <- length(seq.created.estipw.gam.propscout.true$tree.list)

# data.used         = data.used.cont.cont
# est.used          = "IPW"
# type.var          = "cont"
# propsc.mod.out    = T
# propsc.mthd       = "GAM"
# propsc.form.true  = "A ~ s(X1) + s(X2) + s(X3)"
# w                 = cbind(1, data.used.cont.cont[, 3:5])
# propsc.mod.insplt = NULL

final.tree.estipw.gam.propscout.true.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                             tree.list         = seq.created.estipw.gam.propscout.true$tree.list, 
                                                             lambda.used       = 4, 
                                                             val.sample        = data.validation.cont.cont, 
                                                             type.var          = "cont",
                                                             propsc.mod.out    = T, 
                                                             propsc.mthd       = "GAM", 
                                                             propsc.form.true  = "A ~ s(X1) + s(X2) + s(X3)", 
                                                             val.w             = cbind(1, data.validation.cont.cont[, 3:5]),
                                                             propsc.mod.insplt = NULL, 
                                                             min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estipw.gam.propscout.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.gam.propscout.true.cv1[[1]], 
                                                              test.data    = data.cont.cont$test.data,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split)
eval.final.estipw.gam.propscout.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))

# ipw: GAM Model, inside node, True propensity score model
t0 <- Sys.time()
seq.created.estipw.gam.propscinnd.true <- create.sequence(data.used         = data.used.cont.cont,
                                                          est.used          = "IPW",
                                                          type.var          = "cont",
                                                          propsc.mod.out    = F,
                                                          propsc.mthd       = "GAM",
                                                          propsc.form.true  = "A ~ s(X1) + s(X2) + s(X3)",
                                                          w                 = cbind(1, data.used.cont.cont[, 3:5]),
                                                          propsc.mod.insplt = F,
                                                          num.truc.obs      = 30,
                                                          min.node          = 30)
# numtrees.ipw.gam.propscinnd.true <- length(seq.created.estipw.gam.propscinnd.true$tree.list)
# 
# data.used         = data.used.cont.cont
# est.used          = "IPW"
# type.var          = "cont"
# propsc.mod.out    = F
# propsc.mthd       = "GAM"
# propsc.form.true  = "A ~ s(X1) + s(X2) + s(X3)"
# w                 = cbind(1, data.used.cont.cont[, 3:5])
# propsc.mod.insplt = F

final.tree.estipw.gam.propscinnd.true.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                              tree.list         = seq.created.estipw.gam.propscinnd.true$tree.list, 
                                                              lambda.used       = 4, 
                                                              val.sample        = data.validation.cont.cont, 
                                                              type.var          = "cont",
                                                              propsc.mod.out    = F, 
                                                              propsc.mthd       = "GAM", 
                                                              propsc.form.true  = "A ~ s(X1) + s(X2) + s(X3)", 
                                                              val.w             = cbind(1, data.validation.cont.cont[, 3:5]),
                                                              propsc.mod.insplt = F, 
                                                              min.obs.mod       = 30)
t1 <- Sys.time()

eval.final.estipw.gam.propscinnd.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.gam.propscinnd.true.cv1[[1]], 
                                                               test.data    = data.cont.cont$test.data,
                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                               noise.var    = data.cont.cont$noise.var,
                                                               corr.split   = data.cont.cont$corr.split)
eval.final.estipw.gam.propscinnd.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))

# ipw: GAM Model, inside split, True propensity score model
t0 <- Sys.time()
seq.created.estipw.gam.propscinsplt.true <- create.sequence(data.used         = data.used.cont.cont,
                                                            est.used          = "IPW",
                                                            type.var          = "cont",
                                                            propsc.mod.out    = F,
                                                            propsc.mthd       = "GAM",
                                                            propsc.form.true  = "A ~ s(X1) + s(X2) + s(X3)",
                                                            w                 = cbind(1, data.used.cont.cont[, 3:5]),
                                                            propsc.mod.insplt = T,
                                                            num.truc.obs      = 30,
                                                            min.node          = 30)
# numtrees.ipw.gam.propscinsplt.true <- length(seq.created.estipw.gam.propscinsplt.true$tree.list)
# 
# data.used         = data.used.cont.cont
# est.used          = "IPW"
# type.var          = "cont"
# propsc.mod.out    = F
# propsc.mthd       = "GAM"
# propsc.form.true  = "A ~ s(X1) + s(X2) + s(X3)"
# w                 = cbind(1, data.used.cont.cont[, 3:5])
# propsc.mod.insplt = T
# num.truc.obs      = 30

final.tree.estipw.gam.propscinsplt.true.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                                tree.list         = seq.created.estipw.gam.propscinsplt.true$tree.list, 
                                                                lambda.used       = 4, 
                                                                val.sample        = data.validation.cont.cont, 
                                                                type.var          = "cont",
                                                                propsc.mod.out    = F, 
                                                                propsc.mthd       = "GAM", 
                                                                propsc.form.true  = "A ~ s(X1) + s(X2) + s(X3)", 
                                                                val.w             = cbind(1, data.validation.cont.cont[, 3:5]),
                                                                propsc.mod.insplt = T, 
                                                                min.obs.mod       = 30)
t1 <- Sys.time()

eval.final.estipw.gam.propscinsplt.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.gam.propscinsplt.true.cv1[[1]], 
                                                                 test.data    = data.cont.cont$test.data,
                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                 noise.var    = data.cont.cont$noise.var,
                                                                 corr.split   = data.cont.cont$corr.split)
eval.final.estipw.gam.propscinsplt.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))

# ipw: GAM Model, outside node, Misspecified propensity score model
t0 <- Sys.time()
seq.created.estipw.gam.propscout.mis <- create.sequence(data.used         = data.used.cont.cont,
                                                        est.used          = "IPW",
                                                        type.var          = "cont",
                                                        propsc.mod.out    = T,
                                                        propsc.mthd       = "GAM",
                                                        propsc.form.true  = NULL,
                                                        w                 = NULL,
                                                        propsc.mod.insplt = NULL,
                                                        num.truc.obs      = 30,
                                                        min.node          = 20)

final.tree.estipw.gam.propscout.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                            tree.list         = seq.created.estipw.gam.propscout.mis$tree.list, 
                                                            lambda.used       = 4, 
                                                            val.sample        = data.validation.cont.cont, 
                                                            type.var          = "cont",
                                                            propsc.mod.out    = T, 
                                                            propsc.mthd       = "GAM", 
                                                            propsc.form.true  = NULL, 
                                                            val.w             = NULL,
                                                            propsc.mod.insplt = NULL, 
                                                            min.obs.mod       = NULL)
t1 <- Sys.time()

eval.final.estipw.gam.propscout.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.gam.propscout.mis.cv1[[1]], 
                                                             test.data    = data.cont.cont$test.data,
                                                             true.trt.eff = data.cont.cont$true.trt.eff,
                                                             noise.var    = data.cont.cont$noise.var,
                                                             corr.split   = data.cont.cont$corr.split)
eval.final.estipw.gam.propscout.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))

# ipw: GAM Model, inside node, Misspecified propensity score model
t0 <- Sys.time()
seq.created.estipw.gam.propscinnd.mis <- create.sequence(data.used         = data.used.cont.cont,
                                                         est.used          = "IPW",
                                                         type.var          = "cont",
                                                         propsc.mod.out    = F,
                                                         propsc.mthd       = "GAM",
                                                         propsc.form.true  = NULL,
                                                         w                 = NULL,
                                                         propsc.mod.insplt = F,
                                                         num.truc.obs      = 30,
                                                         min.node          = 60)
# numtrees.ipw.gam.propscinnd.mis <- length(seq.created.estipw.gam.propscinnd.mis$tree.list)

final.tree.estipw.gam.propscinnd.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                             tree.list         = seq.created.estipw.gam.propscinnd.mis$tree.list, 
                                                             lambda.used       = 4, 
                                                             val.sample        = data.validation.cont.cont, 
                                                             type.var          = "cont",
                                                             propsc.mod.out    = F, 
                                                             propsc.mthd       = "GAM", 
                                                             propsc.form.true  = NULL, 
                                                             val.w             = NULL,
                                                             propsc.mod.insplt = F, 
                                                             min.obs.mod       = 60)
t1 <- Sys.time()

eval.final.estipw.gam.propscinnd.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.gam.propscinnd.mis.cv1[[1]], 
                                                              test.data    = data.cont.cont$test.data,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split)
eval.final.estipw.gam.propscinnd.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))

# ipw: GAM Model, inside split, Misspecified propensity score model
seq.created.estipw.gam.propscinsplt.mis <- create.sequence(data.used         = data.used.cont.cont,
                                                           est.used          = "IPW",
                                                           type.var          = "cont",
                                                           propsc.mod.out    = F,
                                                           propsc.mthd       = "GAM",
                                                           propsc.form.true  = NULL,
                                                           w                 = NULL,
                                                           propsc.mod.insplt = T,
                                                           num.truc.obs      = 60,
                                                           min.node          = 60)
# numtrees.ipw.gam.propscinsplt.mis <- length(seq.created.estipw.gam.propscinsplt.mis$tree.list)

final.tree.estipw.gam.propscinsplt.mis.cv1 <- EstIpw.CvMethod1(data.used         = data.used.cont.cont, 
                                                               tree.list         = seq.created.estipw.gam.propscinsplt.mis$tree.list, 
                                                               lambda.used       = 4, 
                                                               val.sample        = data.validation.cont.cont, 
                                                               type.var          = "cont",
                                                               propsc.mod.out    = F, 
                                                               propsc.mthd       = "GAM", 
                                                               propsc.form.true  = NULL, 
                                                               val.w             = NULL,
                                                               propsc.mod.insplt = T, 
                                                               min.obs.mod       = 60)
t1 <- Sys.time()

eval.final.estipw.gam.propscinsplt.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estipw.gam.propscinsplt.mis.cv1[[1]], 
                                                                test.data    = data.cont.cont$test.data,
                                                                true.trt.eff = data.cont.cont$true.trt.eff,
                                                                noise.var    = data.cont.cont$noise.var,
                                                                corr.split   = data.cont.cont$corr.split)
eval.final.estipw.gam.propscinsplt.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))








# g: GLM Model inside node
parms.used <- list(trt        = data.used.cont.cont$A,
                   covariates = data.used.cont.cont[, 3:dim(data.used.cont.cont)[2]],
                   response   = data.used.cont.cont$Y,
                   use.var    = "true",
                   mod.in     = T)   

ulist.used <- list(eval = etemp.g, split = stemp.g, init = itemp)
fit <- rpart(Y ~ X1 + X2 + X3 + X4 + X5 + X6, data = data.used.cont.cont,
             method   = ulist.used,
             parms    = parms.used,
             control  = rpart.control(cp = 0, minbucket = 30, maxsurrogate = 0, maxcompete = 0))

# g: GLM Model outside node
tmp <- est.cond.eff(data.used.cont.cont, method = "GLM")
parms.used <- list(trt            = data.used.cont.cont$A,
                   covariates     = data.used.cont.cont[, 3:dim(data.used.cont.cont)[2]],
                   response       = data.used.cont.cont$Y,
                   est.cond.eff.1 = tmp$pred.A.1,
                   est.cond.eff.0 = tmp$pred.A.0,
                   var.rb         = tmp$var.rb,
                   use.var        = "true",
                   mod.in         = F)   

fit <- rpart(Y ~ X1 + X2 + X3 + X4 + X5 + X6, data = data.used.cont.cont,
             method   = ulist.used,
             parms    = parms.used,
             control  = rpart.control(cp = 0, minbucket = 30, maxsurrogate = 0, maxcompete = 0))


# dr: Conditional effect & propensity score both GAM outside node
data.used.cont.cont.noy <- data.used.cont.cont[, -2]
tmp.1 <- est.cond.eff(df = data.used.cont.cont, method = "GAM")
tmp.2 <- est.prop.sc(df.noy = data.used.cont.cont.noy, method = "GAM")
parms.used <- list(trt            = data.used.cont.cont$A,
                   covariates     = data.used.cont.cont[, 3:dim(data.used.cont.cont)[2]],
                   response       = data.used.cont.cont$Y,
                   est.cond.eff.1 = tmp.1$pred.A.1,
                   est.cond.eff.0 = tmp.1$pred.A.0,
                   prop.sc        = tmp.2,
                   use.var        = "true")

ulist.used <- list(eval = etemp.dr, split = stemp.dr, init = itemp)
fit <- rpart(Y ~ X1 + X2 + X3 + X4 + X5 + X6, data = data.used.cont.cont,
             method   = ulist.used,
             parms    = parms.used,
             control  = rpart.control(cp = 0, minbucket = 30, maxsurrogate = 0, maxcompete = 0))





############################### Temporary code ################################
parms <- parms.used
x     <- sort(data.used.cont.cont$X1)
y     <- parms.used$response[order(data.used.cont.cont$X1)]
wt    <- rep(1, length(y))

# create.sequence
data.used         = data.used.cont.cont
est.used          = "IPW"
type.var          = "cont"
propsc.mod.out    = T
propsc.mthd       = "GLM"
# propsc.form.true  = "A ~ X1 + X2 + X3",
# w                 = cbind(1, data.used.cont.cont[, 3:5]),
propsc.form.true  = "A ~ X1"
w                 = cbind(1, data.used.cont.cont[, 3])
propsc.mod.insplt = NULL
num.truc.obs      = 30
min.node          = 20

data = data.used
method   = ulist.used
parms    = parms.used
control  = rpart.control(cp = 0, minsplit = num.truc.obs * 2, minbucket = min.node, maxsurrogate = 0, maxcompete = 0)


# EstIpw.CvMethod1
data.used         = data.used.cont.cont
tree.list         = seq.created.estipw.glm.propscinsplt.true$tree.list
lambda.used       = 4
val.sample        = data.validation.cont.cont
type.var          = "cont"
propsc.mod.out    = F
propsc.mthd       = "GLM"
# propsc.form.true  = "A ~ X1 + X2 + X3"
# val.w             = cbind(1, data.validation.cont.cont[, 3:5]),
propsc.form.true  = "A ~ X1"
val.w             = cbind(1, data.validation.cont.cont[, 3])
propsc.mod.insplt = T
min.obs.mod       = 5

