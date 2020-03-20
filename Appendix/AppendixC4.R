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
########################## 1. ipw: GLM Model, inside node, True propensity score model, Cv2 #########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                              est.used          = "IPW",
                                                              type.var          = "cont",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ X1 + X2 + X3",
                                                              w                 = NULL,
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.true.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont,
                                                              tree.list        = seq.created.estipw.glm.propscinnd.true.cv2$tree.list,
                                                              type.var         = "cont",
                                                              seed             = a[job.number],
                                                              n.cv             = 5,
                                                              propsc.mod.out   = F, 
                                                              propsc.mthd      = "GLM", 
                                                              propsc.form.true = "A ~ X1 + X2 + X3",
                                                              min.obs.mod      = 5)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.true.cv2[[1]],
                                                               test.data    = data.cont.cont$test.data,
                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                               noise.var    = data.cont.cont$noise.var,
                                                               corr.split   = data.cont.cont$corr.split,
                                                               where.split  = data.cont.cont$where.split,
                                                               dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinnd.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinnd.true.cv2$corr.frst.splt <- as.character(seq.created.estipw.glm.propscinnd.true.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("1")

#####################################################################################################################
######################### 2. ipw: GLM Model, inside node, Mis func propensity score model, Cv2 ######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                              est.used          = "IPW",
                                                              type.var          = "cont",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                              w                 = NULL,
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.nois.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont,
                                                              tree.list        = seq.created.estipw.glm.propscinnd.nois.cv2$tree.list,
                                                              type.var         = "cont",
                                                              seed             = a[job.number],
                                                              n.cv             = 5,
                                                              propsc.mod.out   = F, 
                                                              propsc.mthd      = "GLM", 
                                                              propsc.form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                              min.obs.mod      = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.nois.cv2[[1]],
                                                               test.data    = data.cont.cont$test.data,
                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                               noise.var    = data.cont.cont$noise.var,
                                                               corr.split   = data.cont.cont$corr.split,
                                                               where.split  = data.cont.cont$where.split,
                                                               dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinnd.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinnd.nois.cv2$corr.frst.splt <- as.character(seq.created.estipw.glm.propscinnd.nois.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("2")

#####################################################################################################################
##################### 3. ipw: GLM Model, inside node, Unmeasured cov propensity score model, Cv2 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                             est.used          = "IPW",
                                                             type.var          = "cont",
                                                             propsc.mod.out    = F,
                                                             propsc.mthd       = "GLM",
                                                             propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                             w                 = NULL,
                                                             propsc.mod.insplt = F,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscinnd.mis.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont.mis,
                                                             tree.list        = seq.created.estipw.glm.propscinnd.mis.cv2$tree.list,
                                                             type.var         = "cont",
                                                             seed             = a[job.number],
                                                             n.cv             = 5,
                                                             propsc.mod.out   = F, 
                                                             propsc.mthd      = "GLM", 
                                                             propsc.form.true = "A ~ X1 + X3 + X4 + X5 + X6",
                                                             min.obs.mod      = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.mis.cv2[[1]],
                                                              test.data    = test.data.mis,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split,
                                                              where.split  = data.cont.cont$where.split,
                                                              dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinnd.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estipw.glm.propscinnd.mis.cv2$corr.frst.splt <- as.character(seq.created.estipw.glm.propscinnd.mis.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("3")

performance.hetero.ipw <- list(glm.propscinnd.true.cv2 = eval.final.estipw.glm.propscinnd.true.cv2,
                               glm.propscinnd.nois.cv2 = eval.final.estipw.glm.propscinnd.nois.cv2,
                               glm.propscinnd.mis.cv2  = eval.final.estipw.glm.propscinnd.mis.cv2)

#####################################################################################################################
############################# 4. g: GLM Model, inside node, True adjustment model, Cv2 ##############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                         est.used          = "G",
                                                         type.var          = "cont",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 10,
                                                         min.node          = 10)

final.tree.estg.glm.modinnd.true.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                       tree.list         = seq.created.estg.glm.modinnd.true.cv2$tree.list,             
                                                       type.var          = "cont",
                                                       seed              = a[job.number],
                                                       n.cv              = 5,
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                       min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.true.cv2[[1]],
                                                          test.data    = data.cont.cont$test.data,
                                                          true.trt.eff = data.cont.cont$true.trt.eff,
                                                          noise.var    = data.cont.cont$noise.var,
                                                          corr.split   = data.cont.cont$corr.split,
                                                          where.split  = data.cont.cont$where.split,
                                                          dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modinnd.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinnd.true.cv2$corr.frst.splt <- as.character(seq.created.estg.glm.modinnd.true.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("4")

#####################################################################################################################
############################ 5. g: GLM Model, inside node, Mis func adjustment model, Cv2 ###########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                         est.used          = "G",
                                                         type.var          = "cont",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = NULL,
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 15,
                                                         min.node          = 15)

final.tree.estg.glm.modinnd.nois.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                       tree.list         = seq.created.estg.glm.modinnd.nois.cv2$tree.list,             
                                                       type.var          = "cont",
                                                       seed              = a[job.number],
                                                       n.cv              = 5,
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = NULL,
                                                       min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.nois.cv2[[1]],
                                                          test.data    = data.cont.cont$test.data,
                                                          true.trt.eff = data.cont.cont$true.trt.eff,
                                                          noise.var    = data.cont.cont$noise.var,
                                                          corr.split   = data.cont.cont$corr.split,
                                                          where.split  = data.cont.cont$where.split,
                                                          dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modinnd.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinnd.nois.cv2$corr.frst.splt <- as.character(seq.created.estg.glm.modinnd.nois.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("5")

#####################################################################################################################
######################## 6. g: GLM Model, inside node, Unmeasured cov adjustment model, Cv2 #########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                        est.used          = "G",
                                                        type.var          = "cont",
                                                        adj.mod.out       = F,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                        w                 = NULL,
                                                        adj.mod.insplt    = F,
                                                        num.truc.obs      = 15,
                                                        min.node          = 15)

final.tree.estg.glm.modinnd.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                      tree.list         = seq.created.estg.glm.modinnd.mis.cv2$tree.list,             
                                                      type.var          = "cont",
                                                      seed              = a[job.number],
                                                      n.cv              = 5,
                                                      adj.mod.out       = F,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                      min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.mis.cv2[[1]],
                                                         test.data    = test.data.mis,
                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                         noise.var    = data.cont.cont$noise.var,
                                                         corr.split   = data.cont.cont$corr.split,
                                                         where.split  = data.cont.cont$where.split,
                                                         dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modinnd.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinnd.mis.cv2$corr.frst.splt <- as.character(seq.created.estg.glm.modinnd.mis.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("6")

performance.hetero.g <- list(glm.modinnd.true.cv2 = eval.final.estg.glm.modinnd.true.cv2,
                             glm.modinnd.nois.cv2 = eval.final.estg.glm.modinnd.nois.cv2,
                             glm.modinnd.mis.cv2  = eval.final.estg.glm.modinnd.mis.cv2)

#####################################################################################################################
####################### 7. dr: True GLM Model in node, True propensity score model in node, Cv2 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                    est.used          = "DR",
                                                                    type.var          = "cont",
                                                                    propsc.mod.out    = F,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                    propsc.mod.insplt = F,
                                                                    adj.mod.out       = F, 
                                                                    adj.mthd          = "GLM", 
                                                                    adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                    adj.mod.insplt    = F, 
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)

final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                   tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2$tree.list,
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

eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2[[1]],
                                                                     test.data    = data.cont.cont$test.data,
                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                     noise.var    = data.cont.cont$noise.var,
                                                                     corr.split   = data.cont.cont$corr.split,
                                                                     where.split  = data.cont.cont$where.split,
                                                                     dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("7")

#####################################################################################################################
#################### 8. dr: True GLM Model in node, Mis func propensity score model in node, Cv2 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                       est.used          = "DR",
                                                                       type.var          = "cont",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                       propsc.mod.insplt = F,
                                                                       adj.mod.out       = F, 
                                                                       adj.mthd          = "GLM", 
                                                                       adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                       adj.mod.insplt    = F, 
                                                                       num.truc.obs      = 30,
                                                                       min.node          = 20)

final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                      tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$tree.list,
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

eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2[[1]],
                                                                        test.data    = data.cont.cont$test.data,
                                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                                        noise.var    = data.cont.cont$noise.var,
                                                                        corr.split   = data.cont.cont$corr.split,
                                                                        where.split  = data.cont.cont$where.split,
                                                                        dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("8")

#####################################################################################################################
#################### 9. dr: Mis func GLM Model in node, True propensity score model in node, Cv2 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                       est.used          = "DR",
                                                                       type.var          = "cont",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                       propsc.mod.insplt = F,
                                                                       adj.mod.out       = F, 
                                                                       adj.mthd          = "GLM", 
                                                                       adj.form.true     = NULL, 
                                                                       adj.mod.insplt    = F, 
                                                                       num.truc.obs      = 30,
                                                                       min.node          = 20)

final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                      tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$tree.list,
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

eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2[[1]],
                                                                        test.data    = data.cont.cont$test.data,
                                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                                        noise.var    = data.cont.cont$noise.var,
                                                                        corr.split   = data.cont.cont$corr.split,
                                                                        where.split  = data.cont.cont$where.split,
                                                                        dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("9")

#####################################################################################################################
################# 10. dr: Mis func GLM Model in node, Mis func propensity score model in node, Cv2 ##################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                          est.used          = "DR",
                                                                          type.var          = "cont",
                                                                          propsc.mod.out    = F,
                                                                          propsc.mthd       = "GLM",
                                                                          propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                          propsc.mod.insplt = F,
                                                                          adj.mod.out       = F, 
                                                                          adj.mthd          = "GLM", 
                                                                          adj.form.true     = NULL, 
                                                                          adj.mod.insplt    = F, 
                                                                          num.truc.obs      = 30,
                                                                          min.node          = 20)

final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                         tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$tree.list,
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

eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2[[1]],
                                                                           test.data    = data.cont.cont$test.data,
                                                                           true.trt.eff = data.cont.cont$true.trt.eff,
                                                                           noise.var    = data.cont.cont$noise.var,
                                                                           corr.split   = data.cont.cont$corr.split,
                                                                           where.split  = data.cont.cont$where.split,
                                                                           dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("10")

#####################################################################################################################
########### 11. dr: Unmeasured cov GLM Model in node, Unmeasured cov propensity score model in node, Cv2 ############
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                                    est.used          = "DR",
                                                                    type.var          = "cont",
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

final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                                   tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2$tree.list,
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

eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv2[[1]],
                                                                     test.data    = test.data.mis,
                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                     noise.var    = data.cont.cont$noise.var,
                                                                     corr.split   = data.cont.cont$corr.split,
                                                                     where.split  = data.cont.cont$where.split,
                                                                     dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("11")

performance.hetero.drInnd <- list(adjTGlmInnd.propscTGlmInnd.cv2       = eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2,
                                  adjTGlmInnd.propscNoisGlmInnd.cv2    = eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2,
                                  adjNoisGlmInnd.propscTGlmInnd.cv2    = eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2,
                                  adjNoisGlmInnd.propscNoisGlmInnd.cv2 = eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2,
                                  adjFGlmInnd.propscFGlmInnd.cv2       = eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)



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
########################## 12. ipw: GLM Model, inside node, True propensity score model, Cv2 ########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                              est.used          = "IPW",
                                                              type.var          = "cont",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ X1 + X2 + X3",
                                                              w                 = NULL,
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.true.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont,
                                                              tree.list        = seq.created.estipw.glm.propscinnd.true.cv2$tree.list,
                                                              type.var         = "cont",
                                                              seed             = a[job.number],
                                                              n.cv             = 5,
                                                              propsc.mod.out   = F, 
                                                              propsc.mthd      = "GLM", 
                                                              propsc.form.true = "A ~ X1 + X2 + X3",
                                                              min.obs.mod      = 5)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.true.cv2[[1]],
                                                               test.data    = data.cont.cont$test.data,
                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                               noise.var    = data.cont.cont$noise.var,
                                                               corr.split   = data.cont.cont$corr.split,
                                                               where.split  = data.cont.cont$where.split,
                                                               dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinnd.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("12")

#####################################################################################################################
######################## 13. ipw: GLM Model, inside node, Mis func propensity score model, Cv2 ######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                              est.used          = "IPW",
                                                              type.var          = "cont",
                                                              propsc.mod.out    = F,
                                                              propsc.mthd       = "GLM",
                                                              propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                              w                 = NULL,
                                                              propsc.mod.insplt = F,
                                                              num.truc.obs      = 30,
                                                              min.node          = 20)

final.tree.estipw.glm.propscinnd.nois.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont,
                                                              tree.list        = seq.created.estipw.glm.propscinnd.nois.cv2$tree.list,
                                                              type.var         = "cont",
                                                              seed             = a[job.number],
                                                              n.cv             = 5,
                                                              propsc.mod.out   = F, 
                                                              propsc.mthd      = "GLM", 
                                                              propsc.form.true = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                              min.obs.mod      = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.nois.cv2[[1]],
                                                               test.data    = data.cont.cont$test.data,
                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                               noise.var    = data.cont.cont$noise.var,
                                                               corr.split   = data.cont.cont$corr.split,
                                                               where.split  = data.cont.cont$where.split,
                                                               dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinnd.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("13")

#####################################################################################################################
#################### 14. ipw: GLM Model, inside node, Unmeasured cov propensity score model, Cv2 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estipw.glm.propscinnd.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                             est.used          = "IPW",
                                                             type.var          = "cont",
                                                             propsc.mod.out    = F,
                                                             propsc.mthd       = "GLM",
                                                             propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                             w                 = NULL,
                                                             propsc.mod.insplt = F,
                                                             num.truc.obs      = 30,
                                                             min.node          = 20)

final.tree.estipw.glm.propscinnd.mis.cv2 <- EstIpw.CvMethod2(data.used        = data.used.full.cont.cont.mis,
                                                             tree.list        = seq.created.estipw.glm.propscinnd.mis.cv2$tree.list,
                                                             type.var         = "cont",
                                                             seed             = a[job.number],
                                                             n.cv             = 5,
                                                             propsc.mod.out   = F, 
                                                             propsc.mthd      = "GLM", 
                                                             propsc.form.true = "A ~ X1 + X3 + X4 + X5 + X6",
                                                             min.obs.mod      = 10)
t1 <- Sys.time()

eval.final.estipw.glm.propscinnd.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estipw.glm.propscinnd.mis.cv2[[1]],
                                                              test.data    = test.data.mis,
                                                              true.trt.eff = data.cont.cont$true.trt.eff,
                                                              noise.var    = data.cont.cont$noise.var,
                                                              corr.split   = data.cont.cont$corr.split,
                                                              where.split  = data.cont.cont$where.split,
                                                              dir.split    = data.cont.cont$dir.split)
eval.final.estipw.glm.propscinnd.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("14")

performance.homo.ipw <- list(glm.propscinnd.true.cv2 = eval.final.estipw.glm.propscinnd.true.cv2,
                             glm.propscinnd.nois.cv2 = eval.final.estipw.glm.propscinnd.nois.cv2,
                             glm.propscinnd.mis.cv2  = eval.final.estipw.glm.propscinnd.mis.cv2)

#####################################################################################################################
############################# 15. g: GLM Model, inside node, True adjustment model, Cv2 #############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                         est.used          = "G",
                                                         type.var          = "cont",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 10,
                                                         min.node          = 10)

final.tree.estg.glm.modinnd.true.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                       tree.list         = seq.created.estg.glm.modinnd.true.cv2$tree.list,             
                                                       type.var          = "cont",
                                                       seed              = a[job.number],
                                                       n.cv              = 5,
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                       min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.true.cv2[[1]],
                                                          test.data    = data.cont.cont$test.data,
                                                          true.trt.eff = data.cont.cont$true.trt.eff,
                                                          noise.var    = data.cont.cont$noise.var,
                                                          corr.split   = data.cont.cont$corr.split,
                                                          where.split  = data.cont.cont$where.split,
                                                          dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modinnd.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("15")

#####################################################################################################################
########################### 16. g: GLM Model, inside node, Mis func adjustment model, Cv2 ###########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                         est.used          = "G",
                                                         type.var          = "cont",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = NULL,
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 15,
                                                         min.node          = 15)

final.tree.estg.glm.modinnd.nois.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                       tree.list         = seq.created.estg.glm.modinnd.nois.cv2$tree.list,             
                                                       type.var          = "cont",
                                                       seed              = a[job.number],
                                                       n.cv              = 5,
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = NULL,
                                                       min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.nois.cv2[[1]],
                                                          test.data    = data.cont.cont$test.data,
                                                          true.trt.eff = data.cont.cont$true.trt.eff,
                                                          noise.var    = data.cont.cont$noise.var,
                                                          corr.split   = data.cont.cont$corr.split,
                                                          where.split  = data.cont.cont$where.split,
                                                          dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modinnd.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("16")

#####################################################################################################################
######################### 17. g: GLM Model, inside node, Unmeasured cov adjustment model, Cv2 #######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                        est.used          = "G",
                                                        type.var          = "cont",
                                                        adj.mod.out       = F,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                        w                 = NULL,
                                                        adj.mod.insplt    = F,
                                                        num.truc.obs      = 15,
                                                        min.node          = 15)

final.tree.estg.glm.modinnd.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                      tree.list         = seq.created.estg.glm.modinnd.mis.cv2$tree.list,             
                                                      type.var          = "cont",
                                                      seed              = a[job.number],
                                                      n.cv              = 5,
                                                      adj.mod.out       = F,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                      min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.mis.cv2[[1]],
                                                         test.data    = test.data.mis,
                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                         noise.var    = data.cont.cont$noise.var,
                                                         corr.split   = data.cont.cont$corr.split,
                                                         where.split  = data.cont.cont$where.split,
                                                         dir.split    = data.cont.cont$dir.split)

eval.final.estg.glm.modinnd.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("17")

performance.homo.g <- list(glm.modinnd.true.cv2 = eval.final.estg.glm.modinnd.true.cv2,
                           glm.modinnd.nois.cv2 = eval.final.estg.glm.modinnd.nois.cv2,
                           glm.modinnd.mis.cv2  = eval.final.estg.glm.modinnd.mis.cv2)

#####################################################################################################################
###################### 18. dr: True GLM Model in node, True propensity score model in node, Cv2 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                    est.used          = "DR",
                                                                    type.var          = "cont",
                                                                    propsc.mod.out    = F,
                                                                    propsc.mthd       = "GLM",
                                                                    propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                    propsc.mod.insplt = F,
                                                                    adj.mod.out       = F, 
                                                                    adj.mthd          = "GLM", 
                                                                    adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                    adj.mod.insplt    = F, 
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)

final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                   tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2$tree.list,
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

eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2[[1]],
                                                                     test.data    = data.cont.cont$test.data,
                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                     noise.var    = data.cont.cont$noise.var,
                                                                     corr.split   = data.cont.cont$corr.split,
                                                                     where.split  = data.cont.cont$where.split,
                                                                     dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("18")

#####################################################################################################################
################### 19. dr: True GLM Model in node, Mis func propensity score model in node, Cv2 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                       est.used          = "DR",
                                                                       type.var          = "cont",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                       propsc.mod.insplt = F,
                                                                       adj.mod.out       = F, 
                                                                       adj.mthd          = "GLM", 
                                                                       adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                       adj.mod.insplt    = F, 
                                                                       num.truc.obs      = 30,
                                                                       min.node          = 20)

final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                      tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$tree.list,
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

eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2[[1]],
                                                                        test.data    = data.cont.cont$test.data,
                                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                                        noise.var    = data.cont.cont$noise.var,
                                                                        corr.split   = data.cont.cont$corr.split,
                                                                        where.split  = data.cont.cont$where.split,
                                                                        dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("19")

#####################################################################################################################
################## 20. dr: Mis func GLM Model in node, True propensity score model in node, Cv2 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                       est.used          = "DR",
                                                                       type.var          = "cont",
                                                                       propsc.mod.out    = F,
                                                                       propsc.mthd       = "GLM",
                                                                       propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                       propsc.mod.insplt = F,
                                                                       adj.mod.out       = F, 
                                                                       adj.mthd          = "GLM", 
                                                                       adj.form.true     = NULL, 
                                                                       adj.mod.insplt    = F, 
                                                                       num.truc.obs      = 30,
                                                                       min.node          = 20)

final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                      tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$tree.list,
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

eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2[[1]],
                                                                        test.data    = data.cont.cont$test.data,
                                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                                        noise.var    = data.cont.cont$noise.var,
                                                                        corr.split   = data.cont.cont$corr.split,
                                                                        where.split  = data.cont.cont$where.split,
                                                                        dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("20")

#####################################################################################################################
################ 21. dr: Mis func GLM Model in node, Mis func propensity score model in node, Cv2 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                          est.used          = "DR",
                                                                          type.var          = "cont",
                                                                          propsc.mod.out    = F,
                                                                          propsc.mthd       = "GLM",
                                                                          propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                          propsc.mod.insplt = F,
                                                                          adj.mod.out       = F, 
                                                                          adj.mthd          = "GLM", 
                                                                          adj.form.true     = NULL, 
                                                                          adj.mod.insplt    = F, 
                                                                          num.truc.obs      = 30,
                                                                          min.node          = 20)

final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                         tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$tree.list,
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

eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2[[1]],
                                                                           test.data    = data.cont.cont$test.data,
                                                                           true.trt.eff = data.cont.cont$true.trt.eff,
                                                                           noise.var    = data.cont.cont$noise.var,
                                                                           corr.split   = data.cont.cont$corr.split,
                                                                           where.split  = data.cont.cont$where.split,
                                                                           dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("21")

#####################################################################################################################
############ 22. dr: Unmeasured cov GLM Model in node, Unmeasured cov propensity score model in node, Cv2 ###########
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                                    est.used          = "DR",
                                                                    type.var          = "cont",
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

final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                                   tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2$tree.list,
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

eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv2[[1]],
                                                                     test.data    = test.data.mis,
                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                     noise.var    = data.cont.cont$noise.var,
                                                                     corr.split   = data.cont.cont$corr.split,
                                                                     where.split  = data.cont.cont$where.split,
                                                                     dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("22")

performance.homo.drInnd <- list(adjTGlmInnd.propscTGlmInnd.cv2       = eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2,
                                adjTGlmInnd.propscNoisGlmInnd.cv2    = eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2,
                                adjNoisGlmInnd.propscTGlmInnd.cv2    = eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2,
                                adjNoisGlmInnd.propscNoisGlmInnd.cv2 = eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2,
                                adjFGlmInnd.propscFGlmInnd.cv2       = eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)

file.name = paste("../Data/AppendixC4/", toString(job.number), ".RData", sep = "")
save(performance.hetero.ipw, performance.homo.ipw,
     performance.hetero.g, performance.homo.g,
     performance.hetero.drInnd, performance.homo.drInnd,
     file = file.name)

