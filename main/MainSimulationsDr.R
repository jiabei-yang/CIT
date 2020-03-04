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
###################### 1. dr: True GLM Model in node, True propensity score model in node, Cv1 ######################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
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

final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                   tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                   lambda.used       = qchisq(0.95, 1),
                                                                   val.sample        = data.validation.cont.cont,
                                                                   type.var          = "cont",
                                                                   propsc.mod.out    = F,
                                                                   propsc.mthd       = "GLM",
                                                                   propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                   propsc.mod.insplt = F,
                                                                   adj.mod.out       = F,
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                   adj.mod.insplt    = F,
                                                                   min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1[[1]],
                                                                     test.data    = data.cont.cont$test.data,
                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                     noise.var    = data.cont.cont$noise.var,
                                                                     corr.split   = data.cont.cont$corr.split,
                                                                     where.split  = data.cont.cont$where.split,
                                                                     dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("1")

#####################################################################################################################
#################### 2. dr: True GLM Model in node, Mis func propensity score model in node, Cv1 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
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

final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                      tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                      lambda.used       = qchisq(0.95, 1),
                                                                      val.sample        = data.validation.cont.cont,
                                                                      type.var          = "cont",
                                                                      propsc.mod.out    = F,
                                                                      propsc.mthd       = "GLM",
                                                                      propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                      propsc.mod.insplt = F,
                                                                      adj.mod.out       = F,
                                                                      adj.mthd          = "GLM",
                                                                      adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                      adj.mod.insplt    = F,
                                                                      min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1[[1]],
                                                                        test.data    = data.cont.cont$test.data,
                                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                                        noise.var    = data.cont.cont$noise.var,
                                                                        corr.split   = data.cont.cont$corr.split,
                                                                        where.split  = data.cont.cont$where.split,
                                                                        dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("2")

#####################################################################################################################
##################### 3. dr: Mis func GLM Model in node, True propensity score model in node, Cv1 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
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

final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                      tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                      lambda.used       = qchisq(0.95, 1),
                                                                      val.sample        = data.validation.cont.cont,
                                                                      type.var          = "cont",
                                                                      propsc.mod.out    = F,
                                                                      propsc.mthd       = "GLM",
                                                                      propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                      propsc.mod.insplt = F,
                                                                      adj.mod.out       = F,
                                                                      adj.mthd          = "GLM",
                                                                      adj.form.true     = NULL, 
                                                                      adj.mod.insplt    = F,
                                                                      min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1[[1]],
                                                                        test.data    = data.cont.cont$test.data,
                                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                                        noise.var    = data.cont.cont$noise.var,
                                                                        corr.split   = data.cont.cont$corr.split,
                                                                        where.split  = data.cont.cont$where.split,
                                                                        dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("3")

#####################################################################################################################
################## 4. dr: Mis func GLM Model in node, Mis func propensity score model in node, Cv1 ##################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
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

final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                         tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                         lambda.used       = qchisq(0.95, 1),
                                                                         val.sample        = data.validation.cont.cont,
                                                                         type.var          = "cont",
                                                                         propsc.mod.out    = F,
                                                                         propsc.mthd       = "GLM",
                                                                         propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                         propsc.mod.insplt = F,
                                                                         adj.mod.out       = F,
                                                                         adj.mthd          = "GLM",
                                                                         adj.form.true     = NULL, 
                                                                         adj.mod.insplt    = F,
                                                                         min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1[[1]],
                                                                           test.data    = data.cont.cont$test.data,
                                                                           true.trt.eff = data.cont.cont$true.trt.eff,
                                                                           noise.var    = data.cont.cont$noise.var,
                                                                           corr.split   = data.cont.cont$corr.split,
                                                                           where.split  = data.cont.cont$where.split,
                                                                           dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("4")

#####################################################################################################################
############ 5. dr: Unmeasured cov GLM Model in node, Unmeasured cov propensity score model in node, Cv1 ############
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
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

final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                                   tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1$tree.list, 
                                                                   lambda.used       = qchisq(0.95, 1),
                                                                   val.sample        = data.validation.cont.cont.mis,
                                                                   type.var          = "cont",
                                                                   propsc.mod.out    = F,
                                                                   propsc.mthd       = "GLM",
                                                                   propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                   propsc.mod.insplt = F,
                                                                   adj.mod.out       = F,
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                   adj.mod.insplt    = F,
                                                                   min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1[[1]],
                                                                     test.data    = test.data.mis,
                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                     noise.var    = data.cont.cont$noise.var,
                                                                     corr.split   = data.cont.cont$corr.split,
                                                                     where.split  = data.cont.cont$where.split,
                                                                     dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("5")

performance.hetero.drInnd <- list(adjTGlmInnd.propscTGlmInnd.cv1       = eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1,
                                  adjTGlmInnd.propscNoisGlmInnd.cv1    = eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1,
                                  adjNoisGlmInnd.propscTGlmInnd.cv1    = eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1,
                                  adjNoisGlmInnd.propscNoisGlmInnd.cv1 = eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1,
                                  adjFGlmInnd.propscFGlmInnd.cv1       = eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1)

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
###################### 6. dr: True GLM Model in node, True propensity score model in node, Cv1 #####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
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

final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                   tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                   lambda.used       = qchisq(0.95, 1),
                                                                   val.sample        = data.validation.cont.cont,
                                                                   type.var          = "cont",
                                                                   propsc.mod.out    = F,
                                                                   propsc.mthd       = "GLM",
                                                                   propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                   propsc.mod.insplt = F,
                                                                   adj.mod.out       = F,
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                   adj.mod.insplt    = F,
                                                                   min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1[[1]],
                                                                     test.data    = data.cont.cont$test.data,
                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                     noise.var    = data.cont.cont$noise.var,
                                                                     corr.split   = data.cont.cont$corr.split,
                                                                     where.split  = data.cont.cont$where.split,
                                                                     dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("6")

#####################################################################################################################
##################### 7. dr: True GLM Model in node, Mis func propensity score model in node, Cv1 ###################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
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

final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                      tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                      lambda.used       = qchisq(0.95, 1),
                                                                      val.sample        = data.validation.cont.cont,
                                                                      type.var          = "cont",
                                                                      propsc.mod.out    = F,
                                                                      propsc.mthd       = "GLM",
                                                                      propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                      propsc.mod.insplt = F,
                                                                      adj.mod.out       = F,
                                                                      adj.mthd          = "GLM",
                                                                      adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                      adj.mod.insplt    = F,
                                                                      min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1[[1]],
                                                                        test.data    = data.cont.cont$test.data,
                                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                                        noise.var    = data.cont.cont$noise.var,
                                                                        corr.split   = data.cont.cont$corr.split,
                                                                        where.split  = data.cont.cont$where.split,
                                                                        dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("7")

#####################################################################################################################
#################### 8. dr: Mis func GLM Model in node, True propensity score model in node, Cv1 ####################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
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

final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                      tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                      lambda.used       = qchisq(0.95, 1),
                                                                      val.sample        = data.validation.cont.cont,
                                                                      type.var          = "cont",
                                                                      propsc.mod.out    = F,
                                                                      propsc.mthd       = "GLM",
                                                                      propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                      propsc.mod.insplt = F,
                                                                      adj.mod.out       = F,
                                                                      adj.mthd          = "GLM",
                                                                      adj.form.true     = NULL, 
                                                                      adj.mod.insplt    = F,
                                                                      min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1[[1]],
                                                                        test.data    = data.cont.cont$test.data,
                                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                                        noise.var    = data.cont.cont$noise.var,
                                                                        corr.split   = data.cont.cont$corr.split,
                                                                        where.split  = data.cont.cont$where.split,
                                                                        dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("8")

#####################################################################################################################
################## 9. dr: Mis func GLM Model in node, Mis func propensity score model in node, Cv1 ##################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
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

final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                         tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                         lambda.used       = qchisq(0.95, 1),
                                                                         val.sample        = data.validation.cont.cont,
                                                                         type.var          = "cont",
                                                                         propsc.mod.out    = F,
                                                                         propsc.mthd       = "GLM",
                                                                         propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                         propsc.mod.insplt = F,
                                                                         adj.mod.out       = F,
                                                                         adj.mthd          = "GLM",
                                                                         adj.form.true     = NULL, 
                                                                         adj.mod.insplt    = F,
                                                                         min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1[[1]],
                                                                           test.data    = data.cont.cont$test.data,
                                                                           true.trt.eff = data.cont.cont$true.trt.eff,
                                                                           noise.var    = data.cont.cont$noise.var,
                                                                           corr.split   = data.cont.cont$corr.split,
                                                                           where.split  = data.cont.cont$where.split,
                                                                           dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("9")

#####################################################################################################################
########### 10. dr: Unmeasured cov GLM Model in node, Unmeasured cov propensity score model in node, Cv1 ############
#####################################################################################################################
t0 <- Sys.time()
seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
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

final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                                   tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1$tree.list, 
                                                                   lambda.used       = qchisq(0.95, 1),
                                                                   val.sample        = data.validation.cont.cont.mis,
                                                                   type.var          = "cont",
                                                                   propsc.mod.out    = F,
                                                                   propsc.mthd       = "GLM",
                                                                   propsc.form.true  = "A ~ X1 + X3 + X4 + X5 + X6",
                                                                   propsc.mod.insplt = F,
                                                                   adj.mod.out       = F,
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6", 
                                                                   adj.mod.insplt    = F,
                                                                   min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1[[1]],
                                                                     test.data    = test.data.mis,
                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                     noise.var    = data.cont.cont$noise.var,
                                                                     corr.split   = data.cont.cont$corr.split,
                                                                     where.split  = data.cont.cont$where.split,
                                                                     dir.split    = data.cont.cont$dir.split)
eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("10")

performance.homo.drInnd <- list(adjTGlmInnd.propscTGlmInnd.cv1       = eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1,
                                adjTGlmInnd.propscNoisGlmInnd.cv1    = eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1,
                                adjNoisGlmInnd.propscTGlmInnd.cv1    = eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1,
                                adjNoisGlmInnd.propscNoisGlmInnd.cv1 = eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1,
                                adjFGlmInnd.propscFGlmInnd.cv1       = eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1)

