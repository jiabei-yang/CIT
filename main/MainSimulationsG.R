#!/usr/bin/env Rscript
setwd("..")
folder <- paste(getwd(), "/Functions/", sep="")
functions <- list.files(folder)
functions <- paste(folder, functions, sep = "")
for (i in functions){
  source(i)
}

setwd("main/")

# Read in the arguments from command line
option_list = list(
  make_option(c("-a", "--array"), type="integer"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
iters <- opt$array

job.number <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) + iters
load("../seed1000.rda")
set.seed(a[job.number])
print(job.number)

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
############################## 1. g: GLM Model, inside node, True adjustment model, cv1 #############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.true.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                         est.used          = "G",
                                                         type.var          = "cont",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 10,
                                                         min.node          = 10)

final.tree.estg.glm.modinnd.true.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                       tree.list         = seq.created.estg.glm.modinnd.true.cv1$tree.list,
                                                       lambda.used       = qchisq(0.95, 1),
                                                       val.sample        = data.validation.cont.cont,              
                                                       type.var          = "cont",
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                       val.w             = NULL,
                                                       adj.mod.insplt    = F,
                                                       min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.true.cv1[[1]],
                                                          test.data    = data.cont.cont$test.data,
                                                          true.trt.eff = data.cont.cont$true.trt.eff,
                                                          noise.var    = data.cont.cont$noise.var,
                                                          corr.split   = data.cont.cont$corr.split,
                                                          where.split  = data.cont.cont$where.split,
                                                          dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modinnd.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinnd.true.cv1$corr.frst.splt <- as.character(seq.created.estg.glm.modinnd.true.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("1")

#####################################################################################################################
############################ 2. g: GLM Model, inside node, Mis func adjustment model, cv1 ###########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.nois.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                         est.used          = "G",
                                                         type.var          = "cont",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = NULL,
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 15,
                                                         min.node          = 15)

final.tree.estg.glm.modinnd.nois.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                       tree.list         = seq.created.estg.glm.modinnd.nois.cv1$tree.list,
                                                       lambda.used       = qchisq(0.95, 1),
                                                       val.sample        = data.validation.cont.cont,              
                                                       type.var          = "cont",
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = NULL,
                                                       val.w             = NULL,
                                                       adj.mod.insplt    = F,
                                                       min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.nois.cv1[[1]],
                                                          test.data    = data.cont.cont$test.data,
                                                          true.trt.eff = data.cont.cont$true.trt.eff,
                                                          noise.var    = data.cont.cont$noise.var,
                                                          corr.split   = data.cont.cont$corr.split,
                                                          where.split  = data.cont.cont$where.split,
                                                          dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modinnd.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinnd.nois.cv1$corr.frst.splt <- as.character(seq.created.estg.glm.modinnd.nois.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("2")

#####################################################################################################################
######################## 3. g: GLM Model, inside node, Unmeasured cov adjustment model, cv1 #########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                        est.used          = "G",
                                                        type.var          = "cont",
                                                        adj.mod.out       = F,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                        w                 = NULL,
                                                        adj.mod.insplt    = F,
                                                        num.truc.obs      = 15,
                                                        min.node          = 15)

final.tree.estg.glm.modinnd.mis.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                      tree.list         = seq.created.estg.glm.modinnd.mis.cv1$tree.list,
                                                      lambda.used       = qchisq(0.95, 1),
                                                      val.sample        = data.validation.cont.cont.mis,              
                                                      type.var          = "cont",
                                                      adj.mod.out       = F,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                      val.w             = NULL,
                                                      adj.mod.insplt    = F,
                                                      min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.mis.cv1[[1]],
                                                         test.data    = test.data.mis,
                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                         noise.var    = data.cont.cont$noise.var,
                                                         corr.split   = data.cont.cont$corr.split,
                                                         where.split  = data.cont.cont$where.split,
                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modinnd.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
eval.final.estg.glm.modinnd.mis.cv1$corr.frst.splt <- as.character(seq.created.estg.glm.modinnd.mis.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
print("3")

performance.hetero.g <- list(glm.modinnd.true.cv1 = eval.final.estg.glm.modinnd.true.cv1,
                             glm.modinnd.nois.cv1 = eval.final.estg.glm.modinnd.nois.cv1,
                             glm.modinnd.mis.cv1  = eval.final.estg.glm.modinnd.mis.cv1)



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
############################## 4. g: GLM Model, inside node, True adjustment model, cv1 #############################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.true.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                         est.used          = "G",
                                                         type.var          = "cont",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 10,
                                                         min.node          = 10)

final.tree.estg.glm.modinnd.true.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                       tree.list         = seq.created.estg.glm.modinnd.true.cv1$tree.list,
                                                       lambda.used       = qchisq(0.95, 1),
                                                       val.sample        = data.validation.cont.cont,              
                                                       type.var          = "cont",
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                       val.w             = NULL,
                                                       adj.mod.insplt    = F,
                                                       min.obs.mod       = 10)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.true.cv1[[1]],
                                                          test.data    = data.cont.cont$test.data,
                                                          true.trt.eff = data.cont.cont$true.trt.eff,
                                                          noise.var    = data.cont.cont$noise.var,
                                                          corr.split   = data.cont.cont$corr.split,
                                                          where.split  = data.cont.cont$where.split,
                                                          dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modinnd.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("4")

#####################################################################################################################
############################# 5. g: GLM Model, inside node, Mis func adjustment model, cv1 ##########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.nois.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                         est.used          = "G",
                                                         type.var          = "cont",
                                                         adj.mod.out       = F,
                                                         adj.mthd          = "GLM",
                                                         adj.form.true     = NULL,
                                                         w                 = NULL,
                                                         adj.mod.insplt    = F,
                                                         num.truc.obs      = 15,
                                                         min.node          = 15)

final.tree.estg.glm.modinnd.nois.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                       tree.list         = seq.created.estg.glm.modinnd.nois.cv1$tree.list,
                                                       lambda.used       = qchisq(0.95, 1),
                                                       val.sample        = data.validation.cont.cont,              
                                                       type.var          = "cont",
                                                       adj.mod.out       = F,
                                                       adj.mthd          = "GLM",
                                                       adj.form.true     = NULL,
                                                       val.w             = NULL,
                                                       adj.mod.insplt    = F,
                                                       min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.nois.cv1[[1]],
                                                          test.data    = data.cont.cont$test.data,
                                                          true.trt.eff = data.cont.cont$true.trt.eff,
                                                          noise.var    = data.cont.cont$noise.var,
                                                          corr.split   = data.cont.cont$corr.split,
                                                          where.split  = data.cont.cont$where.split,
                                                          dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modinnd.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("5")

#####################################################################################################################
######################### 6. g: GLM Model, inside node, Unmeasured cov adjustment model, cv1 ########################
#####################################################################################################################
t0 <- Sys.time()
seq.created.estg.glm.modinnd.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                        est.used          = "G",
                                                        type.var          = "cont",
                                                        adj.mod.out       = F,
                                                        adj.mthd          = "GLM",
                                                        adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                        w                 = NULL,
                                                        adj.mod.insplt    = F,
                                                        num.truc.obs      = 15,
                                                        min.node          = 15)

final.tree.estg.glm.modinnd.mis.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                      tree.list         = seq.created.estg.glm.modinnd.mis.cv1$tree.list,
                                                      lambda.used       = qchisq(0.95, 1),
                                                      val.sample        = data.validation.cont.cont.mis,              
                                                      type.var          = "cont",
                                                      adj.mod.out       = F,
                                                      adj.mthd          = "GLM",
                                                      adj.form.true     = "Y ~ A + X1 + X3 + X4 + X5 + X6 + A:X1 + A:X3 + A:X4 + A:X5 + A:X6",
                                                      val.w             = NULL,
                                                      adj.mod.insplt    = F,
                                                      min.obs.mod       = 15)
t1 <- Sys.time()

eval.final.estg.glm.modinnd.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.mis.cv1[[1]],
                                                         test.data    = test.data.mis,
                                                         true.trt.eff = data.cont.cont$true.trt.eff,
                                                         noise.var    = data.cont.cont$noise.var,
                                                         corr.split   = data.cont.cont$corr.split,
                                                         where.split  = data.cont.cont$where.split,
                                                         dir.split    = data.cont.cont$dir.split)
eval.final.estg.glm.modinnd.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
print("6")

performance.homo.g <- list(glm.modinnd.true.cv1 = eval.final.estg.glm.modinnd.true.cv1,
                           glm.modinnd.nois.cv1 = eval.final.estg.glm.modinnd.nois.cv1,
                           glm.modinnd.mis.cv1  = eval.final.estg.glm.modinnd.mis.cv1)

file.name = paste("../Data/MainSimulationsG/", toString(job.number), ".RData", sep = "")
save(performance.hetero.g, performance.homo.g,
     file = file.name)

