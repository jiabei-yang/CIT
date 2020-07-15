#!/usr/bin/env Rscript
setwd("../")
folder <- paste(getwd(), "/Functions/", sep="")
functions <- list.files(folder)
functions <- grep("main.R", functions, invert = T, value = T)
functions <- paste(folder, functions, sep = "")
for (i in functions){
  source(i)
}

setwd("Appendix/")

registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))

# Read in the arguments from command line
option_list = list(
  make_option(c("-s", "--start"), type="integer"),
  make_option(c("-e", "--end"), type = "integer"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
start <- opt$start
end   <- opt$end

load("../seed1000.rda")

performance.gOut <- 
  foreach(i = start:end,
          .combine  = "rbind") %dopar% {
            
            set.seed(a[i])
            
            #####################################################################################################################
            ################################################# Heterogeneous #####################################################
            #####################################################################################################################
            data.cont.cont            <- makeData.cont.eff.cont(1000, 1000, coeff.prop.sc = 0.6)
            data.used.full.cont.cont  <- data.cont.cont$data.used
            data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
            # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
            data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ]  
            
            #####################################################################################################################
            ############################ 1. g: GLM Model, outside node, True adjustment model, Cv1 ##############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modout.true.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                    est.used          = "G",
                                                                    type.var          = "cont",
                                                                    adj.mod.loc       = "out",
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)
            
            final.tree.estg.glm.modout.true.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                                  tree.list         = seq.created.estg.glm.modout.true.cv1$tree.list,
                                                                  lambda.used       = qchisq(0.95, 1),
                                                                  val.sample        = data.validation.cont.cont,              
                                                                  type.var          = "cont",
                                                                  adj.mod.loc       = "out",
                                                                  adj.mthd          = "GLM",
                                                                  adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)")
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
            eval.final.estg.glm.modout.true.cv1                <- unlist(eval.final.estg.glm.modout.true.cv1)
            names(eval.final.estg.glm.modout.true.cv1)         <- paste0("hetero.g.glm.modout.true.cv1.", names(eval.final.estg.glm.modout.true.cv1))
            print("1.cv1.T")
            
            #####################################################################################################################
            ############################# 2. g: GLM Model, outside node, Noisy adjustment model, cv1 ############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modout.nois.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                    est.used          = "G",
                                                                    type.var          = "cont",
                                                                    adj.mod.loc       = "out",
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = NULL,
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)
            
            final.tree.estg.glm.modout.nois.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                                  tree.list         = seq.created.estg.glm.modout.nois.cv1$tree.list,
                                                                  lambda.used       = qchisq(0.95, 1),
                                                                  val.sample        = data.validation.cont.cont,              
                                                                  type.var          = "cont",
                                                                  adj.mod.loc       = "out",
                                                                  adj.mthd          = "GLM",
                                                                  adj.form.true     = NULL)
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
            eval.final.estg.glm.modout.nois.cv1                <- unlist(eval.final.estg.glm.modout.nois.cv1)
            names(eval.final.estg.glm.modout.nois.cv1)         <- paste0("hetero.g.glm.modout.nois.cv1.", names(eval.final.estg.glm.modout.nois.cv1))
            print("2.cv1.nois")
            
            #####################################################################################################################
            ########################## 3. g: GLM Model, outside node, Misspecified adjustment model, cv1 ########################
            #####################################################################################################################
            data.used.cont.cont.mis <- data.used.cont.cont %>%
              select(-X2)
            data.validation.cont.cont.mis <- data.validation.cont.cont %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estg.glm.modout.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                                   est.used          = "G",
                                                                   type.var          = "cont",
                                                                   adj.mod.loc       = "out",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = NULL,
                                                                   num.truc.obs      = 30,
                                                                   min.node          = 20)
            
            final.tree.estg.glm.modout.mis.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                                 tree.list         = seq.created.estg.glm.modout.mis.cv1$tree.list,
                                                                 lambda.used       = qchisq(0.95, 1),
                                                                 val.sample        = data.validation.cont.cont.mis,              
                                                                 type.var          = "cont",
                                                                 adj.mod.loc       = "out",
                                                                 adj.mthd          = "GLM",
                                                                 adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modout.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.mis.cv1[[1]],
                                                                    test.data    = data.cont.cont$test.data,
                                                                    true.trt.eff = data.cont.cont$true.trt.eff,
                                                                    noise.var    = data.cont.cont$noise.var,
                                                                    corr.split   = data.cont.cont$corr.split,
                                                                    where.split  = data.cont.cont$where.split,
                                                                    dir.split    = data.cont.cont$dir.split)
            eval.final.estg.glm.modout.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modout.mis.cv1$corr.frst.splt <- as.character(seq.created.estg.glm.modout.mis.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
            eval.final.estg.glm.modout.mis.cv1                <- unlist(eval.final.estg.glm.modout.mis.cv1)
            names(eval.final.estg.glm.modout.mis.cv1)         <- paste0("hetero.g.glm.modout.mis.cv1.", names(eval.final.estg.glm.modout.mis.cv1))
            print("3.cv1.mis")
            
            #####################################################################################################################
            ############################ 4. g: GLM Model, outside node, True adjustment model, Cv2 ##############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modout.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                    est.used          = "G",
                                                                    type.var          = "cont",
                                                                    adj.mod.loc       = "out",
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)
            
            final.tree.estg.glm.modout.true.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                  tree.list         = seq.created.estg.glm.modout.true.cv2$tree.list,             
                                                                  type.var          = "cont",
                                                                  seed              = a[i],
                                                                  n.cv              = 5,
                                                                  adj.mod.loc       = "out",
                                                                  adj.mthd          = "GLM",
                                                                  adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)")
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
            eval.final.estg.glm.modout.true.cv2                <- unlist(eval.final.estg.glm.modout.true.cv2)
            names(eval.final.estg.glm.modout.true.cv2)         <- paste0("hetero.g.glm.modout.true.cv2.", names(eval.final.estg.glm.modout.true.cv2))
            print("4.cv2.T")
            
            #####################################################################################################################
            ############################# 5. g: GLM Model, outside node, Noisy adjustment model, Cv2 ############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modout.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                    est.used          = "G",
                                                                    type.var          = "cont",
                                                                    adj.mod.loc       = "out",
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = NULL,
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)
            
            final.tree.estg.glm.modout.nois.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                  tree.list         = seq.created.estg.glm.modout.nois.cv2$tree.list,             
                                                                  type.var          = "cont",
                                                                  seed              = a[i],
                                                                  n.cv              = 5,
                                                                  adj.mod.loc       = "out",
                                                                  adj.mthd          = "GLM",
                                                                  adj.form.true     = NULL)
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
            eval.final.estg.glm.modout.nois.cv2                <- unlist(eval.final.estg.glm.modout.nois.cv2)
            names(eval.final.estg.glm.modout.nois.cv2)         <- paste0("hetero.g.glm.modout.nois.cv2.", names(eval.final.estg.glm.modout.nois.cv2))
            print("5.cv2.nois")
            
            #####################################################################################################################
            ########################## 6. g: GLM Model, outside node, Misspecified adjustment model, Cv2 ########################
            #####################################################################################################################
            data.used.full.cont.cont.mis <- data.used.full.cont.cont %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estg.glm.modout.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                                   est.used          = "G",
                                                                   type.var          = "cont",
                                                                   adj.mod.loc       = "out",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = NULL,
                                                                   num.truc.obs      = 30,
                                                                   min.node          = 20)
            
            final.tree.estg.glm.modout.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                                 tree.list         = seq.created.estg.glm.modout.mis.cv2$tree.list,             
                                                                 type.var          = "cont",
                                                                 seed              = a[i],
                                                                 n.cv              = 5,
                                                                 adj.mod.loc       = "out",
                                                                 adj.mthd          = "GLM",
                                                                 adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modout.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.mis.cv2[[1]],
                                                                    test.data    = data.cont.cont$test.data,
                                                                    true.trt.eff = data.cont.cont$true.trt.eff,
                                                                    noise.var    = data.cont.cont$noise.var,
                                                                    corr.split   = data.cont.cont$corr.split,
                                                                    where.split  = data.cont.cont$where.split,
                                                                    dir.split    = data.cont.cont$dir.split)
            
            eval.final.estg.glm.modout.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modout.mis.cv2$corr.frst.splt <- as.character(seq.created.estg.glm.modout.mis.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
            eval.final.estg.glm.modout.mis.cv2                <- unlist(eval.final.estg.glm.modout.mis.cv2)
            names(eval.final.estg.glm.modout.mis.cv2)         <- paste0("hetero.g.glm.modout.mis.cv2.", names(eval.final.estg.glm.modout.mis.cv2))
            print("6.cv2.mis")
            
            performance.hetero.g <- c(eval.final.estg.glm.modout.true.cv1,
                                      eval.final.estg.glm.modout.nois.cv1,
                                      eval.final.estg.glm.modout.mis.cv1,
                                      eval.final.estg.glm.modout.true.cv2,
                                      eval.final.estg.glm.modout.nois.cv2,
                                      eval.final.estg.glm.modout.mis.cv2)
            
            #####################################################################################################################
            ################################################### Homogeneous #####################################################
            #####################################################################################################################
            data.cont.cont            <- makeData.cont.noeff.cont(1000, 1000, coeff.prop.sc = 0.6)
            data.used.full.cont.cont  <- data.cont.cont$data.used
            data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
            # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
            data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ]  
            
            #####################################################################################################################
            ############################ 7. g: GLM Model, outside node, True adjustment model, Cv1 ##############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modout.true.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                    est.used          = "G",
                                                                    type.var          = "cont",
                                                                    adj.mod.loc       = "out",
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)
            
            final.tree.estg.glm.modout.true.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                                  tree.list         = seq.created.estg.glm.modout.true.cv1$tree.list,
                                                                  lambda.used       = qchisq(0.95, 1),
                                                                  val.sample        = data.validation.cont.cont,              
                                                                  type.var          = "cont",
                                                                  adj.mod.loc       = "out",
                                                                  adj.mthd          = "GLM",
                                                                  adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)")
            t1 <- Sys.time()
            
            eval.final.estg.glm.modout.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.true.cv1[[1]],
                                                                     test.data    = data.cont.cont$test.data,
                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                     noise.var    = data.cont.cont$noise.var,
                                                                     corr.split   = data.cont.cont$corr.split,
                                                                     where.split  = data.cont.cont$where.split,
                                                                     dir.split    = data.cont.cont$dir.split)
            eval.final.estg.glm.modout.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modout.true.cv1                <- unlist(eval.final.estg.glm.modout.true.cv1)
            names(eval.final.estg.glm.modout.true.cv1)         <- paste0("homo.g.glm.modout.true.cv1.", names(eval.final.estg.glm.modout.true.cv1))
            print("7.cv1.T")
            
            #####################################################################################################################
            ############################# 8. g: GLM Model, outside node, Noisy adjustment model, cv1 ############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modout.nois.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                    est.used          = "G",
                                                                    type.var          = "cont",
                                                                    adj.mod.loc       = "out",
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = NULL,
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)
            
            final.tree.estg.glm.modout.nois.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont,
                                                                  tree.list         = seq.created.estg.glm.modout.nois.cv1$tree.list,
                                                                  lambda.used       = qchisq(0.95, 1),
                                                                  val.sample        = data.validation.cont.cont,              
                                                                  type.var          = "cont",
                                                                  adj.mod.loc       = "out",
                                                                  adj.mthd          = "GLM",
                                                                  adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modout.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.nois.cv1[[1]],
                                                                     test.data    = data.cont.cont$test.data,
                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                     noise.var    = data.cont.cont$noise.var,
                                                                     corr.split   = data.cont.cont$corr.split,
                                                                     where.split  = data.cont.cont$where.split,
                                                                     dir.split    = data.cont.cont$dir.split)
            eval.final.estg.glm.modout.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modout.nois.cv1                <- unlist(eval.final.estg.glm.modout.nois.cv1)
            names(eval.final.estg.glm.modout.nois.cv1)         <- paste0("homo.g.glm.modout.nois.cv1.", names(eval.final.estg.glm.modout.nois.cv1))
            print("8.cv1.nois")
            
            #####################################################################################################################
            ########################## 9. g: GLM Model, outside node, Misspecified adjustment model, cv1 ########################
            #####################################################################################################################
            data.used.cont.cont.mis <- data.used.cont.cont %>%
              select(-X2)
            data.validation.cont.cont.mis <- data.validation.cont.cont %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estg.glm.modout.mis.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                                   est.used          = "G",
                                                                   type.var          = "cont",
                                                                   adj.mod.loc       = "out",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = NULL,
                                                                   num.truc.obs      = 30,
                                                                   min.node          = 20)
            
            final.tree.estg.glm.modout.mis.cv1 <- EstG.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                                 tree.list         = seq.created.estg.glm.modout.mis.cv1$tree.list,
                                                                 lambda.used       = qchisq(0.95, 1),
                                                                 val.sample        = data.validation.cont.cont.mis,              
                                                                 type.var          = "cont",
                                                                 adj.mod.loc       = "out",
                                                                 adj.mthd          = "GLM",
                                                                 adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modout.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.mis.cv1[[1]],
                                                                    test.data    = data.cont.cont$test.data,
                                                                    true.trt.eff = data.cont.cont$true.trt.eff,
                                                                    noise.var    = data.cont.cont$noise.var,
                                                                    corr.split   = data.cont.cont$corr.split,
                                                                    where.split  = data.cont.cont$where.split,
                                                                    dir.split    = data.cont.cont$dir.split)
            eval.final.estg.glm.modout.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modout.mis.cv1                <- unlist(eval.final.estg.glm.modout.mis.cv1)
            names(eval.final.estg.glm.modout.mis.cv1)         <- paste0("homo.g.glm.modout.mis.cv1.", names(eval.final.estg.glm.modout.mis.cv1))
            print("9.cv1.mis")
            
            #####################################################################################################################
            ############################ 10. g: GLM Model, outside node, True adjustment model, Cv2 #############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modout.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                    est.used          = "G",
                                                                    type.var          = "cont",
                                                                    adj.mod.loc       = "out",
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)
            
            final.tree.estg.glm.modout.true.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                  tree.list         = seq.created.estg.glm.modout.true.cv2$tree.list,             
                                                                  type.var          = "cont",
                                                                  seed              = a[i],
                                                                  n.cv              = 5,
                                                                  adj.mod.loc       = "out",
                                                                  adj.mthd          = "GLM",
                                                                  adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)")
            t1 <- Sys.time()
            
            eval.final.estg.glm.modout.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.true.cv2[[1]],
                                                                     test.data    = data.cont.cont$test.data,
                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                     noise.var    = data.cont.cont$noise.var,
                                                                     corr.split   = data.cont.cont$corr.split,
                                                                     where.split  = data.cont.cont$where.split,
                                                                     dir.split    = data.cont.cont$dir.split)
            
            eval.final.estg.glm.modout.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modout.true.cv2                <- unlist(eval.final.estg.glm.modout.true.cv2)
            names(eval.final.estg.glm.modout.true.cv2)         <- paste0("homo.g.glm.modout.true.cv2.", names(eval.final.estg.glm.modout.true.cv2))
            print("10.cv2.T")
            
            #####################################################################################################################
            ############################# 11. g: GLM Model, outside node, Noisy adjustment model, Cv2 ###########################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modout.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                    est.used          = "G",
                                                                    type.var          = "cont",
                                                                    adj.mod.loc       = "out",
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = NULL,
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)
            
            final.tree.estg.glm.modout.nois.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                  tree.list         = seq.created.estg.glm.modout.nois.cv2$tree.list,             
                                                                  type.var          = "cont",
                                                                  seed              = a[i],
                                                                  n.cv              = 5,
                                                                  adj.mod.loc       = "out",
                                                                  adj.mthd          = "GLM",
                                                                  adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modout.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.nois.cv2[[1]],
                                                                     test.data    = data.cont.cont$test.data,
                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                     noise.var    = data.cont.cont$noise.var,
                                                                     corr.split   = data.cont.cont$corr.split,
                                                                     where.split  = data.cont.cont$where.split,
                                                                     dir.split    = data.cont.cont$dir.split)
            
            eval.final.estg.glm.modout.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modout.nois.cv2                <- unlist(eval.final.estg.glm.modout.nois.cv2)
            names(eval.final.estg.glm.modout.nois.cv2)         <- paste0("homo.g.glm.modout.nois.cv2.", names(eval.final.estg.glm.modout.nois.cv2))
            print("11.cv2.nois")
            
            #####################################################################################################################
            ########################## 12. g: GLM Model, outside node, Misspecified adjustment model, Cv2 #######################
            #####################################################################################################################
            data.used.full.cont.cont.mis <- data.used.full.cont.cont %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estg.glm.modout.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                                   est.used          = "G",
                                                                   type.var          = "cont",
                                                                   adj.mod.loc       = "out",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = NULL,
                                                                   num.truc.obs      = 30,
                                                                   min.node          = 20)
            
            final.tree.estg.glm.modout.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                                 tree.list         = seq.created.estg.glm.modout.mis.cv2$tree.list,             
                                                                 type.var          = "cont",
                                                                 seed              = a[i],
                                                                 n.cv              = 5,
                                                                 adj.mod.loc       = "out",
                                                                 adj.mthd          = "GLM",
                                                                 adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modout.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modout.mis.cv2[[1]],
                                                                    test.data    = data.cont.cont$test.data,
                                                                    true.trt.eff = data.cont.cont$true.trt.eff,
                                                                    noise.var    = data.cont.cont$noise.var,
                                                                    corr.split   = data.cont.cont$corr.split,
                                                                    where.split  = data.cont.cont$where.split,
                                                                    dir.split    = data.cont.cont$dir.split)
            
            eval.final.estg.glm.modout.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modout.mis.cv2                <- unlist(eval.final.estg.glm.modout.mis.cv2)
            names(eval.final.estg.glm.modout.mis.cv2)         <- paste0("homo.g.glm.modout.mis.cv2.", names(eval.final.estg.glm.modout.mis.cv2))
            print("12.cv2.mis")
            
            performance.homo.g <- c(eval.final.estg.glm.modout.true.cv1,
                                    eval.final.estg.glm.modout.nois.cv1,
                                    eval.final.estg.glm.modout.mis.cv1,
                                    eval.final.estg.glm.modout.true.cv2,
                                    eval.final.estg.glm.modout.nois.cv2,
                                    eval.final.estg.glm.modout.mis.cv2)
            
            # print must be put before output, otherwise the output will be print output
            if (i%%30 == 0) {print(i)}
            
            c(performance.hetero.g, performance.homo.g)
            
          }

save(performance.gOut, file = paste0("../Data/AppendixC5/OutG.RData"))
