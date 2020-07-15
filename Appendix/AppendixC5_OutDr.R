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

performance.drOut <- 
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
            ###################### 1. dr: True GLM Model outside, True propensity score model outside, Cv1 ######################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmOut.propscTGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                              est.used          = "DR",
                                                                              type.var          = "cont",
                                                                              propsc.mod.loc    = "out",
                                                                              propsc.mthd       = "GLM",
                                                                              propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                              adj.mod.loc       = "out", 
                                                                              adj.mthd          = "GLM", 
                                                                              adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                              num.truc.obs      = 30,
                                                                              min.node          = 20)
            
            final.tree.estdr.adjTGlmOut.propscTGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                             tree.list         = seq.created.estdr.adjTGlmOut.propscTGlmOut.cv1$tree.list, 
                                                                             lambda.used       = qchisq(0.95, 1),
                                                                             val.sample        = data.validation.cont.cont,
                                                                             type.var          = "cont",
                                                                             propsc.mod.loc    = "out",
                                                                             propsc.mthd       = "GLM",
                                                                             propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                             adj.mod.loc       = "out",
                                                                             adj.mthd          = "GLM",
                                                                             adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)")
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
            eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1                <- unlist(eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1)
            names(eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1)         <- paste0("hetero.dr.adjTGlmOut.propscTGlmOut.cv1.", names(eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1))
            print("hetero.cv1.T.T")
            
            #####################################################################################################################
            ######################## 2. dr: True GLM Model outside, Noisy propensity score model outside, Cv1 ###################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                                 est.used          = "DR",
                                                                                 type.var          = "cont",
                                                                                 propsc.mod.loc    = "out",
                                                                                 propsc.mthd       = "GLM",
                                                                                 propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                 adj.mod.loc       = "out", 
                                                                                 adj.mthd          = "GLM", 
                                                                                 adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                                 num.truc.obs      = 30,
                                                                                 min.node          = 20)
            
            final.tree.estdr.adjTGlmOut.propscNoisGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                                tree.list         = seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv1$tree.list, 
                                                                                lambda.used       = qchisq(0.95, 1),
                                                                                val.sample        = data.validation.cont.cont,
                                                                                type.var          = "cont",
                                                                                propsc.mod.loc    = "out",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                adj.mod.loc       = "out",
                                                                                adj.mthd          = "GLM",
                                                                                adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)")
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
            eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1                <- unlist(eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1)
            names(eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1)         <- paste0("hetero.dr.adjTGlmOut.propscNoisGlmOut.cv1.", names(eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1))
            print("hetero.cv1.T.Nois")
            
            #####################################################################################################################
            ##################### 3. dr: Noisy GLM Model outside, True propensity score model outside, Cv1 ######################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                                 est.used          = "DR",
                                                                                 type.var          = "cont",
                                                                                 propsc.mod.loc    = "out",
                                                                                 propsc.mthd       = "GLM",
                                                                                 propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                 adj.mod.loc       = "out", 
                                                                                 adj.mthd          = "GLM", 
                                                                                 adj.form.true     = NULL, 
                                                                                 num.truc.obs      = 30,
                                                                                 min.node          = 20)
            
            final.tree.estdr.adjNoisGlmOut.propscTGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                                tree.list         = seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv1$tree.list, 
                                                                                lambda.used       = qchisq(0.95, 1),
                                                                                val.sample        = data.validation.cont.cont,
                                                                                type.var          = "cont",
                                                                                propsc.mod.loc    = "out",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                adj.mod.loc       = "out",
                                                                                adj.mthd          = "GLM",
                                                                                adj.form.true     = NULL)
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
            eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1                <- unlist(eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1)
            names(eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1)         <- paste0("hetero.dr.adjNoisGlmOut.propscTGlmOut.cv1.", names(eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1))
            print("hetero.cv1.Nois.T")
            
            #####################################################################################################################
            ##################### 4. dr: Noisy GLM Model outside, Noisy propensity score model outside, Cv1 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                                    est.used          = "DR",
                                                                                    type.var          = "cont",
                                                                                    propsc.mod.loc    = "out",
                                                                                    propsc.mthd       = "GLM",
                                                                                    propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                    adj.mod.loc       = "out", 
                                                                                    adj.mthd          = "GLM", 
                                                                                    adj.form.true     = NULL, 
                                                                                    num.truc.obs      = 30,
                                                                                    min.node          = 20)
            
            final.tree.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                                   tree.list         = seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1$tree.list, 
                                                                                   lambda.used       = qchisq(0.95, 1),
                                                                                   val.sample        = data.validation.cont.cont,
                                                                                   type.var          = "cont",
                                                                                   propsc.mod.loc    = "out",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                   adj.mod.loc       = "out",
                                                                                   adj.mthd          = "GLM",
                                                                                   adj.form.true     = NULL)
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
            eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1                <- unlist(eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1)
            names(eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1)         <- paste0("hetero.dr.adjNoisGlmOut.propscNoisGlmOut.cv1.", names(eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1))
            print("hetero.cv1.Nois.Nois")
            
            #####################################################################################################################
            ############## 5. dr: Misspecified GLM Model outside, Misspecified propensity score model outside, Cv1 ##############
            #####################################################################################################################
            data.used.cont.cont.mis <- data.used.cont.cont %>%
              select(-X2)
            data.validation.cont.cont.mis <- data.validation.cont.cont %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estdr.adjFGlmOut.propscFGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                                              est.used          = "DR",
                                                                              type.var          = "cont",
                                                                              propsc.mod.loc    = "out",
                                                                              propsc.mthd       = "GLM",
                                                                              propsc.form.true  = NULL,
                                                                              adj.mod.loc       = "out", 
                                                                              adj.mthd          = "GLM", 
                                                                              adj.form.true     = NULL, 
                                                                              num.truc.obs      = 30,
                                                                              min.node          = 20)
            
            final.tree.estdr.adjFGlmOut.propscFGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                                             tree.list         = seq.created.estdr.adjFGlmOut.propscFGlmOut.cv1$tree.list, 
                                                                             lambda.used       = qchisq(0.95, 1),
                                                                             val.sample        = data.validation.cont.cont.mis,
                                                                             type.var          = "cont",
                                                                             propsc.mod.loc    = "out",
                                                                             propsc.mthd       = "GLM",
                                                                             propsc.form.true  = NULL,
                                                                             adj.mod.loc       = "out",
                                                                             adj.mthd          = "GLM",
                                                                             adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmOut.propscFGlmOut.cv1[[1]],
                                                                               test.data    = data.cont.cont$test.data,
                                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                                               noise.var    = data.cont.cont$noise.var,
                                                                               corr.split   = data.cont.cont$corr.split,
                                                                               where.split  = data.cont.cont$where.split,
                                                                               dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjFGlmOut.propscFGlmOut.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
            eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1                <- unlist(eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1)
            names(eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1)         <- paste0("hetero.dr.adjFGlmOut.propscFGlmOut.cv1.", names(eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1))
            print("hetero.cv1.F.F")
            
            #####################################################################################################################
            ####################### 6. dr: True GLM Model outside, True propensity score model outside, Cv2 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmOut.propscTGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                              est.used          = "DR",
                                                                              type.var          = "cont",
                                                                              propsc.mod.loc    = "out",
                                                                              propsc.mthd       = "GLM",
                                                                              propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                              adj.mod.loc       = "out", 
                                                                              adj.mthd          = "GLM", 
                                                                              adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                              num.truc.obs      = 30,
                                                                              min.node          = 20)
            
            final.tree.estdr.adjTGlmOut.propscTGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                             tree.list         = seq.created.estdr.adjTGlmOut.propscTGlmOut.cv2$tree.list,
                                                                             type.var          = "cont",
                                                                             seed              = a[i], 
                                                                             n.cv              = 5,
                                                                             propsc.mod.loc    = "out",
                                                                             propsc.mthd       = "GLM",
                                                                             propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                             adj.mod.loc       = "out",
                                                                             adj.mthd          = "GLM",
                                                                             adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)")
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
            eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2                <- unlist(eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2)
            names(eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2)         <- paste0("hetero.dr.adjTGlmOut.propscTGlmOut.cv2.", names(eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2))
            print("hetero.cv2.T.T")
            
            #####################################################################################################################
            ####################### 7. dr: True GLM Model outside, Noisy propensity score model outside, Cv2 ####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                                 est.used          = "DR",
                                                                                 type.var          = "cont",
                                                                                 propsc.mod.loc    = "out",
                                                                                 propsc.mthd       = "GLM",
                                                                                 propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                 adj.mod.loc       = "out", 
                                                                                 adj.mthd          = "GLM", 
                                                                                 adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                                 num.truc.obs      = 30,
                                                                                 min.node          = 20)
            
            final.tree.estdr.adjTGlmOut.propscNoisGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                                tree.list         = seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv2$tree.list,
                                                                                type.var          = "cont",
                                                                                seed              = a[i], 
                                                                                n.cv              = 5,
                                                                                propsc.mod.loc    = "out",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                adj.mod.loc       = "out",
                                                                                adj.mthd          = "GLM",
                                                                                adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)")
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
            eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2                <- unlist(eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2)
            names(eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2)         <- paste0("hetero.dr.adjTGlmOut.propscNoisGlmOut.cv2.", names(eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2))
            print("hetero.cv2.T.Nois")
            
            #####################################################################################################################
            ###################### 8. dr: Noisy GLM Model outside, True propensity score model outside, Cv2 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                                 est.used          = "DR",
                                                                                 type.var          = "cont",
                                                                                 propsc.mod.loc    = "out",
                                                                                 propsc.mthd       = "GLM",
                                                                                 propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                 adj.mod.loc       = "out", 
                                                                                 adj.mthd          = "GLM", 
                                                                                 adj.form.true     = NULL, 
                                                                                 num.truc.obs      = 30,
                                                                                 min.node          = 20)
            
            final.tree.estdr.adjNoisGlmOut.propscTGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                                tree.list         = seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv2$tree.list,
                                                                                type.var          = "cont",
                                                                                seed              = a[i], 
                                                                                n.cv              = 5,
                                                                                propsc.mod.loc    = "out",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                adj.mod.loc       = "out",
                                                                                adj.mthd          = "GLM",
                                                                                adj.form.true     = NULL)
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
            eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2                <- unlist(eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2)
            names(eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2)         <- paste0("hetero.dr.adjNoisGlmOut.propscTGlmOut.cv2.", names(eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2))
            print("hetero.cv2.Nois.T")
            
            #####################################################################################################################
            ##################### 9. dr: Noisy GLM Model outside, Noisy propensity score model outside, Cv2 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                                    est.used          = "DR",
                                                                                    type.var          = "cont",
                                                                                    propsc.mod.loc    = "out",
                                                                                    propsc.mthd       = "GLM",
                                                                                    propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                    adj.mod.loc       = "out", 
                                                                                    adj.mthd          = "GLM", 
                                                                                    adj.form.true     = NULL, 
                                                                                    num.truc.obs      = 30,
                                                                                    min.node          = 20)
            
            final.tree.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                                   tree.list         = seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2$tree.list,
                                                                                   type.var          = "cont",
                                                                                   seed              = a[i], 
                                                                                   n.cv              = 5,
                                                                                   propsc.mod.loc    = "out",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                   adj.mod.loc       = "out",
                                                                                   adj.mthd          = "GLM",
                                                                                   adj.form.true     = NULL)
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
            eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2                <- unlist(eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2)
            names(eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2)         <- paste0("hetero.dr.adjNoisGlmOut.propscNoisGlmOut.cv2.", names(eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2))
            print("hetero.cv2.Nois.Nois")
            
            #####################################################################################################################
            ############# 10. dr: Misspecified GLM Model outside, Misspecified propensity score model outside, Cv2 ##############
            #####################################################################################################################
            data.used.full.cont.cont.mis <- data.used.full.cont.cont %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estdr.adjFGlmOut.propscFGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                                              est.used          = "DR",
                                                                              type.var          = "cont",
                                                                              propsc.mod.loc    = "out",
                                                                              propsc.mthd       = "GLM",
                                                                              propsc.form.true  = NULL,
                                                                              adj.mod.loc       = "out", 
                                                                              adj.mthd          = "GLM", 
                                                                              adj.form.true     = NULL, 
                                                                              num.truc.obs      = 30,
                                                                              min.node          = 20)
            
            final.tree.estdr.adjFGlmOut.propscFGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                                             tree.list         = seq.created.estdr.adjFGlmOut.propscFGlmOut.cv2$tree.list,
                                                                             type.var          = "cont",
                                                                             seed              = a[i], 
                                                                             n.cv              = 5,
                                                                             propsc.mod.loc    = "out",
                                                                             propsc.mthd       = "GLM",
                                                                             propsc.form.true  = NULL,
                                                                             adj.mod.loc       = "out",
                                                                             adj.mthd          = "GLM",
                                                                             adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmOut.propscFGlmOut.cv2[[1]],
                                                                               test.data    = data.cont.cont$test.data,
                                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                                               noise.var    = data.cont.cont$noise.var,
                                                                               corr.split   = data.cont.cont$corr.split,
                                                                               where.split  = data.cont.cont$where.split,
                                                                               dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjFGlmOut.propscFGlmOut.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
            eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2                <- unlist(eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2)
            names(eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2)         <- paste0("hetero.dr.adjFGlmOut.propscFGlmOut.cv2.", names(eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2))
            print("hetero.cv2.F.F")
            
            performance.hetero.drOut <- c(eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1,
                                          eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1,
                                          eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1,
                                          eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1,
                                          eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1,
                                          eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2,
                                          eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2,
                                          eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2,
                                          eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2,
                                          eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2)
            
            
            
            #####################################################################################################################
            ################################################### Homogeneous #####################################################
            #####################################################################################################################
            data.cont.cont            <- makeData.cont.noeff.cont(1000, 1000, coeff.prop.sc = 0.6)
            data.used.full.cont.cont  <- data.cont.cont$data.used
            data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
            # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
            data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ] 
            
            #####################################################################################################################
            ###################### 11. dr: True GLM Model outside, True propensity score model outside, Cv1 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmOut.propscTGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                              est.used          = "DR",
                                                                              type.var          = "cont",
                                                                              propsc.mod.loc    = "out",
                                                                              propsc.mthd       = "GLM",
                                                                              propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                              adj.mod.loc       = "out", 
                                                                              adj.mthd          = "GLM", 
                                                                              adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                              num.truc.obs      = 30,
                                                                              min.node          = 20)
            
            final.tree.estdr.adjTGlmOut.propscTGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                             tree.list         = seq.created.estdr.adjTGlmOut.propscTGlmOut.cv1$tree.list, 
                                                                             lambda.used       = qchisq(0.95, 1),
                                                                             val.sample        = data.validation.cont.cont,
                                                                             type.var          = "cont",
                                                                             propsc.mod.loc    = "out",
                                                                             propsc.mthd       = "GLM",
                                                                             propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                             adj.mod.loc       = "out",
                                                                             adj.mthd          = "GLM",
                                                                             adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmOut.propscTGlmOut.cv1[[1]],
                                                                               test.data    = data.cont.cont$test.data,
                                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                                               noise.var    = data.cont.cont$noise.var,
                                                                               corr.split   = data.cont.cont$corr.split,
                                                                               where.split  = data.cont.cont$where.split,
                                                                               dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1                <- unlist(eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1)
            names(eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1)         <- paste0("homo.dr.adjTGlmOut.propscTGlmOut.cv1.", names(eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1))
            print("homo.cv1.T.T")
            
            #####################################################################################################################
            ######################## 12. dr: True GLM Model outside, Noisy propensity score model outside, Cv1 ##################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                                 est.used          = "DR",
                                                                                 type.var          = "cont",
                                                                                 propsc.mod.loc    = "out",
                                                                                 propsc.mthd       = "GLM",
                                                                                 propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                 adj.mod.loc       = "out", 
                                                                                 adj.mthd          = "GLM", 
                                                                                 adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                                 num.truc.obs      = 30,
                                                                                 min.node          = 20)
            
            final.tree.estdr.adjTGlmOut.propscNoisGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                                tree.list         = seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv1$tree.list, 
                                                                                lambda.used       = qchisq(0.95, 1),
                                                                                val.sample        = data.validation.cont.cont,
                                                                                type.var          = "cont",
                                                                                propsc.mod.loc    = "out",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                adj.mod.loc       = "out",
                                                                                adj.mthd          = "GLM",
                                                                                adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmOut.propscNoisGlmOut.cv1[[1]],
                                                                                  test.data    = data.cont.cont$test.data,
                                                                                  true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                  noise.var    = data.cont.cont$noise.var,
                                                                                  corr.split   = data.cont.cont$corr.split,
                                                                                  where.split  = data.cont.cont$where.split,
                                                                                  dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1                <- unlist(eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1)
            names(eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1)         <- paste0("homo.dr.adjTGlmOut.propscNoisGlmOut.cv1.", names(eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1))
            print("homo.cv1.T.Nois")
            
            #####################################################################################################################
            ###################### 13. dr: Noisy GLM Model outside, True propensity score model outside, Cv1 ####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                                 est.used          = "DR",
                                                                                 type.var          = "cont",
                                                                                 propsc.mod.loc    = "out",
                                                                                 propsc.mthd       = "GLM",
                                                                                 propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                 adj.mod.loc       = "out", 
                                                                                 adj.mthd          = "GLM", 
                                                                                 adj.form.true     = NULL, 
                                                                                 num.truc.obs      = 30,
                                                                                 min.node          = 20)
            
            final.tree.estdr.adjNoisGlmOut.propscTGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                                tree.list         = seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv1$tree.list, 
                                                                                lambda.used       = qchisq(0.95, 1),
                                                                                val.sample        = data.validation.cont.cont,
                                                                                type.var          = "cont",
                                                                                propsc.mod.loc    = "out",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                adj.mod.loc       = "out",
                                                                                adj.mthd          = "GLM",
                                                                                adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmOut.propscTGlmOut.cv1[[1]],
                                                                                  test.data    = data.cont.cont$test.data,
                                                                                  true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                  noise.var    = data.cont.cont$noise.var,
                                                                                  corr.split   = data.cont.cont$corr.split,
                                                                                  where.split  = data.cont.cont$where.split,
                                                                                  dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1                <- unlist(eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1)
            names(eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1)         <- paste0("homo.dr.adjNoisGlmOut.propscTGlmOut.cv1.", names(eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1))
            print("homo.cv1.Nois.T")
            
            #####################################################################################################################
            ####################### 14. dr: Noisy GLM Model outside, Noisy propensity score model outside, Cv1 ##################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                                    est.used          = "DR",
                                                                                    type.var          = "cont",
                                                                                    propsc.mod.loc    = "out",
                                                                                    propsc.mthd       = "GLM",
                                                                                    propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                    adj.mod.loc       = "out", 
                                                                                    adj.mthd          = "GLM", 
                                                                                    adj.form.true     = NULL, 
                                                                                    num.truc.obs      = 30,
                                                                                    min.node          = 20)
            
            final.tree.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                                   tree.list         = seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1$tree.list, 
                                                                                   lambda.used       = qchisq(0.95, 1),
                                                                                   val.sample        = data.validation.cont.cont,
                                                                                   type.var          = "cont",
                                                                                   propsc.mod.loc    = "out",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                   adj.mod.loc       = "out",
                                                                                   adj.mthd          = "GLM",
                                                                                   adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1[[1]],
                                                                                     test.data    = data.cont.cont$test.data,
                                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                     noise.var    = data.cont.cont$noise.var,
                                                                                     corr.split   = data.cont.cont$corr.split,
                                                                                     where.split  = data.cont.cont$where.split,
                                                                                     dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1                <- unlist(eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1)
            names(eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1)         <- paste0("homo.dr.adjNoisGlmOut.propscNoisGlmOut.cv1.", names(eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1))
            print("homo.cv1.Nois.Nois")
            
            #####################################################################################################################
            ############## 15. dr: Misspecified GLM Model outside, Misspecified propensity score model outside, Cv1 #############
            #####################################################################################################################
            data.used.cont.cont.mis <- data.used.cont.cont %>%
              select(-X2)
            data.validation.cont.cont.mis <- data.validation.cont.cont %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estdr.adjFGlmOut.propscFGlmOut.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                                              est.used          = "DR",
                                                                              type.var          = "cont",
                                                                              propsc.mod.loc    = "out",
                                                                              propsc.mthd       = "GLM",
                                                                              propsc.form.true  = NULL,
                                                                              adj.mod.loc       = "out", 
                                                                              adj.mthd          = "GLM", 
                                                                              adj.form.true     = NULL, 
                                                                              num.truc.obs      = 30,
                                                                              min.node          = 20)
            
            final.tree.estdr.adjFGlmOut.propscFGlmOut.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                                             tree.list         = seq.created.estdr.adjFGlmOut.propscFGlmOut.cv1$tree.list, 
                                                                             lambda.used       = qchisq(0.95, 1),
                                                                             val.sample        = data.validation.cont.cont.mis,
                                                                             type.var          = "cont",
                                                                             propsc.mod.loc    = "out",
                                                                             propsc.mthd       = "GLM",
                                                                             propsc.form.true  = NULL,
                                                                             adj.mod.loc       = "out",
                                                                             adj.mthd          = "GLM",
                                                                             adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmOut.propscFGlmOut.cv1[[1]],
                                                                               test.data    = data.cont.cont$test.data,
                                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                                               noise.var    = data.cont.cont$noise.var,
                                                                               corr.split   = data.cont.cont$corr.split,
                                                                               where.split  = data.cont.cont$where.split,
                                                                               dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1                <- unlist(eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1)
            names(eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1)         <- paste0("homo.dr.adjFGlmOut.propscFGlmOut.cv1.", names(eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1))
            print("homo.cv1.F.F")
            
            #####################################################################################################################
            ######################## 16. dr: True GLM Model outside, True propensity score model outside, Cv2 ###################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmOut.propscTGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                              est.used          = "DR",
                                                                              type.var          = "cont",
                                                                              propsc.mod.loc    = "out",
                                                                              propsc.mthd       = "GLM",
                                                                              propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                              adj.mod.loc       = "out", 
                                                                              adj.mthd          = "GLM", 
                                                                              adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                              num.truc.obs      = 30,
                                                                              min.node          = 20)
            
            final.tree.estdr.adjTGlmOut.propscTGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                             tree.list         = seq.created.estdr.adjTGlmOut.propscTGlmOut.cv2$tree.list,
                                                                             type.var          = "cont",
                                                                             seed              = a[i], 
                                                                             n.cv              = 5,
                                                                             propsc.mod.loc    = "out",
                                                                             propsc.mthd       = "GLM",
                                                                             propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                             adj.mod.loc       = "out",
                                                                             adj.mthd          = "GLM",
                                                                             adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmOut.propscTGlmOut.cv2[[1]],
                                                                               test.data    = data.cont.cont$test.data,
                                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                                               noise.var    = data.cont.cont$noise.var,
                                                                               corr.split   = data.cont.cont$corr.split,
                                                                               where.split  = data.cont.cont$where.split,
                                                                               dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2                <- unlist(eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2)
            names(eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2)         <- paste0("homo.dr.adjTGlmOut.propscTGlmOut.cv2.", names(eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2))
            print("homo.cv2.T.T")
            
            #####################################################################################################################
            ####################### 17. dr: True GLM Model outside, Noisy propensity score model outside, Cv2 ###################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                                 est.used          = "DR",
                                                                                 type.var          = "cont",
                                                                                 propsc.mod.loc    = "out",
                                                                                 propsc.mthd       = "GLM",
                                                                                 propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                 adj.mod.loc       = "out", 
                                                                                 adj.mthd          = "GLM", 
                                                                                 adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                                 num.truc.obs      = 30,
                                                                                 min.node          = 20)
            
            final.tree.estdr.adjTGlmOut.propscNoisGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                                tree.list         = seq.created.estdr.adjTGlmOut.propscNoisGlmOut.cv2$tree.list,
                                                                                type.var          = "cont",
                                                                                seed              = a[i], 
                                                                                n.cv              = 5,
                                                                                propsc.mod.loc    = "out",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                adj.mod.loc       = "out",
                                                                                adj.mthd          = "GLM",
                                                                                adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmOut.propscNoisGlmOut.cv2[[1]],
                                                                                  test.data    = data.cont.cont$test.data,
                                                                                  true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                  noise.var    = data.cont.cont$noise.var,
                                                                                  corr.split   = data.cont.cont$corr.split,
                                                                                  where.split  = data.cont.cont$where.split,
                                                                                  dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2                <- unlist(eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2)
            names(eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2)         <- paste0("homo.dr.adjTGlmOut.propscNoisGlmOut.cv2.", names(eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2))
            print("homo.cv2.T.Nois")
            
            #####################################################################################################################
            ####################### 18. dr: Noisy GLM Model outside, True propensity score model outside, Cv2 ###################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                                 est.used          = "DR",
                                                                                 type.var          = "cont",
                                                                                 propsc.mod.loc    = "out",
                                                                                 propsc.mthd       = "GLM",
                                                                                 propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                 adj.mod.loc       = "out", 
                                                                                 adj.mthd          = "GLM", 
                                                                                 adj.form.true     = NULL, 
                                                                                 num.truc.obs      = 30,
                                                                                 min.node          = 20)
            
            final.tree.estdr.adjNoisGlmOut.propscTGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                                tree.list         = seq.created.estdr.adjNoisGlmOut.propscTGlmOut.cv2$tree.list,
                                                                                type.var          = "cont",
                                                                                seed              = a[i], 
                                                                                n.cv              = 5,
                                                                                propsc.mod.loc    = "out",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                adj.mod.loc       = "out",
                                                                                adj.mthd          = "GLM",
                                                                                adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmOut.propscTGlmOut.cv2[[1]],
                                                                                  test.data    = data.cont.cont$test.data,
                                                                                  true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                  noise.var    = data.cont.cont$noise.var,
                                                                                  corr.split   = data.cont.cont$corr.split,
                                                                                  where.split  = data.cont.cont$where.split,
                                                                                  dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2                <- unlist(eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2)
            names(eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2)         <- paste0("homo.dr.adjNoisGlmOut.propscTGlmOut.cv2.", names(eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2))
            print("homo.cv2.Nois.T")
            
            #####################################################################################################################
            ###################### 19. dr: Noisy GLM Model outside, Noisy propensity score model outside, Cv2 ###################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                                    est.used          = "DR",
                                                                                    type.var          = "cont",
                                                                                    propsc.mod.loc    = "out",
                                                                                    propsc.mthd       = "GLM",
                                                                                    propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                    adj.mod.loc       = "out", 
                                                                                    adj.mthd          = "GLM", 
                                                                                    adj.form.true     = NULL, 
                                                                                    num.truc.obs      = 30,
                                                                                    min.node          = 20)
            
            final.tree.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                                   tree.list         = seq.created.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2$tree.list,
                                                                                   type.var          = "cont",
                                                                                   seed              = a[i], 
                                                                                   n.cv              = 5,
                                                                                   propsc.mod.loc    = "out",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                   adj.mod.loc       = "out",
                                                                                   adj.mthd          = "GLM",
                                                                                   adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2[[1]],
                                                                                     test.data    = data.cont.cont$test.data,
                                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                     noise.var    = data.cont.cont$noise.var,
                                                                                     corr.split   = data.cont.cont$corr.split,
                                                                                     where.split  = data.cont.cont$where.split,
                                                                                     dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2                <- unlist(eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2)
            names(eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2)         <- paste0("homo.dr.adjNoisGlmOut.propscNoisGlmOut.cv2.", names(eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2))
            print("homo.cv2.Nois.Nois")
            
            #####################################################################################################################
            ############## 20. dr: Misspecified GLM Model outside, Misspecified propensity score model outside, Cv2 #############
            #####################################################################################################################
            data.used.full.cont.cont.mis <- data.used.full.cont.cont %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estdr.adjFGlmOut.propscFGlmOut.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                                              est.used          = "DR",
                                                                              type.var          = "cont",
                                                                              propsc.mod.loc    = "out",
                                                                              propsc.mthd       = "GLM",
                                                                              propsc.form.true  = NULL,
                                                                              adj.mod.loc       = "out", 
                                                                              adj.mthd          = "GLM", 
                                                                              adj.form.true     = NULL, 
                                                                              num.truc.obs      = 30,
                                                                              min.node          = 20)
            
            final.tree.estdr.adjFGlmOut.propscFGlmOut.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                                             tree.list         = seq.created.estdr.adjFGlmOut.propscFGlmOut.cv2$tree.list,
                                                                             type.var          = "cont",
                                                                             seed              = a[i], 
                                                                             n.cv              = 5,
                                                                             propsc.mod.loc    = "out",
                                                                             propsc.mthd       = "GLM",
                                                                             propsc.form.true  = NULL,
                                                                             adj.mod.loc       = "out",
                                                                             adj.mthd          = "GLM",
                                                                             adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmOut.propscFGlmOut.cv2[[1]],
                                                                               test.data    = data.cont.cont$test.data,
                                                                               true.trt.eff = data.cont.cont$true.trt.eff,
                                                                               noise.var    = data.cont.cont$noise.var,
                                                                               corr.split   = data.cont.cont$corr.split,
                                                                               where.split  = data.cont.cont$where.split,
                                                                               dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2                <- unlist(eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2)
            names(eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2)         <- paste0("homo.dr.adjFGlmOut.propscFGlmOut.cv2.", names(eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2))
            print("homo.cv2.F.F")
            
            performance.homo.drOut <- c(eval.final.estdr.adjTGlmOut.propscTGlmOut.cv1,
                                        eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv1,
                                        eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv1,
                                        eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv1,
                                        eval.final.estdr.adjFGlmOut.propscFGlmOut.cv1,
                                        eval.final.estdr.adjTGlmOut.propscTGlmOut.cv2,
                                        eval.final.estdr.adjTGlmOut.propscNoisGlmOut.cv2,
                                        eval.final.estdr.adjNoisGlmOut.propscTGlmOut.cv2,
                                        eval.final.estdr.adjNoisGlmOut.propscNoisGlmOut.cv2,
                                        eval.final.estdr.adjFGlmOut.propscFGlmOut.cv2)
            
            # print must be put before output, otherwise the output will be print output
            if (i%%30 == 0) {print(i)}
            
            c(performance.hetero.drOut, performance.homo.drOut)
            
          }

save(performance.drOut, file = paste0("../Data/AppendixC5/OutDr.RData"))
