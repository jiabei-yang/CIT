#!/usr/bin/env Rscript
setwd("../")
folder <- paste(getwd(), "/Functions/", sep="")
functions <- list.files(folder)
functions <- grep("main.R", functions, invert = T, value = T)
functions <- paste(folder, functions, sep = "")
for (i in functions){
  source(i)
}

setwd("main/")

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

performance.drInnd <- 
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
            ###################### 1. dr: True GLM Model in node, True propensity score model in node, Cv1 ######################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                                est.used          = "DR",
                                                                                type.var          = "cont",
                                                                                propsc.mod.loc    = "node",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                adj.mod.loc       = "node", 
                                                                                adj.mthd          = "GLM", 
                                                                                adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                                num.truc.obs      = 30,
                                                                                min.node          = 20)
            
            final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                               tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                               lambda.used       = qchisq(0.95, 1),
                                                                               val.sample        = data.validation.cont.cont,
                                                                               type.var          = "cont",
                                                                               propsc.mod.loc    = "node",
                                                                               propsc.mthd       = "GLM",
                                                                               propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                               adj.mod.loc       = "node",
                                                                               adj.mthd          = "GLM",
                                                                               adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)")
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
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1                <- unlist(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1)
            names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1)         <- paste0("hetero.dr.adjTGlmInnd.propscTGlmInnd.cv1.", names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1))
            print("1")
            
            #####################################################################################################################
            #################### 2. dr: True GLM Model in node, Mis func propensity score model in node, Cv1 ####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                                   est.used          = "DR",
                                                                                   type.var          = "cont",
                                                                                   propsc.mod.loc    = "node",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                   adj.mod.loc       = "node", 
                                                                                   adj.mthd          = "GLM", 
                                                                                   adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)", 
                                                                                   num.truc.obs      = 30,
                                                                                   min.node          = 20)
            
            final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                                  tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                                  lambda.used       = qchisq(0.95, 1),
                                                                                  val.sample        = data.validation.cont.cont,
                                                                                  type.var          = "cont",
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                  adj.mod.loc       = "node",
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)")
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
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1                <- unlist(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1)
            names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1)         <- paste0("hetero.dr.adjTGlmInnd.propscNoisGlmInnd.cv1.", names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1))
            print("2")
            
            #####################################################################################################################
            ##################### 3. dr: Mis func GLM Model in node, True propensity score model in node, Cv1 ###################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                                   est.used          = "DR",
                                                                                   type.var          = "cont",
                                                                                   propsc.mod.loc    = "node",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                   adj.mod.loc       = "node", 
                                                                                   adj.mthd          = "GLM", 
                                                                                   adj.form.true     = NULL, 
                                                                                   num.truc.obs      = 30,
                                                                                   min.node          = 20)
            
            final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                                  tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                                  lambda.used       = qchisq(0.95, 1),
                                                                                  val.sample        = data.validation.cont.cont,
                                                                                  type.var          = "cont",
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                  adj.mod.loc       = "node",
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = NULL)
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
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1)
            names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1)         <- paste0("hetero.dr.adjNoisGlmInnd.propscTGlmInnd.cv1.", names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1))
            print("3")
            
            #####################################################################################################################
            ################## 4. dr: Mis func GLM Model in node, Mis func propensity score model in node, Cv1 ##################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                                      est.used          = "DR",
                                                                                      type.var          = "cont",
                                                                                      propsc.mod.loc    = "node",
                                                                                      propsc.mthd       = "GLM",
                                                                                      propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                      adj.mod.loc       = "node", 
                                                                                      adj.mthd          = "GLM", 
                                                                                      adj.form.true     = NULL, 
                                                                                      num.truc.obs      = 30,
                                                                                      min.node          = 20)
            
            final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                                     tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                                     lambda.used       = qchisq(0.95, 1),
                                                                                     val.sample        = data.validation.cont.cont,
                                                                                     type.var          = "cont",
                                                                                     propsc.mod.loc    = "node",
                                                                                     propsc.mthd       = "GLM",
                                                                                     propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                     adj.mod.loc       = "node",
                                                                                     adj.mthd          = "GLM",
                                                                                     adj.form.true     = NULL)
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
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1)
            names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1)         <- paste0("hetero.dr.adjNoisGlmInnd.propscNoisGlmInnd.cv1.", names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1))
            print("4")
            
            #####################################################################################################################
            ############ 5. dr: Unmeasured cov GLM Model in node, Unmeasured cov propensity score model in node, Cv1 ############
            #####################################################################################################################
            data.used.cont.cont.mis <- data.used.cont.cont %>%
              select(-X2)
            data.validation.cont.cont.mis <- data.validation.cont.cont %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                                                est.used          = "DR",
                                                                                type.var          = "cont",
                                                                                propsc.mod.loc    = "node",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = NULL,
                                                                                adj.mod.loc       = "node", 
                                                                                adj.mthd          = "GLM", 
                                                                                adj.form.true     = NULL, 
                                                                                num.truc.obs      = 30,
                                                                                min.node          = 20)
            
            final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                                               tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1$tree.list, 
                                                                               lambda.used       = qchisq(0.95, 1),
                                                                               val.sample        = data.validation.cont.cont.mis,
                                                                               type.var          = "cont",
                                                                               propsc.mod.loc    = "node",
                                                                               propsc.mthd       = "GLM",
                                                                               propsc.form.true  = NULL,
                                                                               adj.mod.loc       = "node",
                                                                               adj.mthd          = "GLM",
                                                                               adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1[[1]],
                                                                                 test.data    = data.cont.cont$test.data,
                                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                 noise.var    = data.cont.cont$noise.var,
                                                                                 corr.split   = data.cont.cont$corr.split,
                                                                                 where.split  = data.cont.cont$where.split,
                                                                                 dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1$corr.frst.splt <- as.character(seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1                <- unlist(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1)
            names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1)         <- paste0("hetero.dr.adjFGlmInnd.propscFGlmInnd.cv1.", names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1))
            print("5")
            
            performance.hetero.drInnd <- c(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1, 
                                           eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1,
                                           eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1,
                                           eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1,
                                           eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1)
            
            #####################################################################################################################
            ################################################### Homogeneous #####################################################
            #####################################################################################################################
            data.cont.cont            <- makeData.cont.noeff.cont(1000, 1000, coeff.prop.sc = 0.6)
            data.used.full.cont.cont  <- data.cont.cont$data.used
            data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
            # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
            data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ] 
            
            #####################################################################################################################
            ###################### 6. dr: True GLM Model in node, True propensity score model in node, Cv1 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                                est.used          = "DR",
                                                                                type.var          = "cont",
                                                                                propsc.mod.loc    = "node",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                adj.mod.loc       = "node", 
                                                                                adj.mthd          = "GLM", 
                                                                                adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                                num.truc.obs      = 30,
                                                                                min.node          = 20)
            
            final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                               tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                               lambda.used       = qchisq(0.95, 1),
                                                                               val.sample        = data.validation.cont.cont,
                                                                               type.var          = "cont",
                                                                               propsc.mod.loc    = "node",
                                                                               propsc.mthd       = "GLM",
                                                                               propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                               adj.mod.loc       = "node",
                                                                               adj.mthd          = "GLM",
                                                                               adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1[[1]],
                                                                                 test.data    = data.cont.cont$test.data,
                                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                 noise.var    = data.cont.cont$noise.var,
                                                                                 corr.split   = data.cont.cont$corr.split,
                                                                                 where.split  = data.cont.cont$where.split,
                                                                                 dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1                <- unlist(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1)
            names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1)         <- paste0("homo.dr.adjTGlmInnd.propscTGlmInnd.cv1.", names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1))
            print("6")
            
            #####################################################################################################################
            ##################### 7. dr: True GLM Model in node, Mis func propensity score model in node, Cv1 ###################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                                   est.used          = "DR",
                                                                                   type.var          = "cont",
                                                                                   propsc.mod.loc    = "node",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                   adj.mod.loc       = "node", 
                                                                                   adj.mthd          = "GLM", 
                                                                                   adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)", 
                                                                                   num.truc.obs      = 30,
                                                                                   min.node          = 20)
            
            final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                                  tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                                  lambda.used       = qchisq(0.95, 1),
                                                                                  val.sample        = data.validation.cont.cont,
                                                                                  type.var          = "cont",
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                  adj.mod.loc       = "node",
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1[[1]],
                                                                                    test.data    = data.cont.cont$test.data,
                                                                                    true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                    noise.var    = data.cont.cont$noise.var,
                                                                                    corr.split   = data.cont.cont$corr.split,
                                                                                    where.split  = data.cont.cont$where.split,
                                                                                    dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1                <- unlist(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1)
            names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1)         <- paste0("homo.dr.adjTGlmInnd.propscNoisGlmInnd.cv1.", names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1))
            print("7")
            
            #####################################################################################################################
            #################### 8. dr: Mis func GLM Model in node, True propensity score model in node, Cv1 ####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                                   est.used          = "DR",
                                                                                   type.var          = "cont",
                                                                                   propsc.mod.loc    = "node",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                   adj.mod.loc       = "node", 
                                                                                   adj.mthd          = "GLM", 
                                                                                   adj.form.true     = NULL, 
                                                                                   num.truc.obs      = 30,
                                                                                   min.node          = 20)
            
            final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                                  tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                                  lambda.used       = qchisq(0.95, 1),
                                                                                  val.sample        = data.validation.cont.cont,
                                                                                  type.var          = "cont",
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                  adj.mod.loc       = "node",
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1[[1]],
                                                                                    test.data    = data.cont.cont$test.data,
                                                                                    true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                    noise.var    = data.cont.cont$noise.var,
                                                                                    corr.split   = data.cont.cont$corr.split,
                                                                                    where.split  = data.cont.cont$where.split,
                                                                                    dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1)
            names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1)         <- paste0("homo.dr.adjNoisGlmInnd.propscTGlmInnd.cv1.", names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1))
            print("8")
            
            #####################################################################################################################
            ################## 9. dr: Mis func GLM Model in node, Mis func propensity score model in node, Cv1 ##################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont,
                                                                                      est.used          = "DR",
                                                                                      type.var          = "cont",
                                                                                      propsc.mod.loc    = "node",
                                                                                      propsc.mthd       = "GLM",
                                                                                      propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                      adj.mod.loc       = "node", 
                                                                                      adj.mthd          = "GLM", 
                                                                                      adj.form.true     = NULL, 
                                                                                      num.truc.obs      = 30,
                                                                                      min.node          = 20)
            
            final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont,
                                                                                     tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                                     lambda.used       = qchisq(0.95, 1),
                                                                                     val.sample        = data.validation.cont.cont,
                                                                                     type.var          = "cont",
                                                                                     propsc.mod.loc    = "node",
                                                                                     propsc.mthd       = "GLM",
                                                                                     propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                     adj.mod.loc       = "node",
                                                                                     adj.mthd          = "GLM",
                                                                                     adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1[[1]],
                                                                                       test.data    = data.cont.cont$test.data,
                                                                                       true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                       noise.var    = data.cont.cont$noise.var,
                                                                                       corr.split   = data.cont.cont$corr.split,
                                                                                       where.split  = data.cont.cont$where.split,
                                                                                       dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1)
            names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1)         <- paste0("homo.dr.adjNoisGlmInnd.propscNoisGlmInnd.cv1.", names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1))
            print("9")
            
            #####################################################################################################################
            ########### 10. dr: Unmeasured cov GLM Model in node, Unmeasured cov propensity score model in node, Cv1 ############
            #####################################################################################################################
            data.used.cont.cont.mis <- data.used.cont.cont %>%
              select(-X2)
            data.validation.cont.cont.mis <- data.validation.cont.cont %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- create.sequence(data.used         = data.used.cont.cont.mis,
                                                                                est.used          = "DR",
                                                                                type.var          = "cont",
                                                                                propsc.mod.loc    = "node",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = NULL,
                                                                                adj.mod.loc       = "node", 
                                                                                adj.mthd          = "GLM", 
                                                                                adj.form.true     = NULL, 
                                                                                num.truc.obs      = 30,
                                                                                min.node          = 20)
            
            final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.cont.cont.mis,
                                                                               tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1$tree.list, 
                                                                               lambda.used       = qchisq(0.95, 1),
                                                                               val.sample        = data.validation.cont.cont.mis,
                                                                               type.var          = "cont",
                                                                               propsc.mod.loc    = "node",
                                                                               propsc.mthd       = "GLM",
                                                                               propsc.form.true  = NULL,
                                                                               adj.mod.loc       = "node",
                                                                               adj.mthd          = "GLM",
                                                                               adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1[[1]],
                                                                                 test.data    = data.cont.cont$test.data,
                                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                 noise.var    = data.cont.cont$noise.var,
                                                                                 corr.split   = data.cont.cont$corr.split,
                                                                                 where.split  = data.cont.cont$where.split,
                                                                                 dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1                <- unlist(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1)
            names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1)         <- paste0("homo.dr.adjFGlmInnd.propscFGlmInnd.cv1.", names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1))
            print("10")
            
            performance.homo.drInnd <- c(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1, 
                                         eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1,
                                         eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1,
                                         eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1,
                                         eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1)
            
            # print must be put before output, otherwise the output will be print output
            if (i%%30 == 0) {print(i)}
            
            c(performance.hetero.drInnd, performance.homo.drInnd)
            
          }

save(performance.drInnd, file = paste0("../Data/main/MainDr.RData"))

