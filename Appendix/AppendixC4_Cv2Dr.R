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
            ###################### 1. dr: True GLM Model in node, True propensity score model in node, Cv2 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
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
            
            final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                               tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2$tree.list,
                                                                               type.var          = "cont",
                                                                               seed              = a[i], 
                                                                               n.cv              = 5,
                                                                               propsc.mod.loc    = "node",
                                                                               propsc.mthd       = "GLM",
                                                                               propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                               adj.mod.loc       = "node",
                                                                               adj.mthd          = "GLM",
                                                                               adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)")
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
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2                <- unlist(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2)
            names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2)         <- paste0("hetero.dr.adjTGlmInnd.propscTGlmInnd.cv2.", names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2))
            print("hetero.T.T")
            
            #####################################################################################################################
            #################### 2. dr: True GLM Model in node, Mis func propensity score model in node, Cv1 ####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
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
            
            final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                                  tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$tree.list,
                                                                                  type.var          = "cont",
                                                                                  seed              = a[i], 
                                                                                  n.cv              = 5,
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                  adj.mod.loc       = "node",
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + A:I(X4 > 0) + I(X5^3)")
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
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2                <- unlist(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2)
            names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2)         <- paste0("hetero.dr.adjTGlmInnd.propscNoisGlmInnd.cv2.", names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2))
            print("hetero.T.Nois")
            
            #####################################################################################################################
            ###################### 3. dr: Noisy GLM Model in node, True propensity score model in node, Cv2 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
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
            
            final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                                  tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$tree.list,
                                                                                  type.var          = "cont",
                                                                                  seed              = a[i], 
                                                                                  n.cv              = 5,
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                  adj.mod.loc       = "node", 
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = NULL)
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
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2)
            names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2)         <- paste0("hetero.dr.adjNoisGlmInnd.propscTGlmInnd.cv2.", names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2))
            print("hetero.Nois.T")
            
            #####################################################################################################################
            ##################### 4. dr: Noisy GLM Model in node, Noisy propensity score model in node, Cv2 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
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
            
            final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                                     tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$tree.list,
                                                                                     type.var          = "cont",
                                                                                     seed              = a[i], 
                                                                                     n.cv              = 5,
                                                                                     propsc.mod.loc    = "node",
                                                                                     propsc.mthd       = "GLM",
                                                                                     propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                     adj.mod.loc       = "node",
                                                                                     adj.mthd          = "GLM",
                                                                                     adj.form.true     = NULL)
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
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2)
            names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2)         <- paste0("hetero.dr.adjNoisGlmInnd.propscNoisGlmInnd.cv2.", names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2))
            print("hetero.Nois.Nois")
            
            #####################################################################################################################
            ############## 5. dr: Misspecified GLM Model in node, Misspecified propensity score model in node, Cv2 ##############
            #####################################################################################################################
            data.used.full.cont.cont.mis <- data.used.full.cont.cont %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
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
            
            final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                                               tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2$tree.list,
                                                                               type.var          = "cont",
                                                                               seed              = a[i], 
                                                                               n.cv              = 5,
                                                                               propsc.mod.loc    = "node",
                                                                               propsc.mthd       = "GLM",
                                                                               propsc.form.true  = NULL,
                                                                               adj.mod.loc       = "node",
                                                                               adj.mthd          = "GLM",
                                                                               adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv2[[1]],
                                                                                 test.data    = data.cont.cont$test.data,
                                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                 noise.var    = data.cont.cont$noise.var,
                                                                                 corr.split   = data.cont.cont$corr.split,
                                                                                 where.split  = data.cont.cont$where.split,
                                                                                 dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2$corr.frst.splt <- as.character(seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2                <- unlist(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)
            names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)         <- paste0("hetero.dr.adjFGlmInnd.propscFGlmInnd.cv2.", names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2))
            print("hetero.F.F")
            
            performance.hetero.drInnd <- c(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2, 
                                           eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2,
                                           eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2,
                                           eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2,
                                           eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)
            
            #####################################################################################################################
            ################################################### Homogeneous #####################################################
            #####################################################################################################################
            data.cont.cont            <- makeData.cont.noeff.cont(1000, 1000, coeff.prop.sc = 0.6)
            data.used.full.cont.cont  <- data.cont.cont$data.used
            data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
            # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
            data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ] 
            
            #####################################################################################################################
            ###################### 6. dr: True GLM Model in node, True propensity score model in node, Cv2 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
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
            
            final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                               tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2$tree.list,
                                                                               type.var          = "cont",
                                                                               seed              = a[i], 
                                                                               n.cv              = 5,
                                                                               propsc.mod.loc    = "node",
                                                                               propsc.mthd       = "GLM",
                                                                               propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                               adj.mod.loc       = "node",
                                                                               adj.mthd          = "GLM",
                                                                               adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)")
            t1 <- Sys.time()
            
            
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2[[1]],
                                                                                 test.data    = data.cont.cont$test.data,
                                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                 noise.var    = data.cont.cont$noise.var,
                                                                                 corr.split   = data.cont.cont$corr.split,
                                                                                 where.split  = data.cont.cont$where.split,
                                                                                 dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2                <- unlist(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2)
            names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2)         <- paste0("homo.dr.adjTGlmInnd.propscTGlmInnd.cv2.", names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2))
            print("homo.T.T")
            
            #####################################################################################################################
            #################### 7. dr: True GLM Model in node, Mis func propensity score model in node, Cv1 ####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
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
            
            final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                                  tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$tree.list,
                                                                                  type.var          = "cont",
                                                                                  seed              = a[i], 
                                                                                  n.cv              = 5,
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                  adj.mod.loc       = "node",
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2[[1]],
                                                                                    test.data    = data.cont.cont$test.data,
                                                                                    true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                    noise.var    = data.cont.cont$noise.var,
                                                                                    corr.split   = data.cont.cont$corr.split,
                                                                                    where.split  = data.cont.cont$where.split,
                                                                                    dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2                <- unlist(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2)
            names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2)         <- paste0("homo.dr.adjTGlmInnd.propscNoisGlmInnd.cv2.", names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2))
            print("homo.T.Nois")
            
            #####################################################################################################################
            ###################### 8. dr: Noisy GLM Model in node, True propensity score model in node, Cv2 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
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
            
            final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                                  tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$tree.list,
                                                                                  type.var          = "cont",
                                                                                  seed              = a[i], 
                                                                                  n.cv              = 5,
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ X1 + X2 + X3",
                                                                                  adj.mod.loc       = "node", 
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2[[1]],
                                                                                    test.data    = data.cont.cont$test.data,
                                                                                    true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                    noise.var    = data.cont.cont$noise.var,
                                                                                    corr.split   = data.cont.cont$corr.split,
                                                                                    where.split  = data.cont.cont$where.split,
                                                                                    dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2)
            names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2)         <- paste0("homo.dr.adjNoisGlmInnd.propscTGlmInnd.cv2.", names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2))
            print("homo.Nois.T")
            
            #####################################################################################################################
            ##################### 9. dr: Noisy GLM Model in node, Noisy propensity score model in node, Cv2 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
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
            
            final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                                     tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$tree.list,
                                                                                     type.var          = "cont",
                                                                                     seed              = a[i], 
                                                                                     n.cv              = 5,
                                                                                     propsc.mod.loc    = "node",
                                                                                     propsc.mthd       = "GLM",
                                                                                     propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + exp(X4) + exp(X5) + exp(X6)",
                                                                                     adj.mod.loc       = "node",
                                                                                     adj.mthd          = "GLM",
                                                                                     adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2[[1]],
                                                                                       test.data    = data.cont.cont$test.data,
                                                                                       true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                       noise.var    = data.cont.cont$noise.var,
                                                                                       corr.split   = data.cont.cont$corr.split,
                                                                                       where.split  = data.cont.cont$where.split,
                                                                                       dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2)
            names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2)         <- paste0("homo.dr.adjNoisGlmInnd.propscNoisGlmInnd.cv2.", names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2))
            print("homo.Nois.Nois")
            
            #####################################################################################################################
            ############## 10. dr: Misspecified GLM Model in node, Misspecified propensity score model in node, Cv2 #############
            #####################################################################################################################
            data.used.full.cont.cont.mis <- data.used.full.cont.cont %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
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
            
            final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                                               tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2$tree.list,
                                                                               type.var          = "cont",
                                                                               seed              = a[i], 
                                                                               n.cv              = 5,
                                                                               propsc.mod.loc    = "node",
                                                                               propsc.mthd       = "GLM",
                                                                               propsc.form.true  = NULL,
                                                                               adj.mod.loc       = "node",
                                                                               adj.mthd          = "GLM",
                                                                               adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv2[[1]],
                                                                                 test.data    = data.cont.cont$test.data,
                                                                                 true.trt.eff = data.cont.cont$true.trt.eff,
                                                                                 noise.var    = data.cont.cont$noise.var,
                                                                                 corr.split   = data.cont.cont$corr.split,
                                                                                 where.split  = data.cont.cont$where.split,
                                                                                 dir.split    = data.cont.cont$dir.split)
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2                <- unlist(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)
            names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)         <- paste0("homo.dr.adjFGlmInnd.propscFGlmInnd.cv2.", names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2))
            print("homo.F.F")
            
            performance.homo.drInnd <- c(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2, 
                                         eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2,
                                         eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2,
                                         eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2,
                                         eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)
            
            # print must be put before output, otherwise the output will be print output
            if (i%%30 == 0) {print(i)}
            
            c(performance.hetero.drInnd, performance.homo.drInnd)

          }

save(performance.drInnd, file = "../Data/AppendixC4/Cv2Dr.RData")

