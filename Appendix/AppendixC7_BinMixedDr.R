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
            
            data.bin.mixed            <- makeData.bin.eff.mixed(N             = 1000, 
                                                                n.test        = 1000, 
                                                                p.cont        = 3, 
                                                                p.cate        = 3, 
                                                                n.cate        = 4:6, 
                                                                coeff.prop.sc = 0.3,
                                                                seed          = a[i])
            data.used.full.bin.mixed  <- data.bin.mixed$data.used
            data.used.bin.mixed       <- data.used.full.bin.mixed[1:800, ]
            # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
            data.validation.bin.mixed <- data.used.full.bin.mixed[801:1000, ]  
            
            # data.used.bin.mixed       <- data.used.full.bin.mixed[1:400, ]
            # # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
            # data.validation.bin.mixed <- data.used.full.bin.mixed[401:500, ]  
            
            #####################################################################################################################
            ###################### 1. dr: True GLM Model in node, True propensity score model in node, Cv1 ######################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                                est.used          = "DR",
                                                                                type.var          = "bin",
                                                                                propsc.mod.loc    = "node",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                                adj.mod.loc       = "node", 
                                                                                adj.mthd          = "GLM", 
                                                                                adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))", 
                                                                                num.truc.obs      = 30,
                                                                                min.node          = 20)
            
            final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                               tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                               lambda.used       = qchisq(0.95, 1),
                                                                               val.sample        = data.validation.bin.mixed,
                                                                               type.var          = "bin",
                                                                               propsc.mod.loc    = "node",
                                                                               propsc.mthd       = "GLM",
                                                                               propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                               adj.mod.loc       = "node",
                                                                               adj.mthd          = "GLM",
                                                                               adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1[[1]],
                                                                                 test.data    = data.bin.mixed$test.data,
                                                                                 true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                 noise.var    = data.bin.mixed$noise.var,
                                                                                 corr.split   = data.bin.mixed$corr.split,
                                                                                 where.split  = data.bin.mixed$where.split,
                                                                                 dir.split    = data.bin.mixed$dir.split,
                                                                                 split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1$tree.list[[1]],
                                                                                                       corr.split = data.bin.mixed$corr.split,
                                                                                                       split.cate = data.bin.mixed$split.cate)
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1                <- unlist(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1)
            names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1)         <- paste0("hetero.dr.adjTGlmInnd.propscTGlmInnd.cv1.", names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1))
            print("1.cv1.T.T")
            
            #####################################################################################################################
            #################### 2. dr: True GLM Model in node, Mis func propensity score model in node, Cv1 ####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                                   est.used          = "DR",
                                                                                   type.var          = "bin",
                                                                                   propsc.mod.loc    = "node",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                   adj.mod.loc       = "node", 
                                                                                   adj.mthd          = "GLM", 
                                                                                   adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))", 
                                                                                   num.truc.obs      = 30,
                                                                                   min.node          = 20)
            
            final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                                  tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                                  lambda.used       = qchisq(0.95, 1),
                                                                                  val.sample        = data.validation.bin.mixed,
                                                                                  type.var          = "bin",
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                  adj.mod.loc       = "node",
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1[[1]],
                                                                                    test.data    = data.bin.mixed$test.data,
                                                                                    true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                    noise.var    = data.bin.mixed$noise.var,
                                                                                    corr.split   = data.bin.mixed$corr.split,
                                                                                    where.split  = data.bin.mixed$where.split,
                                                                                    dir.split    = data.bin.mixed$dir.split,
                                                                                    split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$tree.list[[1]],
                                                                                                          corr.split = data.bin.mixed$corr.split,
                                                                                                          split.cate = data.bin.mixed$split.cate)
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1                <- unlist(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1)
            names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1)         <- paste0("hetero.dr.adjTGlmInnd.propscNoisGlmInnd.cv1.", names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1))
            print("2.cv1.T.Nois")
            
            #####################################################################################################################
            ##################### 3. dr: Mis func GLM Model in node, True propensity score model in node, Cv1 ###################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                                   est.used          = "DR",
                                                                                   type.var          = "bin",
                                                                                   propsc.mod.loc    = "node",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                                   adj.mod.loc       = "node", 
                                                                                   adj.mthd          = "GLM", 
                                                                                   adj.form.true     = NULL, 
                                                                                   num.truc.obs      = 30,
                                                                                   min.node          = 20)
            
            final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                                  tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                                  lambda.used       = qchisq(0.95, 1),
                                                                                  val.sample        = data.validation.bin.mixed,
                                                                                  type.var          = "bin",
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                                  adj.mod.loc       = "node",
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1[[1]],
                                                                                    test.data    = data.bin.mixed$test.data,
                                                                                    true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                    noise.var    = data.bin.mixed$noise.var,
                                                                                    corr.split   = data.bin.mixed$corr.split,
                                                                                    where.split  = data.bin.mixed$where.split,
                                                                                    dir.split    = data.bin.mixed$dir.split,
                                                                                    split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$tree.list[[1]],
                                                                                                          corr.split = data.bin.mixed$corr.split,
                                                                                                          split.cate = data.bin.mixed$split.cate)
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1)
            names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1)         <- paste0("hetero.dr.adjNoisGlmInnd.propscTGlmInnd.cv1.", names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1))
            print("3.cv1.Nois.T")
            
            #####################################################################################################################
            ################## 4. dr: Mis func GLM Model in node, Mis func propensity score model in node, Cv1 ##################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                                      est.used          = "DR",
                                                                                      type.var          = "bin",
                                                                                      propsc.mod.loc    = "node",
                                                                                      propsc.mthd       = "GLM",
                                                                                      propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                      adj.mod.loc       = "node", 
                                                                                      adj.mthd          = "GLM", 
                                                                                      adj.form.true     = NULL, 
                                                                                      num.truc.obs      = 30,
                                                                                      min.node          = 20)
            
            final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                                     tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                                     lambda.used       = qchisq(0.95, 1),
                                                                                     val.sample        = data.validation.bin.mixed,
                                                                                     type.var          = "bin",
                                                                                     propsc.mod.loc    = "node",
                                                                                     propsc.mthd       = "GLM",
                                                                                     propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                     adj.mod.loc       = "node",
                                                                                     adj.mthd          = "GLM",
                                                                                     adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1[[1]],
                                                                                       test.data    = data.bin.mixed$test.data,
                                                                                       true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                       noise.var    = data.bin.mixed$noise.var,
                                                                                       corr.split   = data.bin.mixed$corr.split,
                                                                                       where.split  = data.bin.mixed$where.split,
                                                                                       dir.split    = data.bin.mixed$dir.split,
                                                                                       split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$tree.list[[1]],
                                                                                                             corr.split = data.bin.mixed$corr.split,
                                                                                                             split.cate = data.bin.mixed$split.cate)
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1)
            names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1)         <- paste0("hetero.dr.adjNoisGlmInnd.propscNoisGlmInnd.cv1.", names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1))
            print("4.cv1.nois.nois")
            
            #####################################################################################################################
            ############ 5. dr: Unmeasured cov GLM Model in node, Unmeasured cov propensity score model in node, Cv1 ############
            #####################################################################################################################
            data.used.bin.mixed.mis <- data.used.bin.mixed %>%
              select(-X2)
            data.validation.bin.mixed.mis <- data.validation.bin.mixed %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed.mis,
                                                                                est.used          = "DR",
                                                                                type.var          = "bin",
                                                                                propsc.mod.loc    = "node",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = NULL,
                                                                                adj.mod.loc       = "node", 
                                                                                adj.mthd          = "GLM", 
                                                                                adj.form.true     = NULL, 
                                                                                num.truc.obs      = 30,
                                                                                min.node          = 20)
            
            final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed.mis,
                                                                               tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1$tree.list, 
                                                                               lambda.used       = qchisq(0.95, 1),
                                                                               val.sample        = data.validation.bin.mixed.mis,
                                                                               type.var          = "bin",
                                                                               propsc.mod.loc    = "node",
                                                                               propsc.mthd       = "GLM",
                                                                               propsc.form.true  = NULL,
                                                                               adj.mod.loc       = "node",
                                                                               adj.mthd          = "GLM",
                                                                               adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1[[1]],
                                                                                 test.data    = data.bin.mixed$test.data,
                                                                                 true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                 noise.var    = data.bin.mixed$noise.var,
                                                                                 corr.split   = data.bin.mixed$corr.split,
                                                                                 where.split  = data.bin.mixed$where.split,
                                                                                 dir.split    = data.bin.mixed$dir.split,
                                                                                 split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1$tree.list[[1]],
                                                                                                       corr.split = data.bin.mixed$corr.split,
                                                                                                       split.cate = data.bin.mixed$split.cate)
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1                <- unlist(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1)
            names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1)         <- paste0("hetero.dr.adjFGlmInnd.propscFGlmInnd.cv1.", names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1))
            print("5.cv1.F.F")
            
            #####################################################################################################################
            ####################### 6. dr: True GLM Model in node, True propensity score model in node, Cv2 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                                est.used          = "DR",
                                                                                type.var          = "bin",
                                                                                propsc.mod.loc    = "node",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                                adj.mod.loc       = "node", 
                                                                                adj.mthd          = "GLM", 
                                                                                adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))", 
                                                                                num.truc.obs      = 30,
                                                                                min.node          = 20)
            
            final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                               tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2$tree.list,
                                                                               type.var          = "bin",
                                                                               seed              = a[i], 
                                                                               n.cv              = 5,
                                                                               propsc.mod.loc    = "node",
                                                                               propsc.mthd       = "GLM",
                                                                               propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                               adj.mod.loc       = "node",
                                                                               adj.mthd          = "GLM",
                                                                               adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2[[1]],
                                                                                 test.data    = data.bin.mixed$test.data,
                                                                                 true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                 noise.var    = data.bin.mixed$noise.var,
                                                                                 corr.split   = data.bin.mixed$corr.split,
                                                                                 where.split  = data.bin.mixed$where.split,
                                                                                 dir.split    = data.bin.mixed$dir.split,
                                                                                 split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2$tree.list[[1]],
                                                                                                       corr.split = data.bin.mixed$corr.split,
                                                                                                       split.cate = data.bin.mixed$split.cate)
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2                <- unlist(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2)
            names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2)         <- paste0("hetero.dr.adjTGlmInnd.propscTGlmInnd.cv2.", names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2))
            print("6.cv2.T.T")
            
            #####################################################################################################################
            #################### 7. dr: True GLM Model in node, Mis func propensity score model in node, Cv1 ####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                                   est.used          = "DR",
                                                                                   type.var          = "bin",
                                                                                   propsc.mod.loc    = "node",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                   adj.mod.loc       = "node", 
                                                                                   adj.mthd          = "GLM", 
                                                                                   adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))", 
                                                                                   num.truc.obs      = 30,
                                                                                   min.node          = 20)
            
            final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                                  tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$tree.list,
                                                                                  type.var          = "bin",
                                                                                  seed              = a[i], 
                                                                                  n.cv              = 5,
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                  adj.mod.loc       = "node",
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2[[1]],
                                                                                    test.data    = data.bin.mixed$test.data,
                                                                                    true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                    noise.var    = data.bin.mixed$noise.var,
                                                                                    corr.split   = data.bin.mixed$corr.split,
                                                                                    where.split  = data.bin.mixed$where.split,
                                                                                    dir.split    = data.bin.mixed$dir.split,
                                                                                    split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$tree.list[[1]],
                                                                                                          corr.split = data.bin.mixed$corr.split,
                                                                                                          split.cate = data.bin.mixed$split.cate)
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2                <- unlist(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2)
            names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2)         <- paste0("hetero.dr.adjTGlmInnd.propscNoisGlmInnd.cv2.", names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2))
            print("7.cv2.T.Nois")
            
            #####################################################################################################################
            ###################### 8. dr: Noisy GLM Model in node, True propensity score model in node, Cv2 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                                   est.used          = "DR",
                                                                                   type.var          = "bin",
                                                                                   propsc.mod.loc    = "node",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                                   adj.mod.loc       = "node", 
                                                                                   adj.mthd          = "GLM", 
                                                                                   adj.form.true     = NULL, 
                                                                                   num.truc.obs      = 30,
                                                                                   min.node          = 20)
            
            final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                                  tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$tree.list,
                                                                                  type.var          = "bin",
                                                                                  seed              = a[i], 
                                                                                  n.cv              = 5,
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                                  adj.mod.loc       = "node", 
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2[[1]],
                                                                                    test.data    = data.bin.mixed$test.data,
                                                                                    true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                    noise.var    = data.bin.mixed$noise.var,
                                                                                    corr.split   = data.bin.mixed$corr.split,
                                                                                    where.split  = data.bin.mixed$where.split,
                                                                                    dir.split    = data.bin.mixed$dir.split,
                                                                                    split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$tree.list[[1]],
                                                                                                          corr.split = data.bin.mixed$corr.split,
                                                                                                          split.cate = data.bin.mixed$split.cate)
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2)
            names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2)         <- paste0("hetero.dr.adjNoisGlmInnd.propscTGlmInnd.cv2.", names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2))
            print("8.cv2.Nois.T")
            
            #####################################################################################################################
            ##################### 9. dr: Noisy GLM Model in node, Noisy propensity score model in node, Cv2 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                                      est.used          = "DR",
                                                                                      type.var          = "bin",
                                                                                      propsc.mod.loc    = "node",
                                                                                      propsc.mthd       = "GLM",
                                                                                      propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                      adj.mod.loc       = "node", 
                                                                                      adj.mthd          = "GLM", 
                                                                                      adj.form.true     = NULL, 
                                                                                      num.truc.obs      = 30,
                                                                                      min.node          = 20)
            
            final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                                     tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$tree.list,
                                                                                     type.var          = "bin",
                                                                                     seed              = a[i], 
                                                                                     n.cv              = 5,
                                                                                     propsc.mod.loc    = "node",
                                                                                     propsc.mthd       = "GLM",
                                                                                     propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                     adj.mod.loc       = "node",
                                                                                     adj.mthd          = "GLM",
                                                                                     adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2[[1]],
                                                                                       test.data    = data.bin.mixed$test.data,
                                                                                       true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                       noise.var    = data.bin.mixed$noise.var,
                                                                                       corr.split   = data.bin.mixed$corr.split,
                                                                                       where.split  = data.bin.mixed$where.split,
                                                                                       dir.split    = data.bin.mixed$dir.split,
                                                                                       split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$tree.list[[1]],
                                                                                                             corr.split = data.bin.mixed$corr.split,
                                                                                                             split.cate = data.bin.mixed$split.cate)
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2)
            names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2)         <- paste0("hetero.dr.adjNoisGlmInnd.propscNoisGlmInnd.cv2.", names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2))
            print("9.cv2.Nois.Nois")
            
            #####################################################################################################################
            ############## 10. dr: Misspecified GLM Model in node, Misspecified propensity score model in node, Cv2 #############
            #####################################################################################################################
            data.used.full.bin.mixed.mis <- data.used.full.bin.mixed %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed.mis,
                                                                                est.used          = "DR",
                                                                                type.var          = "bin",
                                                                                propsc.mod.loc    = "node",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = NULL,
                                                                                adj.mod.loc       = "node", 
                                                                                adj.mthd          = "GLM", 
                                                                                adj.form.true     = NULL, 
                                                                                num.truc.obs      = 30,
                                                                                min.node          = 20)
            
            final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed.mis,
                                                                               tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2$tree.list,
                                                                               type.var          = "bin",
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
                                                                                 test.data    = data.bin.mixed$test.data,
                                                                                 true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                 noise.var    = data.bin.mixed$noise.var,
                                                                                 corr.split   = data.bin.mixed$corr.split,
                                                                                 where.split  = data.bin.mixed$where.split,
                                                                                 dir.split    = data.bin.mixed$dir.split,
                                                                                 split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2$tree.list[[1]],
                                                                                                       corr.split = data.bin.mixed$corr.split,
                                                                                                       split.cate = data.bin.mixed$split.cate)
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2                <- unlist(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)
            names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)         <- paste0("hetero.dr.adjFGlmInnd.propscFGlmInnd.cv2.", names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2))
            print("10.cv2.F.F")
            
            performance.hetero.drInnd <- c(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1, 
                                           eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1,
                                           eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1,
                                           eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1,
                                           eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1,
                                           eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2, 
                                           eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2,
                                           eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2,
                                           eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2,
                                           eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)
            
            
            
            #####################################################################################################################
            ################################################### Homogeneous #####################################################
            #####################################################################################################################
            data.bin.mixed            <- makeData.bin.noeff.mixed(N             = 1000, 
                                                                  n.test        = 1000, 
                                                                  p.cont        = 3, 
                                                                  p.cate        = 3, 
                                                                  n.cate        = 4:6, 
                                                                  coeff.prop.sc = 0.3, 
                                                                  seed          = a[i])
            data.used.full.bin.mixed  <- data.bin.mixed$data.used
            data.used.bin.mixed       <- data.used.full.bin.mixed[1:800, ]
            # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
            data.validation.bin.mixed <- data.used.full.bin.mixed[801:1000, ]  
            
            #####################################################################################################################
            ###################### 11. dr: True GLM Model in node, True propensity score model in node, Cv1 ######################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                                est.used          = "DR",
                                                                                type.var          = "bin",
                                                                                propsc.mod.loc    = "node",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                                adj.mod.loc       = "node", 
                                                                                adj.mthd          = "GLM", 
                                                                                adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))", 
                                                                                num.truc.obs      = 30,
                                                                                min.node          = 20)
            
            final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                               tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                               lambda.used       = qchisq(0.95, 1),
                                                                               val.sample        = data.validation.bin.mixed,
                                                                               type.var          = "bin",
                                                                               propsc.mod.loc    = "node",
                                                                               propsc.mthd       = "GLM",
                                                                               propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                               adj.mod.loc       = "node",
                                                                               adj.mthd          = "GLM",
                                                                               adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv1[[1]],
                                                                                 test.data    = data.bin.mixed$test.data,
                                                                                 true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                 noise.var    = data.bin.mixed$noise.var,
                                                                                 corr.split   = data.bin.mixed$corr.split,
                                                                                 where.split  = data.bin.mixed$where.split,
                                                                                 dir.split    = data.bin.mixed$dir.split,
                                                                                 split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1                <- unlist(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1)
            names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1)         <- paste0("homo.dr.adjTGlmInnd.propscTGlmInnd.cv1.", names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1))
            print("11.cv1.T.T")
            
            #####################################################################################################################
            #################### 12. dr: True GLM Model in node, Mis func propensity score model in node, Cv1 ###################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                                   est.used          = "DR",
                                                                                   type.var          = "bin",
                                                                                   propsc.mod.loc    = "node",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                   adj.mod.loc       = "node", 
                                                                                   adj.mthd          = "GLM", 
                                                                                   adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))", 
                                                                                   num.truc.obs      = 30,
                                                                                   min.node          = 20)
            
            final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                                  tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                                  lambda.used       = qchisq(0.95, 1),
                                                                                  val.sample        = data.validation.bin.mixed,
                                                                                  type.var          = "bin",
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                  adj.mod.loc       = "node",
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1[[1]],
                                                                                    test.data    = data.bin.mixed$test.data,
                                                                                    true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                    noise.var    = data.bin.mixed$noise.var,
                                                                                    corr.split   = data.bin.mixed$corr.split,
                                                                                    where.split  = data.bin.mixed$where.split,
                                                                                    dir.split    = data.bin.mixed$dir.split,
                                                                                    split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1                <- unlist(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1)
            names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1)         <- paste0("homo.dr.adjTGlmInnd.propscNoisGlmInnd.cv1.", names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1))
            print("12.cv1.T.Nois")
            
            #####################################################################################################################
            ##################### 13. dr: Mis func GLM Model in node, True propensity score model in node, Cv1 ##################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                                   est.used          = "DR",
                                                                                   type.var          = "bin",
                                                                                   propsc.mod.loc    = "node",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                                   adj.mod.loc       = "node", 
                                                                                   adj.mthd          = "GLM", 
                                                                                   adj.form.true     = NULL, 
                                                                                   num.truc.obs      = 30,
                                                                                   min.node          = 20)
            
            final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                                  tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$tree.list, 
                                                                                  lambda.used       = qchisq(0.95, 1),
                                                                                  val.sample        = data.validation.bin.mixed,
                                                                                  type.var          = "bin",
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                                  adj.mod.loc       = "node",
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1[[1]],
                                                                                    test.data    = data.bin.mixed$test.data,
                                                                                    true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                    noise.var    = data.bin.mixed$noise.var,
                                                                                    corr.split   = data.bin.mixed$corr.split,
                                                                                    where.split  = data.bin.mixed$where.split,
                                                                                    dir.split    = data.bin.mixed$dir.split,
                                                                                    split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1)
            names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1)         <- paste0("homo.dr.adjNoisGlmInnd.propscTGlmInnd.cv1.", names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1))
            print("13.cv1.Nois.T")
            
            #####################################################################################################################
            ################## 14. dr: Mis func GLM Model in node, Mis func propensity score model in node, Cv1 #################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                                      est.used          = "DR",
                                                                                      type.var          = "bin",
                                                                                      propsc.mod.loc    = "node",
                                                                                      propsc.mthd       = "GLM",
                                                                                      propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                      adj.mod.loc       = "node", 
                                                                                      adj.mthd          = "GLM", 
                                                                                      adj.form.true     = NULL, 
                                                                                      num.truc.obs      = 30,
                                                                                      min.node          = 20)
            
            final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed,
                                                                                     tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$tree.list, 
                                                                                     lambda.used       = qchisq(0.95, 1),
                                                                                     val.sample        = data.validation.bin.mixed,
                                                                                     type.var          = "bin",
                                                                                     propsc.mod.loc    = "node",
                                                                                     propsc.mthd       = "GLM",
                                                                                     propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                     adj.mod.loc       = "node",
                                                                                     adj.mthd          = "GLM",
                                                                                     adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1[[1]],
                                                                                       test.data    = data.bin.mixed$test.data,
                                                                                       true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                       noise.var    = data.bin.mixed$noise.var,
                                                                                       corr.split   = data.bin.mixed$corr.split,
                                                                                       where.split  = data.bin.mixed$where.split,
                                                                                       dir.split    = data.bin.mixed$dir.split,
                                                                                       split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1)
            names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1)         <- paste0("homo.dr.adjNoisGlmInnd.propscNoisGlmInnd.cv1.", names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1))
            print("14.cv1.nois.nois")
            
            #####################################################################################################################
            ############ 15. dr: Unmeasured cov GLM Model in node, Unmeasured cov propensity score model in node, Cv1 ###########
            #####################################################################################################################
            data.used.bin.mixed.mis <- data.used.bin.mixed %>%
              select(-X2)
            data.validation.bin.mixed.mis <- data.validation.bin.mixed %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- create.sequence(data.used         = data.used.bin.mixed.mis,
                                                                                est.used          = "DR",
                                                                                type.var          = "bin",
                                                                                propsc.mod.loc    = "node",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = NULL,
                                                                                adj.mod.loc       = "node", 
                                                                                adj.mthd          = "GLM", 
                                                                                adj.form.true     = NULL, 
                                                                                num.truc.obs      = 30,
                                                                                min.node          = 20)
            
            final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- EstDr.CvMethod1(data.used         = data.used.bin.mixed.mis,
                                                                               tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv1$tree.list, 
                                                                               lambda.used       = qchisq(0.95, 1),
                                                                               val.sample        = data.validation.bin.mixed.mis,
                                                                               type.var          = "bin",
                                                                               propsc.mod.loc    = "node",
                                                                               propsc.mthd       = "GLM",
                                                                               propsc.form.true  = NULL,
                                                                               adj.mod.loc       = "node",
                                                                               adj.mthd          = "GLM",
                                                                               adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1 <- eval.measures.eff(final.tree   = final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv1[[1]],
                                                                                 test.data    = data.bin.mixed$test.data,
                                                                                 true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                 noise.var    = data.bin.mixed$noise.var,
                                                                                 corr.split   = data.bin.mixed$corr.split,
                                                                                 where.split  = data.bin.mixed$where.split,
                                                                                 dir.split    = data.bin.mixed$dir.split,
                                                                                 split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1                <- unlist(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1)
            names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1)         <- paste0("homo.dr.adjFGlmInnd.propscFGlmInnd.cv1.", names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1))
            print("15.cv1.F.F")
            
            #####################################################################################################################
            ###################### 16. dr: True GLM Model in node, True propensity score model in node, Cv2 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                                est.used          = "DR",
                                                                                type.var          = "bin",
                                                                                propsc.mod.loc    = "node",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                                adj.mod.loc       = "node", 
                                                                                adj.mthd          = "GLM", 
                                                                                adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))", 
                                                                                num.truc.obs      = 30,
                                                                                min.node          = 20)
            
            final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                               tree.list         = seq.created.estdr.adjTGlmInnd.propscTGlmInnd.cv2$tree.list,
                                                                               type.var          = "bin",
                                                                               seed              = a[i], 
                                                                               n.cv              = 5,
                                                                               propsc.mod.loc    = "node",
                                                                               propsc.mthd       = "GLM",
                                                                               propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                               adj.mod.loc       = "node",
                                                                               adj.mthd          = "GLM",
                                                                               adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscTGlmInnd.cv2[[1]],
                                                                                 test.data    = data.bin.mixed$test.data,
                                                                                 true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                 noise.var    = data.bin.mixed$noise.var,
                                                                                 corr.split   = data.bin.mixed$corr.split,
                                                                                 where.split  = data.bin.mixed$where.split,
                                                                                 dir.split    = data.bin.mixed$dir.split,
                                                                                 split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2                <- unlist(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2)
            names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2)         <- paste0("homo.dr.adjTGlmInnd.propscTGlmInnd.cv2.", names(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2))
            print("16.cv2.T.T")
            
            #####################################################################################################################
            #################### 17. dr: True GLM Model in node, Mis func propensity score model in node, Cv1 ###################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                                   est.used          = "DR",
                                                                                   type.var          = "bin",
                                                                                   propsc.mod.loc    = "node",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                   adj.mod.loc       = "node", 
                                                                                   adj.mthd          = "GLM", 
                                                                                   adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))", 
                                                                                   num.truc.obs      = 30,
                                                                                   min.node          = 20)
            
            final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                                  tree.list         = seq.created.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$tree.list,
                                                                                  type.var          = "bin",
                                                                                  seed              = a[i], 
                                                                                  n.cv              = 5,
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                  adj.mod.loc       = "node",
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))")
            t1 <- Sys.time()
            
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2[[1]],
                                                                                    test.data    = data.bin.mixed$test.data,
                                                                                    true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                    noise.var    = data.bin.mixed$noise.var,
                                                                                    corr.split   = data.bin.mixed$corr.split,
                                                                                    where.split  = data.bin.mixed$where.split,
                                                                                    dir.split    = data.bin.mixed$dir.split,
                                                                                    split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2                <- unlist(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2)
            names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2)         <- paste0("homo.dr.adjTGlmInnd.propscNoisGlmInnd.cv2.", names(eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2))
            print("17.cv2.T.Nois")
            
            #####################################################################################################################
            ###################### 18. dr: Noisy GLM Model in node, True propensity score model in node, Cv2 ####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                                   est.used          = "DR",
                                                                                   type.var          = "bin",
                                                                                   propsc.mod.loc    = "node",
                                                                                   propsc.mthd       = "GLM",
                                                                                   propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                                   adj.mod.loc       = "node", 
                                                                                   adj.mthd          = "GLM", 
                                                                                   adj.form.true     = NULL, 
                                                                                   num.truc.obs      = 30,
                                                                                   min.node          = 20)
            
            final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                                  tree.list         = seq.created.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$tree.list,
                                                                                  type.var          = "bin",
                                                                                  seed              = a[i], 
                                                                                  n.cv              = 5,
                                                                                  propsc.mod.loc    = "node",
                                                                                  propsc.mthd       = "GLM",
                                                                                  propsc.form.true  = "A ~ X2 + X3 + ((X6 == 'B') | (X6 == 'C'))",
                                                                                  adj.mod.loc       = "node", 
                                                                                  adj.mthd          = "GLM",
                                                                                  adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2[[1]],
                                                                                    test.data    = data.bin.mixed$test.data,
                                                                                    true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                    noise.var    = data.bin.mixed$noise.var,
                                                                                    corr.split   = data.bin.mixed$corr.split,
                                                                                    where.split  = data.bin.mixed$where.split,
                                                                                    dir.split    = data.bin.mixed$dir.split,
                                                                                    split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2)
            names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2)         <- paste0("homo.dr.adjNoisGlmInnd.propscTGlmInnd.cv2.", names(eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2))
            print("18.cv2.Nois.T")
            
            #####################################################################################################################
            #################### 19. dr: Noisy GLM Model in node, Noisy propensity score model in node, Cv2 #####################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                                      est.used          = "DR",
                                                                                      type.var          = "bin",
                                                                                      propsc.mod.loc    = "node",
                                                                                      propsc.mthd       = "GLM",
                                                                                      propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                      adj.mod.loc       = "node", 
                                                                                      adj.mthd          = "GLM", 
                                                                                      adj.form.true     = NULL, 
                                                                                      num.truc.obs      = 30,
                                                                                      min.node          = 20)
            
            final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                                     tree.list         = seq.created.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$tree.list,
                                                                                     type.var          = "bin",
                                                                                     seed              = a[i], 
                                                                                     n.cv              = 5,
                                                                                     propsc.mod.loc    = "node",
                                                                                     propsc.mthd       = "GLM",
                                                                                     propsc.form.true  = "A ~ exp(X1) + exp(X2) + exp(X3) + X4 + X5 + X6",
                                                                                     adj.mod.loc       = "node",
                                                                                     adj.mthd          = "GLM",
                                                                                     adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2 <- eval.measures.eff(final.tree   = final.tree.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2[[1]],
                                                                                       test.data    = data.bin.mixed$test.data,
                                                                                       true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                       noise.var    = data.bin.mixed$noise.var,
                                                                                       corr.split   = data.bin.mixed$corr.split,
                                                                                       where.split  = data.bin.mixed$where.split,
                                                                                       dir.split    = data.bin.mixed$dir.split,
                                                                                       split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2                <- unlist(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2)
            names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2)         <- paste0("homo.dr.adjNoisGlmInnd.propscNoisGlmInnd.cv2.", names(eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2))
            print("19.cv2.Nois.Nois")
            
            #####################################################################################################################
            ############## 20. dr: Misspecified GLM Model in node, Misspecified propensity score model in node, Cv2 #############
            #####################################################################################################################
            data.used.full.bin.mixed.mis <- data.used.full.bin.mixed %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed.mis,
                                                                                est.used          = "DR",
                                                                                type.var          = "bin",
                                                                                propsc.mod.loc    = "node",
                                                                                propsc.mthd       = "GLM",
                                                                                propsc.form.true  = NULL,
                                                                                adj.mod.loc       = "node", 
                                                                                adj.mthd          = "GLM", 
                                                                                adj.form.true     = NULL, 
                                                                                num.truc.obs      = 30,
                                                                                min.node          = 20)
            
            final.tree.estdr.adjFGlmInnd.propscFGlmInnd.cv2 <- EstDr.CvMethod2(data.used         = data.used.full.bin.mixed.mis,
                                                                               tree.list         = seq.created.estdr.adjFGlmInnd.propscFGlmInnd.cv2$tree.list,
                                                                               type.var          = "bin",
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
                                                                                 test.data    = data.bin.mixed$test.data,
                                                                                 true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                                 noise.var    = data.bin.mixed$noise.var,
                                                                                 corr.split   = data.bin.mixed$corr.split,
                                                                                 where.split  = data.bin.mixed$where.split,
                                                                                 dir.split    = data.bin.mixed$dir.split,
                                                                                 split.cate   = data.bin.mixed$split.cate)
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2                <- unlist(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)
            names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)         <- paste0("homo.dr.adjFGlmInnd.propscFGlmInnd.cv2.", names(eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2))
            print("20.cv2.F.F")
            
            performance.homo.drInnd <- c(eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv1, 
                                           eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv1,
                                           eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv1,
                                           eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv1,
                                           eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv1,
                                           eval.final.estdr.adjTGlmInnd.propscTGlmInnd.cv2, 
                                           eval.final.estdr.adjTGlmInnd.propscNoisGlmInnd.cv2,
                                           eval.final.estdr.adjNoisGlmInnd.propscTGlmInnd.cv2,
                                           eval.final.estdr.adjNoisGlmInnd.propscNoisGlmInnd.cv2,
                                           eval.final.estdr.adjFGlmInnd.propscFGlmInnd.cv2)
            
            # print must be put before output, otherwise the output will be print output
            if (i%%30 == 0) {print(i)}
            
            c(performance.hetero.drInnd, performance.homo.drInnd)

          }

save(performance.drInnd, file = paste0("../Data/AppendixC7/BinMixedDr.RData"))
