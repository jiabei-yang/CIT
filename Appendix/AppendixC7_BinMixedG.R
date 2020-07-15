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

performance.gInnd <- 
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
            
            #####################################################################################################################
            ############################# 1. g: GLM Model, inside node, True adjustment model, cv1 ##############################
            #####################################################################################################################
            t0 <- Sys.time()
            # The difference between true adjustment model and the misspecified model is very dependent on the true adjustment model
            seq.created.estg.glm.modinnd.true.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                     est.used          = "G",
                                                                     type.var          = "bin",
                                                                     adj.mod.loc       = "node",
                                                                     adj.mthd          = "GLM",
                                                                     adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))",
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)
            
            final.tree.estg.glm.modinnd.true.cv1 <- EstG.CvMethod1(data.used         = data.used.bin.mixed,
                                                                   tree.list         = seq.created.estg.glm.modinnd.true.cv1$tree.list,
                                                                   lambda.used       = qchisq(0.95, 1),
                                                                   val.sample        = data.validation.bin.mixed,              
                                                                   type.var          = "bin",
                                                                   adj.mod.loc       = "node",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))")
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.true.cv1[[1]],
                                                                      test.data    = data.bin.mixed$test.data,
                                                                      true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                      noise.var    = data.bin.mixed$noise.var,
                                                                      corr.split   = data.bin.mixed$corr.split,
                                                                      where.split  = data.bin.mixed$where.split,
                                                                      dir.split    = data.bin.mixed$dir.split,
                                                                      split.cate   = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.true.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estg.glm.modinnd.true.cv1$tree.list[[1]],
                                                                                            corr.split = data.bin.mixed$corr.split,
                                                                                            split.cate = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.true.cv1                <- unlist(eval.final.estg.glm.modinnd.true.cv1)
            names(eval.final.estg.glm.modinnd.true.cv1)         <- paste0("hetero.g.glm.modinnd.true.cv1.", names(eval.final.estg.glm.modinnd.true.cv1))
            print("1.cv1.T")
            
            #####################################################################################################################
            ############################ 2. g: GLM Model, inside node, Noisy adjustment model, cv1 ##############################
            #####################################################################################################################
            t0 <- Sys.time()
            # The difference between true adjustment model and the misspecified model is very dependent on the true adjustment model
            seq.created.estg.glm.modinnd.nois.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                     est.used          = "G",
                                                                     type.var          = "bin",
                                                                     adj.mod.loc       = "node",
                                                                     adj.mthd          = "GLM",
                                                                     adj.form.true     = NULL,
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)
            
            final.tree.estg.glm.modinnd.nois.cv1 <- EstG.CvMethod1(data.used         = data.used.bin.mixed,
                                                                   tree.list         = seq.created.estg.glm.modinnd.nois.cv1$tree.list,
                                                                   lambda.used       = qchisq(0.95, 1),
                                                                   val.sample        = data.validation.bin.mixed,              
                                                                   type.var          = "bin",
                                                                   adj.mod.loc       = "node",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.nois.cv1[[1]],
                                                                      test.data    = data.bin.mixed$test.data,
                                                                      true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                      noise.var    = data.bin.mixed$noise.var,
                                                                      corr.split   = data.bin.mixed$corr.split,
                                                                      where.split  = data.bin.mixed$where.split,
                                                                      dir.split    = data.bin.mixed$dir.split,
                                                                      split.cate   = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.nois.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estg.glm.modinnd.nois.cv1$tree.list[[1]],
                                                                                            corr.split = data.bin.mixed$corr.split,
                                                                                            split.cate = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.nois.cv1                <- unlist(eval.final.estg.glm.modinnd.nois.cv1)
            names(eval.final.estg.glm.modinnd.nois.cv1)         <- paste0("hetero.g.glm.modinnd.nois.cv1.", names(eval.final.estg.glm.modinnd.nois.cv1))
            print("2.cv1.Nois")
            
            #####################################################################################################################
            ########################## 3. g: GLM Model, inside node, Misspecified adjustment model, cv1 #########################
            #####################################################################################################################
            data.used.bin.mixed.mis <- data.used.bin.mixed %>%
              select(-X2)
            data.validation.bin.mixed.mis <- data.validation.bin.mixed %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estg.glm.modinnd.mis.cv1 <- create.sequence(data.used         = data.used.bin.mixed.mis,
                                                                    est.used          = "G",
                                                                    type.var          = "bin",
                                                                    adj.mod.loc       = "node",
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = NULL,
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)
            
            final.tree.estg.glm.modinnd.mis.cv1 <- EstG.CvMethod1(data.used         = data.used.bin.mixed.mis,
                                                                  tree.list         = seq.created.estg.glm.modinnd.mis.cv1$tree.list,
                                                                  lambda.used       = qchisq(0.95, 1),
                                                                  val.sample        = data.validation.bin.mixed.mis,              
                                                                  type.var          = "bin",
                                                                  adj.mod.loc       = "node",
                                                                  adj.mthd          = "GLM",
                                                                  adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.mis.cv1[[1]],
                                                                     test.data    = data.bin.mixed$test.data,
                                                                     true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                     noise.var    = data.bin.mixed$noise.var,
                                                                     corr.split   = data.bin.mixed$corr.split,
                                                                     where.split  = data.bin.mixed$where.split,
                                                                     dir.split    = data.bin.mixed$dir.split,
                                                                     split.cate   = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.mis.cv1$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estg.glm.modinnd.mis.cv1$tree.list[[1]],
                                                                                           corr.split = data.bin.mixed$corr.split,
                                                                                           split.cate = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.mis.cv1                <- unlist(eval.final.estg.glm.modinnd.mis.cv1)
            names(eval.final.estg.glm.modinnd.mis.cv1)         <- paste0("hetero.g.glm.modinnd.mis.cv1.", names(eval.final.estg.glm.modinnd.mis.cv1))
            print("3.cv1.Mis")
            
            #####################################################################################################################
            ############################# 4. g: GLM Model, inside node, True adjustment model, Cv2 ##############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modinnd.true.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                     est.used          = "G",
                                                                     type.var          = "bin",
                                                                     adj.mod.loc       = "node",
                                                                     adj.mthd          = "GLM",
                                                                     adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))",
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)
            
            final.tree.estg.glm.modinnd.true.cv2 <- EstG.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                   tree.list         = seq.created.estg.glm.modinnd.true.cv2$tree.list,             
                                                                   type.var          = "bin",
                                                                   seed              = a[i],
                                                                   n.cv              = 5,
                                                                   adj.mod.loc       = "node",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + X2 + A:((X4 == 'B') | (X4 == 'D'))")
            
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.true.cv2[[1]],
                                                                      test.data    = data.bin.mixed$test.data,
                                                                      true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                      noise.var    = data.bin.mixed$noise.var,
                                                                      corr.split   = data.bin.mixed$corr.split,
                                                                      where.split  = data.bin.mixed$where.split,
                                                                      dir.split    = data.bin.mixed$dir.split,
                                                                      split.cate   = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.true.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estg.glm.modinnd.true.cv2$tree.list[[1]],
                                                                                            corr.split = data.bin.mixed$corr.split,
                                                                                            split.cate = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.true.cv2                <- unlist(eval.final.estg.glm.modinnd.true.cv2)
            names(eval.final.estg.glm.modinnd.true.cv2)         <- paste0("hetero.g.glm.modinnd.true.cv2.", names(eval.final.estg.glm.modinnd.true.cv2))
            print("4.cv2.T")
            
            #####################################################################################################################
            ############################# 5. g: GLM Model, inside node, Noisy adjustment model, Cv2 #############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modinnd.nois.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                     est.used          = "G",
                                                                     type.var          = "bin",
                                                                     adj.mod.loc       = "node",
                                                                     adj.mthd          = "GLM",
                                                                     adj.form.true     = NULL,
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)
            
            final.tree.estg.glm.modinnd.nois.cv2 <- EstG.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                   tree.list         = seq.created.estg.glm.modinnd.nois.cv2$tree.list,             
                                                                   type.var          = "bin",
                                                                   seed              = a[i],
                                                                   n.cv              = 5,
                                                                   adj.mod.loc       = "node",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.nois.cv2[[1]],
                                                                      test.data    = data.bin.mixed$test.data,
                                                                      true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                      noise.var    = data.bin.mixed$noise.var,
                                                                      corr.split   = data.bin.mixed$corr.split,
                                                                      where.split  = data.bin.mixed$where.split,
                                                                      dir.split    = data.bin.mixed$dir.split,
                                                                      split.cate   = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.nois.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estg.glm.modinnd.nois.cv2$tree.list[[1]],
                                                                                            corr.split = data.bin.mixed$corr.split,
                                                                                            split.cate = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.nois.cv2                <- unlist(eval.final.estg.glm.modinnd.nois.cv2)
            names(eval.final.estg.glm.modinnd.nois.cv2)         <- paste0("hetero.g.glm.modinnd.nois.cv2.", names(eval.final.estg.glm.modinnd.nois.cv2))
            print("5.cv2.Nois")
            
            #####################################################################################################################
            ########################## 6. g: GLM Model, inside node, Misspecified adjustment model, Cv2 #########################
            #####################################################################################################################
            data.used.full.bin.mixed.mis <- data.used.full.bin.mixed %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estg.glm.modinnd.mis.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed.mis,
                                                                    est.used          = "G",
                                                                    type.var          = "bin",
                                                                    adj.mod.loc       = "node",
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = NULL,
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)
            
            final.tree.estg.glm.modinnd.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.bin.mixed.mis,
                                                                  tree.list         = seq.created.estg.glm.modinnd.mis.cv2$tree.list,             
                                                                  type.var          = "bin",
                                                                  seed              = a[i],
                                                                  n.cv              = 5,
                                                                  adj.mod.loc       = "node",
                                                                  adj.mthd          = "GLM",
                                                                  adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.mis.cv2[[1]],
                                                                     test.data    = data.bin.mixed$test.data,
                                                                     true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                     noise.var    = data.bin.mixed$noise.var,
                                                                     corr.split   = data.bin.mixed$corr.split,
                                                                     where.split  = data.bin.mixed$where.split,
                                                                     dir.split    = data.bin.mixed$dir.split,
                                                                     split.cate   = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.mis.cv2$corr.frst.splt <- eval.cate.corr.frst.splt(large.tree = seq.created.estg.glm.modinnd.mis.cv2$tree.list[[1]],
                                                                                           corr.split = data.bin.mixed$corr.split,
                                                                                           split.cate = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.mis.cv2                <- unlist(eval.final.estg.glm.modinnd.mis.cv2)
            names(eval.final.estg.glm.modinnd.mis.cv2)         <- paste0("hetero.g.glm.modinnd.mis.cv2.", names(eval.final.estg.glm.modinnd.mis.cv2))
            print("6.cv2.Mis")
            
            performance.hetero.g <- c(eval.final.estg.glm.modinnd.true.cv1,
                                      eval.final.estg.glm.modinnd.nois.cv1,
                                      eval.final.estg.glm.modinnd.mis.cv1,
                                      eval.final.estg.glm.modinnd.true.cv2,
                                      eval.final.estg.glm.modinnd.nois.cv2,
                                      eval.final.estg.glm.modinnd.mis.cv2)
            
            
            
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
            ############################# 7. g: GLM Model, inside node, True adjustment model, cv1 ##############################
            #####################################################################################################################
            t0 <- Sys.time()
            # The difference between true adjustment model and the misspecified model is very dependent on the true adjustment model
            seq.created.estg.glm.modinnd.true.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                     est.used          = "G",
                                                                     type.var          = "bin",
                                                                     adj.mod.loc       = "node",
                                                                     adj.mthd          = "GLM",
                                                                     adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))",
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)
            
            final.tree.estg.glm.modinnd.true.cv1 <- EstG.CvMethod1(data.used         = data.used.bin.mixed,
                                                                   tree.list         = seq.created.estg.glm.modinnd.true.cv1$tree.list,
                                                                   lambda.used       = qchisq(0.95, 1),
                                                                   val.sample        = data.validation.bin.mixed,              
                                                                   type.var          = "bin",
                                                                   adj.mod.loc       = "node",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))")
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.true.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.true.cv1[[1]],
                                                                      test.data    = data.bin.mixed$test.data,
                                                                      true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                      noise.var    = data.bin.mixed$noise.var,
                                                                      corr.split   = data.bin.mixed$corr.split,
                                                                      where.split  = data.bin.mixed$where.split,
                                                                      dir.split    = data.bin.mixed$dir.split,
                                                                      split.cate   = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.true.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.true.cv1                <- unlist(eval.final.estg.glm.modinnd.true.cv1)
            names(eval.final.estg.glm.modinnd.true.cv1)         <- paste0("homo.g.glm.modinnd.true.cv1.", names(eval.final.estg.glm.modinnd.true.cv1))
            print("7.cv1.T")
            
            #####################################################################################################################
            ############################ 8. g: GLM Model, inside node, Noisy adjustment model, cv1 ##############################
            #####################################################################################################################
            t0 <- Sys.time()
            # The difference between true adjustment model and the misspecified model is very dependent on the true adjustment model
            seq.created.estg.glm.modinnd.nois.cv1 <- create.sequence(data.used         = data.used.bin.mixed,
                                                                     est.used          = "G",
                                                                     type.var          = "bin",
                                                                     adj.mod.loc       = "node",
                                                                     adj.mthd          = "GLM",
                                                                     adj.form.true     = NULL,
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)
            
            final.tree.estg.glm.modinnd.nois.cv1 <- EstG.CvMethod1(data.used         = data.used.bin.mixed,
                                                                   tree.list         = seq.created.estg.glm.modinnd.nois.cv1$tree.list,
                                                                   lambda.used       = qchisq(0.95, 1),
                                                                   val.sample        = data.validation.bin.mixed,              
                                                                   type.var          = "bin",
                                                                   adj.mod.loc       = "node",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.nois.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.nois.cv1[[1]],
                                                                      test.data    = data.bin.mixed$test.data,
                                                                      true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                      noise.var    = data.bin.mixed$noise.var,
                                                                      corr.split   = data.bin.mixed$corr.split,
                                                                      where.split  = data.bin.mixed$where.split,
                                                                      dir.split    = data.bin.mixed$dir.split,
                                                                      split.cate   = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.nois.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.nois.cv1                <- unlist(eval.final.estg.glm.modinnd.nois.cv1)
            names(eval.final.estg.glm.modinnd.nois.cv1)         <- paste0("homo.g.glm.modinnd.nois.cv1.", names(eval.final.estg.glm.modinnd.nois.cv1))
            print("8.cv1.Nois")
            
            #####################################################################################################################
            ########################## 9. g: GLM Model, inside node, Misspecified adjustment model, cv1 #########################
            #####################################################################################################################
            data.used.bin.mixed.mis <- data.used.bin.mixed %>%
              select(-X2)
            data.validation.bin.mixed.mis <- data.validation.bin.mixed %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estg.glm.modinnd.mis.cv1 <- create.sequence(data.used         = data.used.bin.mixed.mis,
                                                                    est.used          = "G",
                                                                    type.var          = "bin",
                                                                    adj.mod.loc       = "node",
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = NULL,
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)
            
            final.tree.estg.glm.modinnd.mis.cv1 <- EstG.CvMethod1(data.used         = data.used.bin.mixed.mis,
                                                                  tree.list         = seq.created.estg.glm.modinnd.mis.cv1$tree.list,
                                                                  lambda.used       = qchisq(0.95, 1),
                                                                  val.sample        = data.validation.bin.mixed.mis,              
                                                                  type.var          = "bin",
                                                                  adj.mod.loc       = "node",
                                                                  adj.mthd          = "GLM",
                                                                  adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.mis.cv1 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.mis.cv1[[1]],
                                                                     test.data    = data.bin.mixed$test.data,
                                                                     true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                     noise.var    = data.bin.mixed$noise.var,
                                                                     corr.split   = data.bin.mixed$corr.split,
                                                                     where.split  = data.bin.mixed$where.split,
                                                                     dir.split    = data.bin.mixed$dir.split,
                                                                     split.cate   = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.mis.cv1$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.mis.cv1                <- unlist(eval.final.estg.glm.modinnd.mis.cv1)
            names(eval.final.estg.glm.modinnd.mis.cv1)         <- paste0("homo.g.glm.modinnd.mis.cv1.", names(eval.final.estg.glm.modinnd.mis.cv1))
            print("9.cv1.Mis")
            
            #####################################################################################################################
            ############################# 10. g: GLM Model, inside node, True adjustment model, Cv2 #############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modinnd.true.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                     est.used          = "G",
                                                                     type.var          = "bin",
                                                                     adj.mod.loc       = "node",
                                                                     adj.mthd          = "GLM",
                                                                     adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))",
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)
            
            final.tree.estg.glm.modinnd.true.cv2 <- EstG.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                   tree.list         = seq.created.estg.glm.modinnd.true.cv2$tree.list,             
                                                                   type.var          = "bin",
                                                                   seed              = a[i],
                                                                   n.cv              = 5,
                                                                   adj.mod.loc       = "node",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + X2 + ((X4 == 'B') | (X4 == 'D'))")
            
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.true.cv2[[1]],
                                                                      test.data    = data.bin.mixed$test.data,
                                                                      true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                      noise.var    = data.bin.mixed$noise.var,
                                                                      corr.split   = data.bin.mixed$corr.split,
                                                                      where.split  = data.bin.mixed$where.split,
                                                                      dir.split    = data.bin.mixed$dir.split,
                                                                      split.cate   = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.true.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.true.cv2                <- unlist(eval.final.estg.glm.modinnd.true.cv2)
            names(eval.final.estg.glm.modinnd.true.cv2)         <- paste0("homo.g.glm.modinnd.true.cv2.", names(eval.final.estg.glm.modinnd.true.cv2))
            print("10.cv2.T")
            
            #####################################################################################################################
            ############################# 11. g: GLM Model, inside node, Noisy adjustment model, Cv2 ############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modinnd.nois.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed,
                                                                     est.used          = "G",
                                                                     type.var          = "bin",
                                                                     adj.mod.loc       = "node",
                                                                     adj.mthd          = "GLM",
                                                                     adj.form.true     = NULL,
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)
            
            final.tree.estg.glm.modinnd.nois.cv2 <- EstG.CvMethod2(data.used         = data.used.full.bin.mixed,
                                                                   tree.list         = seq.created.estg.glm.modinnd.nois.cv2$tree.list,             
                                                                   type.var          = "bin",
                                                                   seed              = a[i],
                                                                   n.cv              = 5,
                                                                   adj.mod.loc       = "node",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.nois.cv2[[1]],
                                                                      test.data    = data.bin.mixed$test.data,
                                                                      true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                      noise.var    = data.bin.mixed$noise.var,
                                                                      corr.split   = data.bin.mixed$corr.split,
                                                                      where.split  = data.bin.mixed$where.split,
                                                                      dir.split    = data.bin.mixed$dir.split,
                                                                      split.cate   = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.nois.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.nois.cv2                <- unlist(eval.final.estg.glm.modinnd.nois.cv2)
            names(eval.final.estg.glm.modinnd.nois.cv2)         <- paste0("homo.g.glm.modinnd.nois.cv2.", names(eval.final.estg.glm.modinnd.nois.cv2))
            print("11.cv2.Nois")
            
            #####################################################################################################################
            ########################## 12. g: GLM Model, inside node, Misspecified adjustment model, Cv2 ########################
            #####################################################################################################################
            data.used.full.bin.mixed.mis <- data.used.full.bin.mixed %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estg.glm.modinnd.mis.cv2 <- create.sequence(data.used         = data.used.full.bin.mixed.mis,
                                                                    est.used          = "G",
                                                                    type.var          = "bin",
                                                                    adj.mod.loc       = "node",
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = NULL,
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)
            
            final.tree.estg.glm.modinnd.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.bin.mixed.mis,
                                                                  tree.list         = seq.created.estg.glm.modinnd.mis.cv2$tree.list,             
                                                                  type.var          = "bin",
                                                                  seed              = a[i],
                                                                  n.cv              = 5,
                                                                  adj.mod.loc       = "node",
                                                                  adj.mthd          = "GLM",
                                                                  adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.mis.cv2[[1]],
                                                                     test.data    = data.bin.mixed$test.data,
                                                                     true.trt.eff = data.bin.mixed$true.trt.eff,
                                                                     noise.var    = data.bin.mixed$noise.var,
                                                                     corr.split   = data.bin.mixed$corr.split,
                                                                     where.split  = data.bin.mixed$where.split,
                                                                     dir.split    = data.bin.mixed$dir.split,
                                                                     split.cate   = data.bin.mixed$split.cate)
            eval.final.estg.glm.modinnd.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.mis.cv2                <- unlist(eval.final.estg.glm.modinnd.mis.cv2)
            names(eval.final.estg.glm.modinnd.mis.cv2)         <- paste0("homo.g.glm.modinnd.mis.cv2.", names(eval.final.estg.glm.modinnd.mis.cv2))
            print("12.cv2.Mis")
            
            performance.homo.g <- c(eval.final.estg.glm.modinnd.true.cv1,
                                    eval.final.estg.glm.modinnd.nois.cv1,
                                    eval.final.estg.glm.modinnd.mis.cv1,
                                    eval.final.estg.glm.modinnd.true.cv2,
                                    eval.final.estg.glm.modinnd.nois.cv2,
                                    eval.final.estg.glm.modinnd.mis.cv2)
            
            # print must be put before output, otherwise the output will be print output
            if (i%%30 == 0) {print(i)}
            
            c(performance.hetero.g, performance.homo.g)
            
          }

save(performance.gInnd, file = paste0("../Data/AppendixC7/BinMixedG.RData"))

