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
            
            #####################################################################################################################
            ################################################# Heterogeneous #####################################################
            #####################################################################################################################
            data.cont.cont            <- makeData.cont.eff.cont(1000, 1000, coeff.prop.sc = 0.6)
            data.used.full.cont.cont  <- data.cont.cont$data.used
            data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
            # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
            data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ]  
            
            #####################################################################################################################
            ############################# 1. g: GLM Model, inside node, True adjustment model, Cv2 ##############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modinnd.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                     est.used          = "G",
                                                                     type.var          = "cont",
                                                                     adj.mod.loc       = "node",
                                                                     adj.mthd          = "GLM",
                                                                     adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)",
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)
            
            final.tree.estg.glm.modinnd.true.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                   tree.list         = seq.created.estg.glm.modinnd.true.cv2$tree.list,             
                                                                   type.var          = "cont",
                                                                   seed              = a[i],
                                                                   n.cv              = 5,
                                                                   adj.mod.loc       = "node",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X5^3) + A:I(X4 > 0)")
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
            eval.final.estg.glm.modinnd.true.cv2                <- unlist(eval.final.estg.glm.modinnd.true.cv2)
            names(eval.final.estg.glm.modinnd.true.cv2)         <- paste0("hetero.g.glm.modinnd.true.cv2.", names(eval.final.estg.glm.modinnd.true.cv2))
            print("hetero.T")
            
            #####################################################################################################################
            ############################# 2. g: GLM Model, inside node, Noisy adjustment model, Cv2 #############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modinnd.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                     est.used          = "G",
                                                                     type.var          = "cont",
                                                                     adj.mod.loc       = "node",
                                                                     adj.mthd          = "GLM",
                                                                     adj.form.true     = NULL,
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)
            
            final.tree.estg.glm.modinnd.nois.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                   tree.list         = seq.created.estg.glm.modinnd.nois.cv2$tree.list,             
                                                                   type.var          = "cont",
                                                                   seed              = a[i],
                                                                   n.cv              = 5,
                                                                   adj.mod.loc       = "node",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = NULL)
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
            eval.final.estg.glm.modinnd.nois.cv2                <- unlist(eval.final.estg.glm.modinnd.nois.cv2)
            names(eval.final.estg.glm.modinnd.nois.cv2)         <- paste0("hetero.g.glm.modinnd.nois.cv2.", names(eval.final.estg.glm.modinnd.nois.cv2))
            print("hetero.Nois")
            
            #####################################################################################################################
            ########################## 3. g: GLM Model, inside node, Misspecified adjustment model, Cv2 #########################
            #####################################################################################################################
            data.used.full.cont.cont.mis <- data.used.full.cont.cont %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estg.glm.modinnd.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                                    est.used          = "G",
                                                                    type.var          = "cont",
                                                                    adj.mod.loc       = "node",
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = NULL,
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)
            
            final.tree.estg.glm.modinnd.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                                  tree.list         = seq.created.estg.glm.modinnd.mis.cv2$tree.list,             
                                                                  type.var          = "cont",
                                                                  seed              = a[i],
                                                                  n.cv              = 5,
                                                                  adj.mod.loc       = "node",
                                                                  adj.mthd          = "GLM",
                                                                  adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.mis.cv2[[1]],
                                                                     test.data    = data.cont.cont$test.data,
                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                     noise.var    = data.cont.cont$noise.var,
                                                                     corr.split   = data.cont.cont$corr.split,
                                                                     where.split  = data.cont.cont$where.split,
                                                                     dir.split    = data.cont.cont$dir.split)
            
            eval.final.estg.glm.modinnd.mis.cv2$t <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.mis.cv2$corr.frst.splt <- as.character(seq.created.estg.glm.modinnd.mis.cv2$tree.list[[1]]$frame$var[1]) == data.cont.cont$corr.split
            eval.final.estg.glm.modinnd.mis.cv2                <- unlist(eval.final.estg.glm.modinnd.mis.cv2)
            names(eval.final.estg.glm.modinnd.mis.cv2)         <- paste0("hetero.g.glm.modinnd.mis.cv2.", names(eval.final.estg.glm.modinnd.mis.cv2))
            print("hetero.Mis")
            
            performance.hetero.gInnd <- c(eval.final.estg.glm.modinnd.true.cv2, 
                                          eval.final.estg.glm.modinnd.nois.cv2,
                                          eval.final.estg.glm.modinnd.mis.cv2)
            
            #####################################################################################################################
            ################################################### Homogeneous #####################################################
            #####################################################################################################################
            data.cont.cont            <- makeData.cont.noeff.cont(1000, 1000, coeff.prop.sc = 0.6)
            data.used.full.cont.cont  <- data.cont.cont$data.used
            data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
            # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
            data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ]  
            
            #####################################################################################################################
            ############################# 4. g: GLM Model, inside node, True adjustment model, Cv2 ##############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modinnd.true.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                     est.used          = "G",
                                                                     type.var          = "cont",
                                                                     adj.mod.loc       = "node",
                                                                     adj.mthd          = "GLM",
                                                                     adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)",
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)
            
            final.tree.estg.glm.modinnd.true.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                   tree.list         = seq.created.estg.glm.modinnd.true.cv2$tree.list,             
                                                                   type.var          = "cont",
                                                                   seed              = a[i],
                                                                   n.cv              = 5,
                                                                   adj.mod.loc       = "node",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = "Y ~ A + I(X1 < 0) + exp(X2) + I(X4 > 0) + I(X5^3)")
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.true.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.true.cv2[[1]],
                                                                      test.data    = data.cont.cont$test.data,
                                                                      true.trt.eff = data.cont.cont$true.trt.eff,
                                                                      noise.var    = data.cont.cont$noise.var,
                                                                      corr.split   = data.cont.cont$corr.split,
                                                                      where.split  = data.cont.cont$where.split,
                                                                      dir.split    = data.cont.cont$dir.split)
            
            eval.final.estg.glm.modinnd.true.cv2$t      <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.true.cv2        <- unlist(eval.final.estg.glm.modinnd.true.cv2)
            names(eval.final.estg.glm.modinnd.true.cv2) <- paste0("homo.g.glm.modinnd.true.cv2.", names(eval.final.estg.glm.modinnd.true.cv2))
            print("homo.T")
            
            #####################################################################################################################
            ############################# 5. g: GLM Model, inside node, Noisy adjustment model, Cv2 #############################
            #####################################################################################################################
            t0 <- Sys.time()
            seq.created.estg.glm.modinnd.nois.cv2 <- create.sequence(data.used         = data.used.full.cont.cont,
                                                                     est.used          = "G",
                                                                     type.var          = "cont",
                                                                     adj.mod.loc       = "node",
                                                                     adj.mthd          = "GLM",
                                                                     adj.form.true     = NULL,
                                                                     num.truc.obs      = 30,
                                                                     min.node          = 20)
            
            final.tree.estg.glm.modinnd.nois.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont,
                                                                   tree.list         = seq.created.estg.glm.modinnd.nois.cv2$tree.list,             
                                                                   type.var          = "cont",
                                                                   seed              = a[i],
                                                                   n.cv              = 5,
                                                                   adj.mod.loc       = "node",
                                                                   adj.mthd          = "GLM",
                                                                   adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.nois.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.nois.cv2[[1]],
                                                                      test.data    = data.cont.cont$test.data,
                                                                      true.trt.eff = data.cont.cont$true.trt.eff,
                                                                      noise.var    = data.cont.cont$noise.var,
                                                                      corr.split   = data.cont.cont$corr.split,
                                                                      where.split  = data.cont.cont$where.split,
                                                                      dir.split    = data.cont.cont$dir.split)
            
            eval.final.estg.glm.modinnd.nois.cv2$t      <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.nois.cv2        <- unlist(eval.final.estg.glm.modinnd.nois.cv2)
            names(eval.final.estg.glm.modinnd.nois.cv2) <- paste0("homo.g.glm.modinnd.nois.cv2.", names(eval.final.estg.glm.modinnd.nois.cv2))
            print("homo.Nois")
            
            #####################################################################################################################
            ########################## 6. g: GLM Model, inside node, Misspecified adjustment model, Cv2 #########################
            #####################################################################################################################
            data.used.full.cont.cont.mis <- data.used.full.cont.cont %>%
              select(-X2)
            
            t0 <- Sys.time()
            seq.created.estg.glm.modinnd.mis.cv2 <- create.sequence(data.used         = data.used.full.cont.cont.mis,
                                                                    est.used          = "G",
                                                                    type.var          = "cont",
                                                                    adj.mod.loc       = "node",
                                                                    adj.mthd          = "GLM",
                                                                    adj.form.true     = NULL,
                                                                    num.truc.obs      = 30,
                                                                    min.node          = 20)
            
            final.tree.estg.glm.modinnd.mis.cv2 <- EstG.CvMethod2(data.used         = data.used.full.cont.cont.mis,
                                                                  tree.list         = seq.created.estg.glm.modinnd.mis.cv2$tree.list,             
                                                                  type.var          = "cont",
                                                                  seed              = a[i],
                                                                  n.cv              = 5,
                                                                  adj.mod.loc       = "node",
                                                                  adj.mthd          = "GLM",
                                                                  adj.form.true     = NULL)
            t1 <- Sys.time()
            
            eval.final.estg.glm.modinnd.mis.cv2 <- eval.measures.eff(final.tree   = final.tree.estg.glm.modinnd.mis.cv2[[1]],
                                                                     test.data    = data.cont.cont$test.data,
                                                                     true.trt.eff = data.cont.cont$true.trt.eff,
                                                                     noise.var    = data.cont.cont$noise.var,
                                                                     corr.split   = data.cont.cont$corr.split,
                                                                     where.split  = data.cont.cont$where.split,
                                                                     dir.split    = data.cont.cont$dir.split)
            
            eval.final.estg.glm.modinnd.mis.cv2$t      <- as.numeric(difftime(t1, t0, units = "secs"))
            eval.final.estg.glm.modinnd.mis.cv2        <- unlist(eval.final.estg.glm.modinnd.mis.cv2)
            names(eval.final.estg.glm.modinnd.mis.cv2) <- paste0("homo.g.glm.modinnd.mis.cv2.", names(eval.final.estg.glm.modinnd.mis.cv2))
            print("homo.Mis")
            
            performance.homo.gInnd <- c(eval.final.estg.glm.modinnd.true.cv2, 
                                        eval.final.estg.glm.modinnd.nois.cv2,
                                        eval.final.estg.glm.modinnd.mis.cv2)
            
            # print must be put before output, otherwise the output will be print output
            if (i%%30 == 0) {print(i)}
            
            c(performance.hetero.gInnd, performance.homo.gInnd)
            
          }

save(performance.gInnd, file = "../Data/AppendixC4/Cv2G.RData")

