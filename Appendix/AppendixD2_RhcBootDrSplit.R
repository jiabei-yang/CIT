#!/usr/bin/env Rscript
load("../seed1000.rda")

setwd("../")
folder <- paste(getwd(), "/Functions/", sep="")
functions <- list.files(folder)
functions <- grep("main.R", functions, invert = T, value = T)
functions <- paste(folder, functions, sep = "")
for (i in functions){
  source(i)
}

setwd("Appendix/")

# Rscript running part
registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))

# Read in the arguments from command line
option_list = list(
  make_option(c("-s", "--start"), type="integer"),
  make_option(c("-e", "--end"), type = "integer"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
start <- opt$start
end   <- opt$end

# Data Preparation
rhc.orig <- read.csv("../Data/RHC-cleanup/rhc.csv")

# Remove the first column since it contains unnecessary row numbers
rhc.used <- rhc.orig[, -1]

rhc.used <- rhc.used %>%
  mutate(SadmdDte = as.Date(sadmdte,  origin = "1960-01-01"),
         DschDte  = as.Date(dschdte,  origin = "1960-01-01"),
         DthDte   = as.Date(dthdte,   origin = "1960-01-01"),
         LstctDte = as.Date(lstctdte, origin = "1960-01-01"))

# Factor
for (i in 9:20) {
  rhc.used[, i] <- as.factor(rhc.used[, i])
}

rhc.used <- rhc.used %>%
  mutate(adld3p = as.factor(adld3p))

# Replace NA in cat2 since they have practical meaning
rhc.used <- rhc.used %>%
  mutate(imputed.cat2 = as.factor(ifelse(is.na(cat2), "NA", as.character(cat2))))

# Change factor treatment and outcome to numeric treatment
rhc.used <- rhc.used %>%
  mutate(num.swang1 = ifelse(swang1 == "No RHC", 0, 1))
rhc.used <- rhc.used %>%
  mutate(num.death = ifelse(death == "No", 0, 1))
rhc.used <- rhc.used %>%
  mutate(num.dth30 = ifelse(dth30 == "No", 0, 1))

rhc.used <- rhc.used %>%
  mutate(sqrt.wblc1 = sqrt(wblc1),
         log.alb1   = log(alb1),
         log.bili1  = log(bili1),
         log.crea1  = log(crea1),
         sqrt.urin1 = sqrt(urin1))

# Create wt0 variable
num.col <- ncol(rhc.used)
rhc.used <- rhc.used %>%
  mutate(wt0 = ifelse(wtkilo1 == 0, 1, 0))

rhc.model <- data.frame(A = rhc.used$num.swang1,
                        Y = rhc.used$num.dth30,
                        rhc.used[, c("cat1", "imputed.cat2", "ca", 
                                     "cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx", 
                                     "renalhx", "liverhx", "gibledhx", "malighx", "immunhx",
                                     "transhx", "amihx", "age", "sex", "edu", 
                                     "surv2md1", "das2d3pc", "aps1", "scoma1", "meanbp1", 
                                     "wblc1", "hrt1", "resp1", "temp1", "pafi1", 
                                     "alb1", "hema1", "bili1", "crea1", "sod1", 
                                     "pot1", "paco21", "ph1", "wtkilo1", "dnr1",
                                     "ninsclas", "resp", "card", "neuro", "gastr",
                                     "renal", "meta", "hema", "seps", "trauma",
                                     "ortho", "race", "income", "wt0")])

rhc.boot.dr.result <- 
  foreach(i = start:end,
          .combine = "rbind") %dopar% {
            
            set.seed(a[i])
            
            ind.tree  <- sample(1:nrow(rhc.model), ceiling(nrow(rhc.model) * 0.5), replace = F)
            ind.train <- sample(ind.tree, round(length(ind.tree) * 0.8), replace = F)
            ind.val   <- ind.tree[!(ind.tree %in% ind.train)]
            
            rhc.model.train <- rhc.model[ind.train, ]
            rhc.model.val   <- rhc.model[ind.val, ]
            
            t0 <- Sys.time()
            seq.created.orig.cv1 <- create.sequence(data.used         = rhc.model.train,
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
            t1 <- Sys.time()
            
            final.tree.orig.cv1 <- EstDr.CvMethod1(data.used         = rhc.model.train,
                                                   tree.list         = seq.created.orig.cv1$tree.list, 
                                                   lambda.used       = qchisq(0.95, 1),
                                                   val.sample        = rhc.model.val,
                                                   type.var          = "bin",
                                                   propsc.mod.loc    = "node",
                                                   propsc.mthd       = "GLM",
                                                   propsc.form.true  = NULL,
                                                   adj.mod.loc       = "node",
                                                   adj.mthd          = "GLM",
                                                   adj.form.true     = NULL)
            t2 <- Sys.time()
            
            # time
            t.step.12 <- as.numeric(difftime(t1, t0, units = "secs"))
            t.step.3  <- as.numeric(difftime(t2, t1, units = "secs")) 
            
            # size of final tree
            size.tree <- sum(final.tree.orig.cv1[[1]]$frame$var == "<leaf>")
            
            all.nodes <- rownames(seq.created.orig.cv1$tree.list[[1]]$frame)
            
            ind.node.1  <- which(all.nodes == 1)
            frst.splits <- as.character(seq.created.orig.cv1$tree.list[[1]]$frame$var[ind.node.1])
            ind.splits  <- which(seq.created.orig.cv1$tree.list[[1]]$splits[, 3] == seq.created.orig.cv1$tree.list[[1]]$frame$split.stat[ind.node.1])
            dir.split   <- seq.created.orig.cv1$tree.list[[1]]$splits[ind.splits, 2]
            ind.csplit  <- seq.created.orig.cv1$tree.list[[1]]$splits[ind.splits, 4]
            ncol.csplit <- ncol(seq.created.orig.cv1$tree.list[[1]]$csplit)
            
            # record the corresponding splitting point
            # need both 1) index of csplit be an integer and 2) ncat greater than or equal to 2
            if (((ind.csplit %% 1) == 0) & (dir.split >= 2)) {
              # placeholder for continuous split
              frst.splits <- c(frst.splits, rep(NA, 2))
              frst.splits <- c(frst.splits, seq.created.orig.cv1$tree.list[[1]]$csplit[ind.csplit, ])
            } else {
              frst.splits <- c(frst.splits, seq.created.orig.cv1$tree.list[[1]]$splits[ind.splits, c(2, 4)])
              frst.splits <- c(frst.splits, rep(NA, ncol.csplit))
            }
            
            
            for (split.i in 2:7) {
              
              if (split.i %in% all.nodes) {
                
                ind.node  <- which(all.nodes == split.i)
                split.var <- as.character(seq.created.orig.cv1$tree.list[[1]]$frame$var[ind.node])
                
                # if terminal node, put NA as placeholders
                if (split.var == "<leaf>") {
                  frst.splits <- c(frst.splits, rep(NA, 1 + 2 + ncol.csplit))
                  
                } else { # assign split information if not terminal node
                  frst.splits <- c(frst.splits, split.var)
                  ind.splits  <- which(seq.created.orig.cv1$tree.list[[1]]$splits[, 3] == seq.created.orig.cv1$tree.list[[1]]$frame$split.stat[ind.node])
                  dir.split   <- seq.created.orig.cv1$tree.list[[1]]$splits[ind.splits, 2]
                  ind.csplit  <- seq.created.orig.cv1$tree.list[[1]]$splits[ind.splits, 4]
                  
                  # record the corresponding splitting point
                  # need both 1) index of csplit be an integer and 2) ncat greater than or equal to 2
                  if (((ind.csplit %% 1) == 0) & (dir.split >= 2)) {
                    # placeholder for continuous split
                    frst.splits <- c(frst.splits, rep(NA, 2))
                    frst.splits <- c(frst.splits, seq.created.orig.cv1$tree.list[[1]]$csplit[ind.csplit, ])
                  } else {
                    frst.splits <- c(frst.splits, seq.created.orig.cv1$tree.list[[1]]$splits[ind.splits, c(2, 4)])
                    frst.splits <- c(frst.splits, rep(NA, ncol.csplit))
                  }
                  
                }
                
              } else {
                frst.splits <- c(frst.splits, rep(NA, 1 + 2 + ncol.csplit))
              }
              
            }
            
            # print must be put before output, otherwise the output will be print output
            # if (i%%20 == 0) {print(i)}
            print(c(i, length(c(t.step.12, t.step.3, size.tree, frst.splits, ind.tree))))
                        
            c(t.step.12, t.step.3, size.tree, frst.splits, ind.tree)
            
          }

save(rhc.boot.dr.result, file = paste0("../Data/AppendixD2/RhcBootDrSplit.RData"))

