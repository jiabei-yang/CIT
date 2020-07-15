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

nrow.rhc.model <- nrow(rhc.model)
boot.R <- 10^3

ci.cit <- 
  foreach(i = start:end,
          .combine  = "rbind") %dopar% {
            
            set.seed(a[i])
            
            ind.tree  <- sample(1:nrow.rhc.model, ceiling(nrow.rhc.model * 0.5), replace = F)
            ind.est   <- (1:nrow.rhc.model)[!((1:nrow.rhc.model) %in% ind.tree)]
            
            rhc.model.est <- rhc.model[ind.est, ]
            
            trt.eff <- NULL
            for (boot.i in 1:boot.R) {
              
              set.seed(a[boot.i])
              
              rhc.model.est.used <- rhc.model.est[sample.int(nrow(rhc.model.est), replace = T), ]
              
              tmp <- gen.fullrank.ipw(df.noy           = rhc.model.est.used[, !colnames(rhc.model.est.used) %in% c("Y")], 
                                      propsc.form.true = NULL)
              tmp.propsc <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                                        method    = "GLM",
                                        form.true = tmp$propsc.form.true.updated)
              tmp.propsc <- tmp.propsc$prop.sc
              
              tmp <- gen.fullrank.g(df            = rhc.model.est.used,
                                    adj.form.true = NULL)
              tmp.g <- withWarnings(est.cond.eff(df        = tmp$df.fullrank,
                                                 method    = "GLM",
                                                 form.true = tmp$adj.form.true.updated, 
                                                 type.var  = "bin"))
              
              ipw.trt.eff <- mean(rhc.model.est.used$Y * rhc.model.est.used$A / tmp.propsc) - 
                mean(rhc.model.est.used$Y * (1 - rhc.model.est.used$A) / (1 - tmp.propsc))
              
              g.trt.eff <- mean(tmp.g$value$pred.A.1 - tmp.g$value$pred.A.0)
              
              dr.trt.eff <- mean(rhc.model.est.used$A * (rhc.model.est.used$Y - tmp.g$value$pred.A.1) / tmp.propsc + tmp.g$value$pred.A.1) - 
                mean((1 - rhc.model.est.used$A) * (rhc.model.est.used$Y - tmp.g$value$pred.A.0) / (1 - tmp.propsc) + tmp.g$value$pred.A.0)
              
              trt.eff <- rbind(trt.eff, c(ipw.trt.eff, g.trt.eff, dr.trt.eff))
              
            } # for boot.i loop
            
            res <- as.vector(apply(trt.eff, 2, quantile, probs = c(0.025, 0.5, 0.975)))
            
            if ((i %% 20) == 0) {print(i)}
            
            res
            
          }

save(ci.cit, file = paste0("../Data/AppendixD3/RhcSpltRootCi.RData"))

 