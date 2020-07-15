#!/usr/bin/env Rscript
# load("../Data/Revision/archive/RhcDrCitOpt20200509.RData")
# load("../Data/Revision/RhcWhole/RhcOptDr20200706.RData")
load("../seed1000.rda")

setwd("../")
folder <- paste(getwd(), "/Functions/", sep="")
functions <- list.files(folder)
functions <- grep("main.R", functions, invert = T, value = T)
functions <- paste(folder, functions, sep = "")
for (i in functions){
  source(i)
}

setwd("main/")

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

trt.eff <- NULL
for (i in 1:1000) {
  
  # set seed
  set.seed(a[i]) 
  
  ind.boot <- sample.int(nrow(rhc.model), replace = T)
  rhc.model.boot <- rhc.model[ind.boot, ]
  
  tmp <- gen.fullrank.ipw(df.noy           = rhc.model.boot[, !colnames(rhc.model.boot) %in% c("Y")], 
                          propsc.form.true = NULL)
  tmp.propsc <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                            method    = "GLM",
                            form.true = tmp$propsc.form.true.updated)
  tmp.propsc <- tmp.propsc$prop.sc
  
  tmp <- gen.fullrank.g(df            = rhc.model.boot,
                        adj.form.true = NULL)
  tmp.g <- withWarnings(est.cond.eff(df        = tmp$df.fullrank,
                                     method    = "GLM",
                                     form.true = tmp$adj.form.true.updated, 
                                     type.var  = "bin"))
  
  ipw.trt.eff <- mean(rhc.model.boot$Y * rhc.model.boot$A / tmp.propsc) - 
    mean(rhc.model.boot$Y * (1 - rhc.model.boot$A) / (1 - tmp.propsc))
  
  g.trt.eff <- mean(tmp.g$value$pred.A.1 - tmp.g$value$pred.A.0)
  
  dr.trt.eff <- mean(rhc.model.boot$A * (rhc.model.boot$Y - tmp.g$value$pred.A.1) / tmp.propsc + tmp.g$value$pred.A.1) - 
    mean((1 - rhc.model.boot$A) * (rhc.model.boot$Y - tmp.g$value$pred.A.0) / (1 - tmp.propsc) + tmp.g$value$pred.A.0)
  
  if ((i %% 20) == 0) {print(i)}
  
  trt.eff <- rbind(trt.eff, c(ipw.trt.eff, g.trt.eff, dr.trt.eff))
}

save(trt.eff, file = "../Data/main/RhcRootBi.RData")
# load("../Data/main/RhcRootBi.RData")

# CI for the treatment effect estimate
trt.eff.ci <- data.frame(apply(trt.eff, 2, quantile, probs = c(0.025, 0.5, 0.975)))
colnames(trt.eff.ci) <- c("IPW", "G", "DR")
round(trt.eff.ci, 3)




