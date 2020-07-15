library(causalTree)
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

# set seed
set.seed(a[1]) 

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

ind.train <- sample(1:nrow(rhc.model), round(nrow(rhc.model) * 0.8), replace = F)
rhc.model.train <- rhc.model[ind.train, ]
rhc.model.val   <- rhc.model[!(1:nrow(rhc.model) %in% ind.train), ]

# # Causal Tree Original
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = rhc.model[, !colnames(rhc.model) %in% c("Y")],
                          method    = "GLM",
                          form.true = NULL)
tmp.propsc$prop.sc <- ifelse(rhc.model$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

formula.tree <- as.formula(paste("Y ~ ", paste(colnames(rhc.model)[3:ncol(rhc.model)], collapse= "+")))
orig.ct.propsc <- causalTree(formula.tree,
                             data         = rhc.model,
                             weights      = 1 / tmp.propsc$prop.sc,
                             treatment    = rhc.model$A,
                             split.Rule   = "CT",
                             cv.option    = "CT",
                             split.Honest = T,
                             cv.Honest    = T,
                             split.Bucket = F,
                             xval         = 5,
                             cp           = 0,
                             minsize      = 20)
cptable.propsc <- orig.ct.propsc$cptable[,1][which.min(orig.ct.propsc$cptable[,4])]
orig.final.tree.propsc <- prune(orig.ct.propsc, cptable.propsc)
t1 <- Sys.time()

t.orig.ct <- as.numeric(difftime(t1, t0, units = "secs"))
save(orig.ct.propsc, orig.final.tree.propsc,
     t.orig.ct,
     file = "../Data/main/RhcOptCtOrig.RData")

# Causal Tree Best Main
t0 <- Sys.time()
tmp.propsc <- est.prop.sc(df.noy    = rhc.model[, !colnames(rhc.model) %in% c("Y")],
                          method    = "GLM",
                          form.true = NULL)
tmp.propsc$prop.sc <- ifelse(rhc.model$A == 1, tmp.propsc$prop.sc, 1 - tmp.propsc$prop.sc)

formula.tree <- as.formula(paste("Y ~ ", paste(colnames(rhc.model)[3:ncol(rhc.model)], collapse= "+")))
best.ct.propsc <- causalTree(formula.tree, 
                             data         = rhc.model,
                             weights      = 1 / tmp.propsc$prop.sc,
                             treatment    = rhc.model$A,
                             split.Rule   = "tstats", 
                             split.Honest = T, 
                             cv.option    = "matching", 
                             split.Bucket = F, 
                             xval         = 5, 
                             cp           = 0, 
                             minsize      = 20)
cptable.propsc <- best.ct.propsc$cptable[,1][which.min(best.ct.propsc$cptable[,4])]
best.final.tree.propsc <- prune(best.ct.propsc, cptable.propsc)
t1 <- Sys.time()

t.best.ct <- as.numeric(difftime(t1, t0, units = "secs"))
save(best.ct.propsc, best.final.tree.propsc, 
     t.best.ct, 
     file = "../Data/main/RhcOptCtBest.RData")
