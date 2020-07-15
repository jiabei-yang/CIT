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

#####################################################################################################################
######################################################## IPW ########################################################
#####################################################################################################################
t0 <- Sys.time()
seq.created.ipwInnd.cv1 <- create.sequence(data.used         = rhc.model.train,
                                           est.used          = "IPW",
                                           type.var          = "bin",
                                           propsc.mod.loc    = "node",
                                           propsc.mthd       = "GLM", 
                                           propsc.form.true  = NULL,
                                           num.truc.obs      = 30,
                                           min.node          = 20)
t1 <- Sys.time()

final.tree.ipwInnd.cv1 <- EstIpw.CvMethod1(data.used         = rhc.model.train, 
                                           tree.list         = seq.created.ipwInnd.cv1$tree.list, 
                                           lambda.used       = qchisq(0.95, 1), 
                                           val.sample        = rhc.model.val, 
                                           type.var          = "bin",
                                           propsc.mod.loc    = "node", 
                                           propsc.mthd       = "GLM", 
                                           propsc.form.true  = NULL)
t2 <- Sys.time()

# time
t.ipwInnd.step.12 <- as.numeric(difftime(t1, t0, units = "secs"))
t.ipwInnd.step.3  <- as.numeric(difftime(t2, t1, units = "secs")) 
print("IPW")

save(seq.created.ipwInnd.cv1, final.tree.ipwInnd.cv1,
     t.ipwInnd.step.12, t.ipwInnd.step.3,
     file = "../Data/main/RhcOptIpw.RData")

