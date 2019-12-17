job.number <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
load("seed1000.rda")
set.seed(a[job.number])

setwd("..")
folder <- paste(getwd(), "/Functions/", sep="")
functions <- list.files(folder)
functions <- paste(folder, functions, sep = "")
for (i in functions){
  source(i)
}

setwd("main/")

# Causal Tree
# library(devtools) 
# install_github("susanathey/causalTree")
library(causalTree)

coeff.prop.sc <- 0.6
N.training    <- 10^3
N.testing     <- 10^3

#####################################################################################################################
################################################# Heterogeneous #####################################################
#####################################################################################################################
data.cont.cont            <- makeData.cont.eff.cont(N             = N.training, 
                                                    n.test        = N.testing, 
                                                    coeff.prop.sc = coeff.prop.sc)
data.used.full.cont.cont  <- data.cont.cont$data.used
data.used.cont.cont       <- data.used.full.cont.cont[1:800, ]
# val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
data.validation.cont.cont <- data.used.full.cont.cont[801:1000, ]  


