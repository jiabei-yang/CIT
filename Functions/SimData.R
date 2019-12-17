# Continuous outcome, continuous covariates
# Heterogeneous
makeData.cont.eff.cont <- function(N, n.test, p = 6, coeff.prop.sc){
  
  # Covariates
  # p columns of continuous X
  sgm <- diag(p)
  sgm <- sgm + 0.3
  sgm <- sgm - diag(p) * 0.3
  X <- mvrnorm(N, mu = rep(0, p), Sigma = sgm)
  X <- data.frame(X)
  
  # Treatment
  # prop.sc <- 1 / (1 + exp(- 0.6 * X[, 1] + 0.6 * X[, 2] - 0.6 * X[, 3]))
  prop.sc <- 1 / (1 + exp(- coeff.prop.sc * X[, 1] + coeff.prop.sc * X[, 2] - coeff.prop.sc * X[, 3]))
  A <- rbinom(N, 1, prop.sc)
  
  # Outcome
  Y <- 2 + 2 * A + 2 * (X[, 1] < 0) + exp(X[, 2]) + 3 * A * (X[, 4] > 0) + (X[, 5])^3 + rnorm(N)
  # Y <- 2 + 2 * A + A * (X[, 1] < 0.5) * (X[, 4] > 0) + exp(X[, 2]) + 3 * A * (X[, 4] > 0) + (X[, 5])^3 + rnorm(N)
  data.used <- data.frame(A, Y, X)
  colnames(data.used)[3:dim(data.used)[2]] <- paste("X", 1:p, sep = "")
  
  X.test <- mvrnorm(n.test, mu = rep(0, p), Sigma = sgm)
  X.test <- data.frame(X.test)
  test.data <- X.test
  colnames(test.data)[1:dim(test.data)[2]] <- paste("X", 1:p, sep = "")
  
  # True treatment effect
  true.trt.eff <- 2 + 3 * (test.data[, 4] > 0)
  # true.trt.eff <- 2 + (X[, 1] < 0.5) * (X[, 4] > 0) + 3 * (test.data[, 4] > 0)
  return(list(data.used    = data.used, 
              test.data    = test.data,
              true.trt.eff = true.trt.eff,
              noise.var    = c("X1", "X2", "X3", "X5", "X6"),
              corr.split   = c("X4"),              # corr.split needs to be ordered by treatment effect difference
              where.split  = list(c(1)),     # where to make splits (rownumbers in rpart$frame) corresponds to corr.split
              dir.split    = list(c(NULL))))     # direction of making splits, except for the splits into terminal nodes,
                                                                 # other splits' direction will matter; also corresponds to corr.split

  # return(list(data.used    = data.used, 
  #             test.data    = test.data,
  #             true.trt.eff = true.trt.eff,
  #             noise.var    = c("X2", "X3", "X5", "X6"),
  #             corr.split   = c("X4", "X1"),              # corr.split needs to be ordered by treatment effect difference
  #             where.split  = list(c(1, 2), c(1, 3)),     # where to make splits (rownumbers in rpart$frame) corresponds to corr.split
  #             dir.split    = list(c(NULL))))     # direction of making splits, except for the splits into terminal nodes,

  
}

# Continuous outcome, continuous covariates
# Homogeneous
makeData.cont.noeff.cont <- function(N, n.test, p = 6, coeff.prop.sc){
  
  # Covariates
  # p columns of continuous X
  sgm <- diag(p)
  sgm <- sgm + 0.3
  sgm <- sgm - diag(p) * 0.3
  X <- mvrnorm(N, mu = rep(0, p), Sigma = sgm)
  X <- data.frame(X)
  
  # Treatment
  # prop.sc <- 1 / (1 + exp(- 0.6 * X[, 1]))
  prop.sc <- 1 / (1 + exp(- coeff.prop.sc * X[, 1] + coeff.prop.sc * X[, 2] - coeff.prop.sc * X[, 3]))
  # prop.sc <- 0.5
  A <- rbinom(N, 1, prop.sc)
  
  # Outcome
  Y <- 2 + 2 * A + 2 * (X[, 1] < 0) + exp(X[, 2]) + 3 * (X[, 4] > 0) + (X[, 5])^3 + rnorm(N)
  data.used <- data.frame(A, Y, X)
  colnames(data.used)[3:dim(data.used)[2]] <- paste("X", 1:p, sep = "")
  
  X.test <- mvrnorm(n.test, mu = rep(0, p), Sigma = sgm)
  X.test <- data.frame(X.test)
  test.data <- X.test
  colnames(test.data)[1:dim(test.data)[2]] <- paste("X", 1:p, sep = "")
  
  # True treatment effect
  true.trt.eff <- rep(2, n.test)
  return(list(data.used      = data.used, 
              test.data      = test.data,
              true.trt.eff   = true.trt.eff,
              noise.var      = c("X1", "X2", "X3", "X4", "X5", "X6"),
              corr.split     = NULL,                 # corr.split needs to be ordered by treatment effect difference
              where.split    = NULL,
              dir.split      = NULL)) 
  
}

# Continuous outcome, mixed covariates
# Heterogeneous
makeData.cont.eff.mixed = function(N, n.test, p.cont = 3, p.cate = 3, n.cate = 4:6, coeff.prop.sc){
  
  # Covariates
  # Two columns of continuous X
  # Two column of categorical X, one with 4 levels, one with 5 levels
  sgm <- diag(p.cont)
  sgm <- sgm + 0.3
  sgm <- sgm - diag(p.cont) * 0.3
  X <- mvrnorm(N, mu = rep(0, p.cont), Sigma = sgm)
  
  for (i in 1:p.cate) {
    tmp.cate <- base::sample(1:n.cate[i], size = N, replace = T)
    X <- data.frame(X, as.factor(tmp.cate))
    
    colnames(X)[p.cont + i] <- paste("X", p.cont + i, sep = "")
    levels(X[, ncol(X)]) <- LETTERS[1:n.cate[i]]
  }
  
  # Treatment
  prop.sc <- 1 / (1 + exp(- coeff.prop.sc * X[, p.cont - 1] + coeff.prop.sc * X[, p.cont] - coeff.prop.sc * (X[, p.cont + 3] %in% LETTERS[2:3])))
  A <- rbinom(N, 1, prop.sc)
  
  # Outcome
  Y <- 2 + 2 * A + 2 * (X[, 1] < 0) + exp(X[, 2]) + 3 * A * (X[, p.cont + 1] %in% LETTERS[c(2, 4)]) + rnorm(N) 
    # A * (X[, p.cont + 1] %in% LETTERS[c(2, 4)]) * (X[, p.cont + 2] %in% LETTERS[1:2]) + rnorm(N)
  data.used <- data.frame(A, Y, X)
  
  X.test <- mvrnorm(n.test, mu = rep(0, p.cont), Sigma = sgm)
  for (i in 1:p.cate) {
    tmp.cate <- base::sample(1:n.cate[i], size = n.test, replace = T)
    X.test <- data.frame(X.test, as.factor(tmp.cate))
    
    colnames(X.test)[p.cont + i] <- paste("X", p.cont + i, sep = "")
    levels(X.test[, ncol(X.test)]) <- LETTERS[1:n.cate[i]]
  }
  test.data <- X.test
  
  # True treatment effect
  true.trt.eff <- 2 + 3 * (test.data[, p.cont + 1] %in% LETTERS[c(2, 4)]) 
  return(list(data.used    = data.used, 
              test.data    = test.data,
              true.trt.eff = true.trt.eff,
              noise.var    = c("X1", "X2", "X3", "X5", "X6"),
              corr.split   = c("X4"),              # corr.split needs to be ordered by treatment effect difference
              where.split  = list(c(1)),     # where to make splits (rownumbers in rpart$frame) corresponds to corr.split
              dir.split    = list(c(NULL))))             # direction of making splits, except for the splits into terminal nodes,
                                                         # other splits' direction will matter; also corresponds to corr.split
  
}

# Continuous outcome, mixed covariates
# Homogeneous
makeData.cont.noeff.mixed = function(N, n.test, p.cont = 3, p.cate = 3, n.cate = 4:6, coeff.prop.sc){
  
  # Covariates
  # Two columns of continuous X
  # Two column of categorical X, one with 4 levels, one with 5 levels
  sgm <- diag(p.cont)
  sgm <- sgm + 0.3
  sgm <- sgm - diag(p.cont) * 0.3
  X <- mvrnorm(N, mu = rep(0, p.cont), Sigma = sgm)
  
  for (i in 1:p.cate) {
    tmp.cate <- base::sample(1:n.cate[i], size = N, replace = T)
    X <- data.frame(X, as.factor(tmp.cate))
    
    colnames(X)[p.cont + i] <- paste("X", p.cont + i, sep = "")
    levels(X[, ncol(X)]) <- LETTERS[1:n.cate[i]]
  }
  
  # Treatment
  prop.sc <- 1 / (1 + exp(- coeff.prop.sc * X[, p.cont - 1] + coeff.prop.sc * X[, p.cont] - coeff.prop.sc * (X[, p.cont + 3] %in% LETTERS[2:3])))
  A <- rbinom(N, 1, prop.sc)
  
  # Outcome
  Y <- 2 + 2 * A + 2 * (X[, 1] < 0) + exp(X[, 2]) + 3 * (X[, p.cont + 1] %in% LETTERS[c(2, 4)]) + rnorm(N) 
  # A * (X[, p.cont + 1] %in% LETTERS[c(2, 4)]) * (X[, p.cont + 2] %in% LETTERS[1:2]) + rnorm(N)
  data.used <- data.frame(A, Y, X)
  
  X.test <- mvrnorm(n.test, mu = rep(0, p.cont), Sigma = sgm)
  for (i in 1:p.cate) {
    tmp.cate <- base::sample(1:n.cate[i], size = n.test, replace = T)
    X.test <- data.frame(X.test, as.factor(tmp.cate))
    
    colnames(X.test)[p.cont + i] <- paste("X", p.cont + i, sep = "")
    levels(X.test[, ncol(X.test)]) <- LETTERS[1:n.cate[i]]
  }
  test.data <- X.test
  
  # True treatment effect
  true.trt.eff <- rep(2, n.test)
  return(list(data.used    = data.used, 
              test.data    = test.data,
              true.trt.eff = true.trt.eff,
              noise.var    = c("X1", "X2", "X3", "X4", "X5", "X6"),
              corr.split   = NULL,              # corr.split needs to be ordered by treatment effect difference
              where.split  = NULL,              # where to make splits (rownumbers in rpart$frame) corresponds to corr.split
              dir.split    = NULL))             # direction of making splits, except for the splits into terminal nodes,
  # other splits' direction will matter; also corresponds to corr.split
  
}

# Binary outcome, mixed covariates
# Heterogeneous
makeData.bin.eff.mixed = function(N, n.test, p.cont = 3, p.cate = 3, n.cate = 4:6, coeff.prop.sc, seed = NULL){
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Covariates
  # Two columns of continuous X
  # Two column of categorical X, one with 4 levels, one with 5 levels
  sgm <- diag(p.cont)
  sgm <- sgm + 0.3
  sgm <- sgm - diag(p.cont) * 0.3
  X <- mvrnorm(N, mu = rep(0, p.cont), Sigma = sgm)
  
  for (i in 1:p.cate) {
    tmp.cate <- base::sample(1:n.cate[i], size = N, replace = T)
    X <- data.frame(X, as.factor(tmp.cate))
    
    colnames(X)[p.cont + i] <- paste("X", p.cont + i, sep = "")
    levels(X[, ncol(X)]) <- LETTERS[1:n.cate[i]]
  }
  
  # Treatment
  prop.sc <- 1 / (1 + exp(- coeff.prop.sc * X[, p.cont - 1] + coeff.prop.sc * X[, p.cont] - coeff.prop.sc * (X[, p.cont + 3] %in% LETTERS[2:3])))
  A <- rbinom(N, 1, prop.sc)
  
  # Outcome
  prob.y <- 0.1 + 0.1 * A + 1 / (1 + exp(- 0.2 * X[, 2])) - 0.4 * A * (X[, p.cont + 1] %in% LETTERS[c(2, 4)]) 
  Y      <- rbinom(N, 1, prob = prob.y)
  # A * (X[, p.cont + 1] %in% LETTERS[c(2, 4)]) * (X[, p.cont + 2] %in% LETTERS[1:2]) + rnorm(N)
  data.used <- data.frame(A, Y, X)
  
  X.test <- mvrnorm(n.test, mu = rep(0, p.cont), Sigma = sgm)
  for (i in 1:p.cate) {
    tmp.cate <- base::sample(1:n.cate[i], size = n.test, replace = T)
    X.test <- data.frame(X.test, as.factor(tmp.cate))
    
    colnames(X.test)[p.cont + i] <- paste("X", p.cont + i, sep = "")
    levels(X.test[, ncol(X.test)]) <- LETTERS[1:n.cate[i]]
  }
  test.data <- X.test
  
  # True treatment effect
  true.trt.eff <- 0.1 - 0.4 * (test.data[, p.cont + 1] %in% LETTERS[c(2, 4)])  
  return(list(prob.y       = prob.y,
              data.used    = data.used, 
              test.data    = test.data,
              true.trt.eff = true.trt.eff,
              noise.var    = c("X1", "X2", "X3", "X5", "X6"),
              corr.split   = c("X4"),              # corr.split needs to be ordered by treatment effect difference
              where.split  = list(c(1)),           # where to make splits (rownumbers in rpart$frame), each vector's length is the same as corr.split, each vector gives one possibility of row number of making splits in rpart$frame
              dir.split    = list(c(NULL)),        # direction of making splits, except for the splits into terminal nodes, other splits' direction will matter; also corresponds to corr.split
              split.cate   = list(c(2, 4))))       # Categories that should be in the same node, each vector corresponds to one element in corr.split            
                                                         
  
}

# Binary outcome, mixed covariates
# Homogeneous
makeData.bin.noeff.mixed = function(N, n.test, p.cont = 3, p.cate = 3, n.cate = 4:6, coeff.prop.sc, seed = NULL){
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Covariates
  # Two columns of continuous X
  # Two column of categorical X, one with 4 levels, one with 5 levels
  sgm <- diag(p.cont)
  sgm <- sgm + 0.3
  sgm <- sgm - diag(p.cont) * 0.3
  X <- mvrnorm(N, mu = rep(0, p.cont), Sigma = sgm)
  
  for (i in 1:p.cate) {
    tmp.cate <- base::sample(1:n.cate[i], size = N, replace = T)
    X <- data.frame(X, as.factor(tmp.cate))
    
    colnames(X)[p.cont + i] <- paste("X", p.cont + i, sep = "")
    levels(X[, ncol(X)]) <- LETTERS[1:n.cate[i]]
  }
  
  # Treatment
  prop.sc <- 1 / (1 + exp(- coeff.prop.sc * X[, p.cont - 1] + coeff.prop.sc * X[, p.cont] - coeff.prop.sc * (X[, p.cont + 3] %in% LETTERS[2:3])))
  A <- rbinom(N, 1, prop.sc)
  
  # Outcome
  prob.y <- 0.15 + 0.1 * A + 1 / (1 + exp(- 0.2 * X[, 2])) - 0.4 * (X[, p.cont + 1] %in% LETTERS[c(2, 4)])
  # prob.y[prob.y < 0] <- 0
  # prob.y[prob.y > 1] <- 1
  Y      <- rbinom(N, 1, prob = prob.y)
  # A * (X[, p.cont + 1] %in% LETTERS[c(2, 4)]) * (X[, p.cont + 2] %in% LETTERS[1:2]) + rnorm(N)
  data.used <- data.frame(A, Y, X)
  
  X.test <- mvrnorm(n.test, mu = rep(0, p.cont), Sigma = sgm)
  for (i in 1:p.cate) {
    tmp.cate <- base::sample(1:n.cate[i], size = n.test, replace = T)
    X.test <- data.frame(X.test, as.factor(tmp.cate))
    
    colnames(X.test)[p.cont + i] <- paste("X", p.cont + i, sep = "")
    levels(X.test[, ncol(X.test)]) <- LETTERS[1:n.cate[i]]
  }
  test.data <- X.test
  
  # True treatment effect
  true.trt.eff <- rep(0.1, n.test)
  return(list(prob.y       = prob.y,
              data.used    = data.used, 
              test.data    = test.data,
              true.trt.eff = true.trt.eff,
              noise.var    = c("X1", "X2", "X3", "X4", "X5", "X6"),
              corr.split   = NULL,              # corr.split needs to be ordered by treatment effect difference
              where.split  = NULL,           # where to make splits (rownumbers in rpart$frame), each vector's length is the same as corr.split, each vector gives one possibility of row number of making splits in rpart$frame
              dir.split    = NULL,        # direction of making splits, except for the splits into terminal nodes, other splits' direction will matter; also corresponds to corr.split
              split.cate   = NULL))       # Categories that should be in the same node, each vector corresponds to one element in corr.split            
  
  
}
