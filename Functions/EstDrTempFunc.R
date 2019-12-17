calc.dr.trt.eff <- function(data.node, x, ux,
                            prop.sc, propsc.mthd, propsc.form.true,                       # propsc.mthd, propsc.form.true are placeholders for insplit models
                            est.cond.eff.0, est.cond.eff.1, adj.mthd, adj.form.true, type.var = "cont"){     # adj.mthd, adj.form.true are placeholders for insplit models
  
  # calculates ipw mean response difference between treatment and control for each level of x 
  # when x is categorical and the outcome is continuous.
  # data.node: a data frame containing observations at the node.
  # x:         vector of the \texttt{x} values at the node.
  # ux:        unique levels of the categorical covariate x at the node.
  
  trt.eff <- NULL
  
  # Crude mean response difference between treatment and control for each level of x
  for (i in 1:length(ux)){
    
    data.sub.node      <- data.node[x == ux[i], ]
    est.cond.eff.0.sub <- est.cond.eff.0[x == ux[i]]
    est.cond.eff.1.sub <- est.cond.eff.1[x == ux[i]]
    prop.sc.sub        <- prop.sc[x == ux[i]]
    
    mu.1 <- mean(data.sub.node$A * (data.sub.node$Y - est.cond.eff.1.sub) / prop.sc.sub + est.cond.eff.1.sub)
    mu.0 <- mean((1 - data.sub.node$A) * (data.sub.node$Y - est.cond.eff.0.sub) / (1 - prop.sc.sub) + est.cond.eff.0.sub)
    
    trt.eff     <- c(trt.eff, mu.1 - mu.0)
    
  } # for loop
  
  return(trt.eff)
  
}

# IN SPLIT NOT WORKING; INNODE INSPLIT EVAL FUNCTION
stemp.dr <- function(y, wt, x, parms, continuous){

  n <- length(y)
  
  sub.ind      <- y
  sub.x        <- parms$covariates[sub.ind, ]
  A            <- parms$trt[sub.ind]
  Y            <- parms$response[sub.ind]
  data.node    <- data.frame(Y, A, sub.x)
  num.truc.obs <- parms$num.truc.obs
  
  prop.sc           <- parms$prop.sc[sub.ind]
  propsc.mthd       <- parms$propsc.mthd
  propsc.form.true  <- parms$propsc.form.true
  propsc.mod.insplt <- parms$propsc.mod.insplt
  
  est.cond.eff.0 <- parms$est.cond.eff.0[sub.ind]
  est.cond.eff.1 <- parms$est.cond.eff.1[sub.ind]
  adj.mod.insplt <- parms$adj.mod.insplt
  adj.mthd       <- parms$adj.mthd
  adj.form.true  <- parms$adj.form.true
  
  type.var <- parms$type.var
  
  if (is.null(prop.sc)) {

    if (!propsc.mod.insplt) {
      
      data.node.noy <- data.node[, !colnames(data.node) %in% c("Y")]
      tmp <- gen.fullrank.ipw(df.noy           = data.node.noy, 
                              propsc.form.true = propsc.form.true)
      tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                         method    = propsc.mthd,
                         form.true = tmp$propsc.form.true.updated)
      prop.sc <- tmp$prop.sc
      
    }
  }
  
  if (is.null(est.cond.eff.0)) {
    
    if (!adj.mod.insplt) {
      
      tmp <- gen.fullrank.g(df            = data.node,
                            adj.form.true = adj.form.true)
      tmp <- est.cond.eff(df        = tmp$df.fullrank,
                          method    = adj.mthd,
                          form.true = tmp$adj.form.true.updated,
                          type.var  = type.var)
  
      est.cond.eff.1 <- tmp$pred.A.1
      est.cond.eff.0 <- tmp$pred.A.0
    }
    
  }
  
  if (continuous){
    
    # Skip the first 10 and last 10 splits
    # goodness <- NULL
    # direction <- NULL
    goodness  <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    if (is.null(est.cond.eff.1)) { 
      
      if (is.null(prop.sc)) { # Both models in split
        
        # Ensure the model is at least fitted on num.truc.obs obs
        for (i in (num.truc.obs:(n-num.truc.obs))) { 
          
          data.node.l <- data.node[1:i, ]
          data.node.r <- data.node[(i+1):n, ]
          
          # Propensity score model
          data.node.noy.l <- data.node.l[, !colnames(data.node.l) %in% c("Y")]
          tmp <- gen.fullrank.ipw(df.noy           = data.node.noy.l, 
                                  propsc.form.true = propsc.form.true)
          tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                             method    = propsc.mthd,
                             form.true = tmp$propsc.form.true.updated)
          prop.sc.l <- tmp$prop.sc
          
          data.node.noy.r <- data.node.r[, !colnames(data.node.r) %in% c("Y")]
          tmp <- gen.fullrank.ipw(df.noy           = data.node.noy.r, 
                                  propsc.form.true = propsc.form.true)
          tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                             method    = propsc.mthd,
                             form.true = tmp$propsc.form.true.updated)
          prop.sc.r <- tmp$prop.sc
          
          prop.sc <- c(prop.sc.l, prop.sc.r)
          
          # Outcome model
          tmp <- gen.fullrank.g(df            = data.node.l,
                                adj.form.true = adj.form.true)
          tmp <- est.cond.eff(df        = tmp$df.fullrank,
                              method    = adj.mthd,
                              form.true = tmp$adj.form.true.updated,
                              type.var  = type.var)
          
          est.cond.eff.1.l <- tmp$pred.A.1
          est.cond.eff.0.l <- tmp$pred.A.0
          
          tmp <- gen.fullrank.g(df            = data.node.r,
                                adj.form.true = adj.form.true)
          tmp <- est.cond.eff(df        = tmp$df.fullrank,
                              method    = adj.mthd,
                              form.true = tmp$adj.form.true.updated,
                              type.var  = type.var)
          
          est.cond.eff.1.r <- tmp$pred.A.1
          est.cond.eff.0.r <- tmp$pred.A.0
          
          est.cond.eff.1 <- c(est.cond.eff.1.l, est.cond.eff.1.r)
          est.cond.eff.0 <- c(est.cond.eff.0.l, est.cond.eff.0.r)
          
          # Calculating the estimator
          mu.1l <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[1:i])
          mu.0l <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[1:i])
          mu.1r <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[(i+1):n])
          mu.0r <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[(i+1):n])
          
          t.dr <- (mu.1l - mu.0l) - (mu.1r - mu.0r)
          p.l <- i/n
          p.r <- (n-i)/n
          I_i <- (x %in% x[1:i]) / p.l * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
                                            ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
            (x %in% x[(i+1):n]) / p.r * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
                                           ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
            t.dr
          
          var <- (mean(I_i^2) - (p.r * (mu.1l - mu.0l) + p.l * (mu.1r - mu.0r))^2 / (p.l * p.r)) / n
          
          # Skip this loop if there is negative estimated variance
          if (var < 0) {
            next
          }
          
          goodness[i] <- (t.dr / sqrt(var))^2
          direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))      
        }
        
      } else { # prop.sc not in split, outcome model in split
        
        # Ensure the model is at least fitted on num.truc.obs obs
        for (i in (num.truc.obs:(n-num.truc.obs))) { 
          
          data.node.l <- data.node[1:i, ]
          data.node.r <- data.node[(i+1):n, ]
          
          # Outcome model
          tmp <- gen.fullrank.g(df            = data.node.l,
                                adj.form.true = adj.form.true)
          tmp <- est.cond.eff(df        = tmp$df.fullrank,
                              method    = adj.mthd,
                              form.true = tmp$adj.form.true.updated,
                              type.var  = type.var)
          
          est.cond.eff.1.l <- tmp$pred.A.1
          est.cond.eff.0.l <- tmp$pred.A.0
          
          tmp <- gen.fullrank.g(df            = data.node.r,
                                adj.form.true = adj.form.true)
          tmp <- est.cond.eff(df        = tmp$df.fullrank,
                              method    = adj.mthd,
                              form.true = tmp$adj.form.true.updated,
                              type.var  = type.var)
          
          est.cond.eff.1.r <- tmp$pred.A.1
          est.cond.eff.0.r <- tmp$pred.A.0
          
          est.cond.eff.1 <- c(est.cond.eff.1.l, est.cond.eff.1.r)
          est.cond.eff.0 <- c(est.cond.eff.0.l, est.cond.eff.0.r)
          
          # Calculating the estimator
          mu.1l <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[1:i])
          mu.0l <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[1:i])
          mu.1r <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[(i+1):n])
          mu.0r <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[(i+1):n])
          
          t.dr <- (mu.1l - mu.0l) - (mu.1r - mu.0r)
          p.l <- i/n
          p.r <- (n-i)/n
          I_i <- (x %in% x[1:i]) / p.l * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
                                            ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
            (x %in% x[(i+1):n]) / p.r * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
                                           ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
            t.dr
          
          var <- (mean(I_i^2) - (p.r * (mu.1l - mu.0l) + p.l * (mu.1r - mu.0r))^2 / (p.l * p.r)) / n
          
          # Skip this loop if there is negative estimated variance
          if (var < 0) {
            next
          }
          
          goodness[i] <- (t.dr / sqrt(var))^2
          direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))      
        }
        
        
      }
      
      
    } else { # Outcome model in node or outside
      
      if (is.null(prop.sc)) { # propensity score model in split
        
        # Ensure the model is at least fitted on num.truc.obs obs
        for (i in (num.truc.obs:(n-num.truc.obs))) { 
          
          data.node.l <- data.node[1:i, ]
          data.node.r <- data.node[(i+1):n, ]
          
          # Propensity score model
          data.node.noy.l <- data.node.l[, !colnames(data.node.l) %in% c("Y")]
          tmp <- gen.fullrank.ipw(df.noy           = data.node.noy.l, 
                                  propsc.form.true = propsc.form.true)
          tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                             method    = propsc.mthd,
                             form.true = tmp$propsc.form.true.updated)
          prop.sc.l <- tmp$prop.sc
          
          data.node.noy.r <- data.node.r[, !colnames(data.node.r) %in% c("Y")]
          tmp <- gen.fullrank.ipw(df.noy           = data.node.noy.r, 
                                  propsc.form.true = propsc.form.true)
          tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                             method    = propsc.mthd,
                             form.true = tmp$propsc.form.true.updated)
          prop.sc.r <- tmp$prop.sc
          
          prop.sc <- c(prop.sc.l, prop.sc.r)
          
          # Calculating the estimator
          mu.1l <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[1:i])
          mu.0l <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[1:i])
          mu.1r <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[(i+1):n])
          mu.0r <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[(i+1):n])
          
          t.dr <- (mu.1l - mu.0l) - (mu.1r - mu.0r)
          p.l <- i/n
          p.r <- (n-i)/n
          I_i <- (x %in% x[1:i]) / p.l * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
                                            ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
            (x %in% x[(i+1):n]) / p.r * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
                                           ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
            t.dr
          
          var <- (mean(I_i^2) - (p.r * (mu.1l - mu.0l) + p.l * (mu.1r - mu.0r))^2 / (p.l * p.r)) / n
          
          # Skip this loop if there is negative estimated variance
          if (var < 0) {
            next
          }
          
          goodness[i] <- (t.dr / sqrt(var))^2
          direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))      
        }
        
      } else { # Both models outside split
        
        # Ensure the model is at least fitted on num.truc.obs obs
        for (i in (num.truc.obs:(n-num.truc.obs))) { 
          
          # Calculating the estimator
          mu.1l <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[1:i])
          mu.0l <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[1:i])
          mu.1r <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[(i+1):n])
          mu.0r <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[(i+1):n])
          
          t.dr <- (mu.1l - mu.0l) - (mu.1r - mu.0r)
          p.l <- i/n
          p.r <- (n-i)/n
          I_i <- (x %in% x[1:i]) / p.l * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
                                            ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
            (x %in% x[(i+1):n]) / p.r * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
                                           ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
            t.dr
          
          var <- (mean(I_i^2) - (p.r * (mu.1l - mu.0l) + p.l * (mu.1r - mu.0r))^2 / (p.l * p.r)) / n
          
          # Skip this loop if there is negative estimated variance
          if (var < 0) {
            next
          }
          
          goodness[i] <- (t.dr / sqrt(var))^2
          direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))      
        }
        
      } # Both models outside split

    } # if (is.null(est.cond.eff.1)) {} else{}

    
    
  } else {
    
    ux <- sort(unique(x))
    
    # Order the levels of the categorical covariates by their treatment effect
    # Currently code only works for outside node and innode model fitting, not for in split
    trt.eff.dr <- calc.dr.trt.eff(data.node, x, ux,
                                 prop.sc          = prop.sc,
                                 propsc.mthd      = NULL,
                                 propsc.form.true = NULL,
                                 est.cond.eff.0   = est.cond.eff.0,
                                 est.cond.eff.1   = est.cond.eff.1,
                                 adj.mthd         = NULL,
                                 adj.form.true    = NULL,
                                 type.var         = NULL)
    
    
    ord.ux <- order(trt.eff.dr)
    goodness  <- rep(0, length(ux) - 1)
    
    # Splits
    for (i in 1:(length(ux) - 1)){
      
      ind.l <- x %in% ux[ord.ux[1:i]]
      ind.r <- x %in% ux[ord.ux[(i+1):length(ux)]] 
        
      mu.1l <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[ind.l])
      mu.0l <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[ind.l])
      mu.1r <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[ind.r])
      mu.0r <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[ind.r])
      
      t.dr <- (mu.1l - mu.0l) - (mu.1r - mu.0r)
      p.l <- sum(ind.l) / n
      p.r <- sum(ind.r) / n
      I_i <- ind.l / p.l * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
                              ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
        ind.r / p.r * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
                         ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
        t.dr
      
      var <- (mean(I_i^2) - (p.r * (mu.1l - mu.0l) + p.l * (mu.1r - mu.0r))^2 / (p.l * p.r)) / n
      
      # Skip this loop if there is negative estimated variance
      if (var < 0) {
        next
      }
      
      goodness[i] <- (t.dr / sqrt(var))^2
      
    }
    
    direction <- ux[ord.ux] 
    
  }
  
  list(goodness  = goodness,
       direction = direction)
}



etemp.dr <- function(y, wt, parms) {
  
  n <- length(y)
  
  sub.ind <- y
  sub.x   <- parms$covariates[sub.ind, ]
  A       <- parms$trt[sub.ind]
  Y       <- parms$response[sub.ind]
  data.node      <- data.frame(Y, A, sub.x)
  
  prop.sc           <- parms$prop.sc[sub.ind]
  propsc.mthd       <- parms$propsc.mthd
  propsc.form.true  <- parms$propsc.form.true
  
  est.cond.eff.0 <- parms$est.cond.eff.0[sub.ind]
  est.cond.eff.1 <- parms$est.cond.eff.1[sub.ind]
  adj.mthd       <- parms$adj.mthd
  adj.form.true  <- parms$adj.form.true
  
  type.var <- parms$type.var
  
  # Need to move this part of code here for return value earlier
  # TO DO: Weird how this does not work when rss is used.
  # not crucial to the performance 
  wmean <- sum(Y*wt)/sum(wt)
  rss <- sum(wt*(Y-wmean)^2)
  
  if (is.null(prop.sc)) {
    
    # If there is only one observation in this node, avg.trt.eff cannot be calculated so return earlier
    if (n == 1) {
      return(list(label = NA, deviance = rss))
    }
    
    data.node.noy <- data.node[, !colnames(data.node) %in% c("Y")]
    tmp <- gen.fullrank.ipw(df.noy           = data.node.noy, 
                            propsc.form.true = propsc.form.true)
    tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                       method    = propsc.mthd,
                       form.true = tmp$propsc.form.true.updated)
    prop.sc <- tmp$prop.sc
    
  }

  if (is.null(est.cond.eff.0)) {
    
    # If there is only one observation in this node, avg.trt.eff cannot be calculated so return earlier
    if (n == 1) {
      return(list(label = NA, deviance = rss))
    }
    
    tmp <- gen.fullrank.g(df            = data.node,
                          adj.form.true = adj.form.true)
    tmp <- est.cond.eff(df        = tmp$df.fullrank,
                        method    = adj.mthd,
                        form.true = tmp$adj.form.true.updated,
                        type.var  = type.var)
    
    est.cond.eff.1 <- tmp$pred.A.1
    est.cond.eff.0 <- tmp$pred.A.0
    
  } 
  
  # Calculating the estimator  
  mu.1 <- mean(A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)
  mu.0 <- mean((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)

  avg.trt.effct <- mu.1 - mu.0
  
  list(label = avg.trt.effct, deviance = rss)
}
