calc.dr.trt.eff <- function(data.node, x, ux,
                            prop.sc, propsc.mod.loc, propsc.mthd, propsc.form.true, whl.propsc,                      # propsc.mthd, propsc.form.true are placeholders for insplit models
                            est.cond.eff.0, est.cond.eff.1, adj.mod.loc, adj.mthd, adj.form.true, type.var,
                            whl.est.cond.eff.0, whl.est.cond.eff.1, avg.trt.effct){     # adj.mthd, adj.form.true are placeholders for insplit models
  
  # calculates ipw mean response difference between treatment and control for each level of x 
  # when x is categorical and the outcome is continuous.
  # data.node: a data frame containing observations at the node.
  # x:         vector of the \texttt{x} values at the node.
  # ux:        unique levels of the categorical covariate x at the node.
  
  # data.node, x, ux,
  # prop.sc            = prop.sc
  # propsc.mod.loc     = parms$propsc.mod.loc
  # propsc.mthd        = parms$propsc.mthd
  # propsc.form.true   = parms$propsc.form.true
  # whl.propsc         = whl.propsc
  # est.cond.eff.0     = est.cond.eff.0
  # est.cond.eff.1     = est.cond.eff.1
  # adj.mod.loc        = parms$adj.mod.loc
  # adj.mthd           = parms$adj.mthd
  # adj.form.true      = parms$adj.form.true
  # type.var           = parms$type.var
  # whl.est.cond.eff.0 = whl.est.cond.eff.0
  # whl.est.cond.eff.1 = whl.est.cond.eff.1
  # avg.trt.effct      = parms$avg.trt.effct
  
  trt.eff <- NULL
  
  # Crude mean response difference between treatment and control for each level of x
  for (i in 1:length(ux)){
    
    data.sub.node      <- data.node[x == ux[i], ]
    
    # when there is only treated/untreated unit in the split,
    # use grand average treatment effect
    if (length(unique(data.sub.node$A)) > 1) {
      
      if (propsc.mod.loc != "split") {
        prop.sc.sub        <- prop.sc[x == ux[i]]
        
      } else {
        
        data.sub.node.noy <- data.sub.node[, !colnames(data.sub.node) %in% c("Y")]
        tmp <- gen.fullrank.ipw(df.noy           = data.sub.node.noy, 
                                propsc.form.true = propsc.form.true)
        tmp <- withWarnings(est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                                        method    = propsc.mthd,
                                        form.true = tmp$propsc.form.true.updated))
        
        if (!is.null(tmp$warnings)) {
          
          cond.warnings <- T
          for (warnings.i in 1:length(tmp$warnings)) {
            cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          # when there is only rank-deficient warning, use the local prediction
          if (cond.warnings) {
            prop.sc.sub <- tmp$value$prop.sc
          } else {
            prop.sc.sub <- whl.propsc[x == ux[i]]
          }
          
        } else if (identical(sort(unique(tmp$value$prop.sc)), c(0.1, 0.9))) {
          prop.sc.sub <- whl.propsc[x == ux[i]]
          
        } else {
          prop.sc.sub <- tmp$value$prop.sc
        }
      }
      
      if (adj.mod.loc != "split") {
        est.cond.eff.0.sub <- est.cond.eff.0[x == ux[i]]
        est.cond.eff.1.sub <- est.cond.eff.1[x == ux[i]]
        
      } else {
        
        tmp <- gen.fullrank.g(df            = data.sub.node,
                              adj.form.true = adj.form.true)
        tmp <- withWarnings(est.cond.eff(df        = tmp$df.fullrank,
                                         method    = adj.mthd,
                                         form.true = tmp$adj.form.true.updated, 
                                         type.var  = type.var))
        
        # when the outcome model fitting produces warning
        # use outside model fitting results for conditional means
        if (!is.null(tmp$warnings)) {
          
          # when there is only rank-deficient warning, use the local prediction
          cond.warnings <- T
          for (warnings.i in 1:length(tmp$warnings)) {
            cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          if (cond.warnings) {
            est.cond.eff.1.sub <- tmp$value$pred.A.1
            est.cond.eff.0.sub <- tmp$value$pred.A.0
          } else {
            est.cond.eff.1.sub <- whl.est.cond.eff.1[x == ux[i]]
            est.cond.eff.0.sub <- whl.est.cond.eff.0[x == ux[i]]
          }
          
        } else {
          est.cond.eff.1.sub <- tmp$value$pred.A.1
          est.cond.eff.0.sub <- tmp$value$pred.A.0
        }
      }
      
      mu.1 <- mean(data.sub.node$A * (data.sub.node$Y - est.cond.eff.1.sub) / prop.sc.sub + est.cond.eff.1.sub)
      mu.0 <- mean((1 - data.sub.node$A) * (data.sub.node$Y - est.cond.eff.0.sub) / (1 - prop.sc.sub) + est.cond.eff.0.sub)
      mu.sub <- mu.1 - mu.0
      
    } else { # if (length(unique(data.sub.node$A)) > 1) {
      mu.sub <- avg.trt.effct
    }
    
    trt.eff     <- c(trt.eff, mu.sub)
    
  } # for loop
  
  return(trt.eff)
  
}

# IN SPLIT for categorical covariate NOT WORKING;
# In split for categorical covariate not needed since insplit only needed for continuous covariate
stemp.dr <- function(y, wt, x, parms, continuous){
  
  n <- length(y)
  
  sub.ind      <- y
  sub.x        <- parms$covariates[sub.ind, ]
  A            <- parms$trt[sub.ind]
  Y            <- parms$response[sub.ind]
  data.node    <- data.frame(Y, A, sub.x)
  
  whl.propsc         <- parms$whl.propsc[sub.ind]
  whl.est.cond.eff.0 <- parms$whl.est.cond.eff.0[sub.ind]
  whl.est.cond.eff.1 <- parms$whl.est.cond.eff.1[sub.ind]
  
  # created for categorical splits, sort categories
  est.cond.eff.1 <- NULL
  est.cond.eff.0 <- NULL
  prop.sc <- NULL
  
  # Propensity score models if "out" or "node"
  if (parms$propsc.mod.loc == "out") {
    prop.sc <- whl.propsc
    
  } else if (parms$propsc.mod.loc == "node") {
    
    data.node.noy <- data.node[, !colnames(data.node) %in% c("Y")]
    
    tmp <- gen.fullrank.ipw(df.noy           = data.node.noy, 
                            propsc.form.true = parms$propsc.form.true)
    tmp <- withWarnings(est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                                    method    = parms$propsc.mthd,
                                    form.true = tmp$propsc.form.true.updated))
    
    if (!is.null(tmp$warnings)) {
      
      cond.warnings <- T
      for (warnings.i in 1:length(tmp$warnings)) {
        cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
      }
      
      # when there is only rank-deficient warning, use the local prediction
      if (cond.warnings) {
        prop.sc <- tmp$value$prop.sc
      } else {
        prop.sc <- whl.propsc
      }
      
      # impute with the global propensity scores if unstable estimates
    } else if (identical(sort(unique(tmp$value$prop.sc)), c(0.1, 0.9))) {
      prop.sc <- whl.propsc
      
    }else {
      prop.sc <- tmp$value$prop.sc
    }
  }
  
  # Outcome model if "out" or "node"
  if (parms$adj.mod.loc == "out") {
    est.cond.eff.1 <- whl.est.cond.eff.1
    est.cond.eff.0 <- whl.est.cond.eff.0
    
  } else if (parms$adj.mod.loc == "node") {
    tmp <- gen.fullrank.g(df            = data.node,
                          adj.form.true = parms$adj.form.true)
    tmp <- withWarnings(est.cond.eff(df        = tmp$df.fullrank,
                                     method    = parms$adj.mthd,
                                     form.true = tmp$adj.form.true.updated,
                                     type.var  = parms$type.var))
    
    # when the outcome model fitting produces warning
    # use outside model fitting results for conditional means
    if (!is.null(tmp$warnings)) {
      
      # when there is only rank-deficient warning, use the local prediction
      cond.warnings <- T
      for (warnings.i in 1:length(tmp$warnings)) {
        cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
      }
      
      if (cond.warnings) {
        est.cond.eff.1 <- tmp$value$pred.A.1
        est.cond.eff.0 <- tmp$value$pred.A.0
      } else {
        est.cond.eff.1 <- whl.est.cond.eff.1
        est.cond.eff.0 <- whl.est.cond.eff.0
      }
      
    } else {
      est.cond.eff.1 <- tmp$value$pred.A.1
      est.cond.eff.0 <- tmp$value$pred.A.0
    }
    
  }
  
  if (continuous) {
    
    # Skip the first 10 and last 10 splits
    # goodness <- NULL
    # direction <- NULL
    goodness  <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    # Ensure the model is at least fitted on num.truc.obs obs
    for (i in (parms$num.truc.obs:(n-parms$num.truc.obs))) { 
      
      data.node.l <- data.node[1:i, ]
      data.node.r <- data.node[(i+1):n, ]
      
      # Propensity score model
      if (parms$propsc.mod.loc == "split") {
        
        # left
        data.node.noy.l <- data.node.l[, !colnames(data.node.l) %in% c("Y")]
        tmp <- gen.fullrank.ipw(df.noy           = data.node.noy.l, 
                                propsc.form.true = parms$propsc.form.true)
        tmp <- withWarnings(est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                                        method    = parms$propsc.mthd,
                                        form.true = tmp$propsc.form.true.updated))
        
        # if 1) propensity score model fitting produces warning or 2) only produce thresholded propensity scores
        # use outside model fitting results
        if (!is.null(tmp$warnings)) {
          
          cond.warnings <- T
          for (warnings.i in 1:length(tmp$warnings)) {
            cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          # when there is only rank-deficient warning, use the local prediction
          if (cond.warnings) {
            prop.sc.l <- tmp$value$prop.sc
          } else {
            # since whl.propsc already subsetted and sorted
            prop.sc.l <- whl.propsc[1:i]
          }
          
        } else if (identical(sort(unique(tmp$value$prop.sc)), c(0.1, 0.9))) {
          prop.sc.l <- whl.propsc[1:i]
        } else {
          prop.sc.l <- tmp$value$prop.sc
        }
        
        
        # right
        data.node.noy.r <- data.node.r[, !colnames(data.node.r) %in% c("Y")]
        tmp <- gen.fullrank.ipw(df.noy           = data.node.noy.r, 
                                propsc.form.true = parms$propsc.form.true)
        tmp <- withWarnings(est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                                        method    = parms$propsc.mthd,
                                        form.true = tmp$propsc.form.true.updated))
        
        if (!is.null(tmp$warnings)) {
          
          cond.warnings <- T
          for (warnings.i in 1:length(tmp$warnings)) {
            cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          # when there is only rank-deficient warning, use the local prediction
          if (cond.warnings) {
            prop.sc.r <- tmp$value$prop.sc
          } else {
            # since whl.propsc already subsetted and sorted
            prop.sc.r <- whl.propsc[(i+1):n]
          }
          
        } else if (identical(sort(unique(tmp$value$prop.sc)), c(0.1, 0.9))) { 
          prop.sc.r <- whl.propsc[(i+1):n]
        } else {
          prop.sc.r <- tmp$value$prop.sc
        }
        
        prop.sc <- c(prop.sc.l, prop.sc.r)
      } 
      
      # Outcome model
      if (parms$adj.mod.loc == "split") {
        
        # left
        tmp <- gen.fullrank.g(df            = data.node.l,
                              adj.form.true = parms$adj.form.true)
        tmp <- withWarnings(est.cond.eff(df        = tmp$df.fullrank,
                                         method    = parms$adj.mthd,
                                         form.true = tmp$adj.form.true.updated,
                                         type.var  = parms$type.var))
        
        # when the outcome model fitting produces warning
        # use outside model fitting results for conditional means
        if (!is.null(tmp$warnings)) {
          
          # when there is only rank-deficient warning, use the local prediction
          cond.warnings <- T
          for (warnings.i in 1:length(tmp$warnings)) {
            cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          if (cond.warnings) {
            est.cond.eff.1.l <- tmp$value$pred.A.1
            est.cond.eff.0.l <- tmp$value$pred.A.0
          } else {
            est.cond.eff.0.l <- whl.est.cond.eff.0[1:i]
            est.cond.eff.1.l <- whl.est.cond.eff.1[1:i]
          }
          
        } else {
          est.cond.eff.0.l <- tmp$value$pred.A.0
          est.cond.eff.1.l <- tmp$value$pred.A.1
        }
        
        # right
        tmp <- gen.fullrank.g(df            = data.node.r,
                              adj.form.true = parms$adj.form.true)
        tmp <- withWarnings(est.cond.eff(df        = tmp$df.fullrank,
                                         method    = parms$adj.mthd,
                                         form.true = tmp$adj.form.true.updated,
                                         type.var  = parms$type.var))
        
        if (!is.null(tmp$warnings)) {
          
          # when there is only rank-deficient warning, use the local prediction
          cond.warnings <- T
          for (warnings.i in 1:length(tmp$warnings)) {
            cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          if (cond.warnings) {
            est.cond.eff.0.r <- tmp$value$pred.A.0
            est.cond.eff.1.r <- tmp$value$pred.A.1
          } else {
            est.cond.eff.0.r <- whl.est.cond.eff.0[(i+1):n]
            est.cond.eff.1.r <- whl.est.cond.eff.1[(i+1):n]
          }
          
        } else {
          est.cond.eff.0.r <- tmp$value$pred.A.0
          est.cond.eff.1.r <- tmp$value$pred.A.1
        }
        
        est.cond.eff.1 <- c(est.cond.eff.1.l, est.cond.eff.1.r)
        est.cond.eff.0 <- c(est.cond.eff.0.l, est.cond.eff.0.r)
        
      }
      
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
      
      # variance smaller if both are glm models
      if ((parms$adj.mthd == "GLM") & (parms$propsc.mthd == "GLM")) {
        var <- (mean(I_i^2) - (p.r * (mu.1l - mu.0l) + p.l * (mu.1r - mu.0r))^2 / (p.l * p.r)) / n
        
        # Skip this loop if there is negative estimated variance
        if (var < 0) {
          next
        }
      } else {
        var <- mean(I_i^2) / n
      }
      
      goodness[i] <- (t.dr / sqrt(var))^2
      direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))   
      
    }
    
    # if (is.null(est.cond.eff.1)) { 
    #   
    #   if (is.null(prop.sc)) { # Both models in split
    #     
    #     # Ensure the model is at least fitted on num.truc.obs obs
    #     for (i in (parms$num.truc.obs:(n-parms$num.truc.obs))) { 
    #       
    #       data.node.l <- data.node[1:i, ]
    #       data.node.r <- data.node[(i+1):n, ]
    #       
    #       # Propensity score model
    #       data.node.noy.l <- data.node.l[, !colnames(data.node.l) %in% c("Y")]
    #       tmp <- gen.fullrank.ipw(df.noy           = data.node.noy.l, 
    #                               propsc.form.true = parms$propsc.form.true)
    #       tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
    #                          method    = parms$propsc.mthd,
    #                          form.true = tmp$propsc.form.true.updated)
    #       prop.sc.l <- tmp$prop.sc
    #       
    #       data.node.noy.r <- data.node.r[, !colnames(data.node.r) %in% c("Y")]
    #       tmp <- gen.fullrank.ipw(df.noy           = data.node.noy.r, 
    #                               propsc.form.true = parms$propsc.form.true)
    #       tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
    #                          method    = parms$propsc.mthd,
    #                          form.true = tmp$propsc.form.true.updated)
    #       prop.sc.r <- tmp$prop.sc
    #       
    #       prop.sc <- c(prop.sc.l, prop.sc.r)
    #       
    #       # Outcome model
    #       tmp <- gen.fullrank.g(df            = data.node.l,
    #                             adj.form.true = parms$adj.form.true)
    #       tmp <- est.cond.eff(df        = tmp$df.fullrank,
    #                           method    = parms$adj.mthd,
    #                           form.true = tmp$adj.form.true.updated,
    #                           type.var  = parms$type.var)
    #       
    #       est.cond.eff.1.l <- tmp$pred.A.1
    #       est.cond.eff.0.l <- tmp$pred.A.0
    #       
    #       tmp <- gen.fullrank.g(df            = data.node.r,
    #                             adj.form.true = parms$adj.form.true)
    #       tmp <- est.cond.eff(df        = tmp$df.fullrank,
    #                           method    = parms$adj.mthd,
    #                           form.true = tmp$adj.form.true.updated,
    #                           type.var  = parms$type.var)
    #       
    #       est.cond.eff.1.r <- tmp$pred.A.1
    #       est.cond.eff.0.r <- tmp$pred.A.0
    #       
    #       est.cond.eff.1 <- c(est.cond.eff.1.l, est.cond.eff.1.r)
    #       est.cond.eff.0 <- c(est.cond.eff.0.l, est.cond.eff.0.r)
    #       
    #       # Calculating the estimator
    #       mu.1l <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[1:i])
    #       mu.0l <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[1:i])
    #       mu.1r <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[(i+1):n])
    #       mu.0r <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[(i+1):n])
    #       
    #       t.dr <- (mu.1l - mu.0l) - (mu.1r - mu.0r)
    #       p.l <- i/n
    #       p.r <- (n-i)/n
    #       I_i <- (x %in% x[1:i]) / p.l * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
    #                                         ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
    #         (x %in% x[(i+1):n]) / p.r * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
    #                                        ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
    #         t.dr
    #       
    #       var <- (mean(I_i^2) - (p.r * (mu.1l - mu.0l) + p.l * (mu.1r - mu.0r))^2 / (p.l * p.r)) / n
    #       
    #       # Skip this loop if there is negative estimated variance
    #       if (var < 0) {
    #         next
    #       }
    #       
    #       goodness[i] <- (t.dr / sqrt(var))^2
    #       direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))      
    #     }
    #     
    #   } else { # prop.sc not in split, outcome model in split
    #     
    #     # Ensure the model is at least fitted on num.truc.obs obs
    #     for (i in (parms$num.truc.obs:(n-parms$num.truc.obs))) { 
    #       
    #       data.node.l <- data.node[1:i, ]
    #       data.node.r <- data.node[(i+1):n, ]
    #       
    #       # Outcome model
    #       tmp <- gen.fullrank.g(df            = data.node.l,
    #                             adj.form.true = parms$adj.form.true)
    #       tmp <- est.cond.eff(df        = tmp$df.fullrank,
    #                           method    = parms$adj.mthd,
    #                           form.true = tmp$adj.form.true.updated,
    #                           type.var  = parms$type.var)
    #       
    #       est.cond.eff.1.l <- tmp$pred.A.1
    #       est.cond.eff.0.l <- tmp$pred.A.0
    #       
    #       tmp <- gen.fullrank.g(df            = data.node.r,
    #                             adj.form.true = parms$adj.form.true)
    #       tmp <- est.cond.eff(df        = tmp$df.fullrank,
    #                           method    = parms$adj.mthd,
    #                           form.true = tmp$adj.form.true.updated,
    #                           type.var  = parms$type.var)
    #       
    #       est.cond.eff.1.r <- tmp$pred.A.1
    #       est.cond.eff.0.r <- tmp$pred.A.0
    #       
    #       est.cond.eff.1 <- c(est.cond.eff.1.l, est.cond.eff.1.r)
    #       est.cond.eff.0 <- c(est.cond.eff.0.l, est.cond.eff.0.r)
    #       
    #       # Calculating the estimator
    #       mu.1l <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[1:i])
    #       mu.0l <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[1:i])
    #       mu.1r <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[(i+1):n])
    #       mu.0r <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[(i+1):n])
    #       
    #       t.dr <- (mu.1l - mu.0l) - (mu.1r - mu.0r)
    #       p.l <- i/n
    #       p.r <- (n-i)/n
    #       I_i <- (x %in% x[1:i]) / p.l * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
    #                                         ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
    #         (x %in% x[(i+1):n]) / p.r * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
    #                                        ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
    #         t.dr
    #       
    #       var <- (mean(I_i^2) - (p.r * (mu.1l - mu.0l) + p.l * (mu.1r - mu.0r))^2 / (p.l * p.r)) / n
    #       
    #       # Skip this loop if there is negative estimated variance
    #       if (var < 0) {
    #         next
    #       }
    #       
    #       goodness[i] <- (t.dr / sqrt(var))^2
    #       direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))      
    #     }
    #     
    #     
    #   }
    #   
    #   
    # } else { # Outcome model in node or outside
    #   
    #   if (is.null(prop.sc)) { # propensity score model in split
    #     
    #     # Ensure the model is at least fitted on num.truc.obs obs
    #     for (i in (parms$num.truc.obs:(n-parms$num.truc.obs))) { 
    #       
    #       data.node.l <- data.node[1:i, ]
    #       data.node.r <- data.node[(i+1):n, ]
    #       
    #       # Propensity score model
    #       data.node.noy.l <- data.node.l[, !colnames(data.node.l) %in% c("Y")]
    #       tmp <- gen.fullrank.ipw(df.noy           = data.node.noy.l, 
    #                               propsc.form.true = parms$propsc.form.true)
    #       tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
    #                          method    = parms$propsc.mthd,
    #                          form.true = tmp$propsc.form.true.updated)
    #       prop.sc.l <- tmp$prop.sc
    #       
    #       data.node.noy.r <- data.node.r[, !colnames(data.node.r) %in% c("Y")]
    #       tmp <- gen.fullrank.ipw(df.noy           = data.node.noy.r, 
    #                               propsc.form.true = parms$propsc.form.true)
    #       tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
    #                          method    = parms$propsc.mthd,
    #                          form.true = tmp$propsc.form.true.updated)
    #       prop.sc.r <- tmp$prop.sc
    #       
    #       prop.sc <- c(prop.sc.l, prop.sc.r)
    #       
    #       # Calculating the estimator
    #       mu.1l <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[1:i])
    #       mu.0l <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[1:i])
    #       mu.1r <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[(i+1):n])
    #       mu.0r <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[(i+1):n])
    #       
    #       t.dr <- (mu.1l - mu.0l) - (mu.1r - mu.0r)
    #       p.l <- i/n
    #       p.r <- (n-i)/n
    #       I_i <- (x %in% x[1:i]) / p.l * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
    #                                         ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
    #         (x %in% x[(i+1):n]) / p.r * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
    #                                        ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
    #         t.dr
    #       
    #       var <- (mean(I_i^2) - (p.r * (mu.1l - mu.0l) + p.l * (mu.1r - mu.0r))^2 / (p.l * p.r)) / n
    #       
    #       # Skip this loop if there is negative estimated variance
    #       if (var < 0) {
    #         next
    #       }
    #       
    #       goodness[i] <- (t.dr / sqrt(var))^2
    #       direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))      
    #     }
    #     
    #   } else { # Both models outside split
    #     
    #     # Ensure the model is at least fitted on num.truc.obs obs
    #     for (i in (parms$num.truc.obs:(n-parms$num.truc.obs))) { 
    #       
    #       # Calculating the estimator
    #       mu.1l <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[1:i])
    #       mu.0l <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[1:i])
    #       mu.1r <- mean((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)[(i+1):n])
    #       mu.0r <- mean(((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)[(i+1):n])
    #       
    #       t.dr <- (mu.1l - mu.0l) - (mu.1r - mu.0r)
    #       p.l <- i/n
    #       p.r <- (n-i)/n
    #       I_i <- (x %in% x[1:i]) / p.l * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
    #                                         ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
    #         (x %in% x[(i+1):n]) / p.r * ((A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1) - 
    #                                        ((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)) - 
    #         t.dr
    #       
    #       var <- (mean(I_i^2) - (p.r * (mu.1l - mu.0l) + p.l * (mu.1r - mu.0r))^2 / (p.l * p.r)) / n
    #       
    #       # Skip this loop if there is negative estimated variance
    #       if (var < 0) {
    #         next
    #       }
    #       
    #       goodness[i] <- (t.dr / sqrt(var))^2
    #       direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))      
    #     }
    #     
    #   } # Both models outside split
    #   
    # } # if (is.null(est.cond.eff.1)) {} else{}
    
  } else {
    
    ux <- sort(unique(x))
    
    # Order the levels of the categorical covariates by their treatment effect
    # Currently code only works for outside node and innode model fitting, not for in split
    trt.eff.dr <- calc.dr.trt.eff(data.node, x, ux,
                                  prop.sc            = prop.sc,
                                  propsc.mod.loc     = parms$propsc.mod.loc,
                                  propsc.mthd        = parms$propsc.mthd,
                                  propsc.form.true   = parms$propsc.form.true,
                                  whl.propsc         = whl.propsc,
                                  est.cond.eff.0     = est.cond.eff.0,
                                  est.cond.eff.1     = est.cond.eff.1,
                                  adj.mod.loc        = parms$adj.mod.loc,
                                  adj.mthd           = parms$adj.mthd,
                                  adj.form.true      = parms$adj.form.true,
                                  type.var           = parms$type.var,
                                  whl.est.cond.eff.0 = whl.est.cond.eff.0,
                                  whl.est.cond.eff.1 = whl.est.cond.eff.1,
                                  avg.trt.effct      = parms$avg.trt.effct)
    
    
    ord.ux <- order(trt.eff.dr)
    goodness  <- rep(0, length(ux) - 1)
    
    # Splits
    for (i in 1:(length(ux) - 1)){
      
      ind.l <- x %in% ux[ord.ux[1:i]]
      ind.r <- x %in% ux[ord.ux[(i+1):length(ux)]] 
      
      data.node.l <- data.node[ind.l, ]
      data.node.r <- data.node[ind.r, ]
      
      # Propensity score model
      if (parms$propsc.mod.loc == "split") {
        
        # left
        data.node.noy.l <- data.node.l[, !colnames(data.node.l) %in% c("Y")]
        tmp <- gen.fullrank.ipw(df.noy           = data.node.noy.l, 
                                propsc.form.true = parms$propsc.form.true)
        tmp <- withWarnings(est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                                        method    = parms$propsc.mthd,
                                        form.true = tmp$propsc.form.true.updated))
        
        # if 1) propensity score model fitting produces warning or 2) only produce thresholded propensity scores
        # use outside model fitting results
        if (!is.null(tmp$warnings)) {
          
          cond.warnings <- T
          for (warnings.i in 1:length(tmp$warnings)) {
            cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          # when there is only rank-deficient warning, use the local prediction
          if (cond.warnings) {
            prop.sc.l <- tmp$value$prop.sc
          } else {
            # since whl.propsc already subsetted and sorted
            prop.sc.l <- whl.propsc[ind.l]
          }
          
        } else if (identical(sort(unique(tmp$value$prop.sc)), c(0.1, 0.9))) {
          prop.sc.l <- whl.propsc[ind.l]
          
        } else {
          prop.sc.l <- tmp$value$prop.sc
        }
        
        # right
        data.node.noy.r <- data.node.r[, !colnames(data.node.r) %in% c("Y")]
        tmp <- gen.fullrank.ipw(df.noy           = data.node.noy.r, 
                                propsc.form.true = parms$propsc.form.true)
        tmp <- withWarnings(est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                                        method    = parms$propsc.mthd,
                                        form.true = tmp$propsc.form.true.updated))
        
        if (!is.null(tmp$warnings)) {
          
          cond.warnings <- T
          for (warnings.i in 1:length(tmp$warnings)) {
            cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          # when there is only rank-deficient warning, use the local prediction
          if (cond.warnings) {
            prop.sc.r <- tmp$value$prop.sc
          } else {
            # since whl.propsc already subsetted and sorted
            prop.sc.r <- whl.propsc[ind.r]
          }
          
        } else if (identical(sort(unique(tmp$value$prop.sc)), c(0.1, 0.9))) {
          prop.sc.r <- whl.propsc[ind.r]
          
        } else {
          prop.sc.r <- tmp$value$prop.sc
        }
        
        prop.sc <- rep(0, n)
        prop.sc[ind.l] <- prop.sc.l
        prop.sc[ind.r] <- prop.sc.r
        
      }
      
      # Outcome model
      if (parms$adj.mod.loc == "split") {
        # left
        tmp <- gen.fullrank.g(df            = data.node.l,
                              adj.form.true = parms$adj.form.true)
        tmp <- withWarnings(est.cond.eff(df        = tmp$df.fullrank,
                                         method    = parms$adj.mthd,
                                         form.true = tmp$adj.form.true.updated,
                                         type.var  = parms$type.var))
        
        # when the outcome model fitting produces warning
        # use outside model fitting results for conditional means
        if (!is.null(tmp$warnings)) {
          
          # when there is only rank-deficient warning, use the local prediction
          cond.warnings <- T
          for (warnings.i in 1:length(tmp$warnings)) {
            cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          if (cond.warnings) {
            est.cond.eff.1.l <- tmp$value$pred.A.1
            est.cond.eff.0.l <- tmp$value$pred.A.0
          } else {
            est.cond.eff.0.l <- whl.est.cond.eff.0[ind.l]
            est.cond.eff.1.l <- whl.est.cond.eff.1[ind.l]
          }
          
        } else {
          est.cond.eff.0.l <- tmp$value$pred.A.0
          est.cond.eff.1.l <- tmp$value$pred.A.1
        }
        
        # right
        tmp <- gen.fullrank.g(df            = data.node.r,
                              adj.form.true = parms$adj.form.true)
        tmp <- withWarnings(est.cond.eff(df        = tmp$df.fullrank,
                                         method    = parms$adj.mthd,
                                         form.true = tmp$adj.form.true.updated,
                                         type.var  = parms$type.var))
        
        if (!is.null(tmp$warnings)) {
          
          # when there is only rank-deficient warning, use the local prediction
          cond.warnings <- T
          for (warnings.i in 1:length(tmp$warnings)) {
            cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          if (cond.warnings) {
            est.cond.eff.0.r <- tmp$value$pred.A.0
            est.cond.eff.1.r <- tmp$value$pred.A.1
          } else {
            est.cond.eff.0.r <- whl.est.cond.eff.0[ind.r]
            est.cond.eff.1.r <- whl.est.cond.eff.1[ind.r]
          }
          
        } else {
          est.cond.eff.0.r <- tmp$value$pred.A.0
          est.cond.eff.1.r <- tmp$value$pred.A.1
        }
        
        est.cond.eff.1 <- rep(0, n)
        est.cond.eff.0 <- rep(0, n)
        
        est.cond.eff.1[ind.l] <- est.cond.eff.1.l
        est.cond.eff.1[ind.r] <- est.cond.eff.1.r
        est.cond.eff.0[ind.l] <- est.cond.eff.0.l
        est.cond.eff.0[ind.r] <- est.cond.eff.0.r
        
      }
      
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
      
      # variance smaller if both are glm models
      if ((parms$adj.mthd == "GLM") & (parms$propsc.mthd == "GLM")) {
        var <- (mean(I_i^2) - (p.r * (mu.1l - mu.0l) + p.l * (mu.1r - mu.0r))^2 / (p.l * p.r)) / n
        
        # Skip this loop if there is negative estimated variance
        if (var < 0) {
          next
        }
      } else {
        var <- mean(I_i^2) / n
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
  
  whl.propsc         <- parms$whl.propsc[sub.ind]
  whl.est.cond.eff.0 <- parms$whl.est.cond.eff.0[sub.ind]
  whl.est.cond.eff.1 <- parms$whl.est.cond.eff.1[sub.ind]
  
  # Need to move this part of code here for return value earlier
  # TO DO: Weird how this does not work when rss is used.
  # not crucial to the performance 
  wmean <- sum(Y*wt)/sum(wt)
  rss <- sum(wt*(Y-wmean)^2)
  
  # when there are both treated and untreated unit in split, use avg.trt.effct
  if (length(unique(data.node$A)) == 1) {
    return(list(label = parms$avg.trt.effct, deviance = rss))
    
  } else {
    
    if (parms$propsc.mod.loc == "out") {          
      prop.sc <- whl.propsc
      
    } else {
      
      data.node.noy <- data.node[, !colnames(data.node) %in% c("Y")]
      tmp <- gen.fullrank.ipw(df.noy           = data.node.noy, 
                              propsc.form.true = parms$propsc.form.true)
      tmp <- withWarnings(est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                                      method    = parms$propsc.mthd,
                                      form.true = tmp$propsc.form.true.updated))
      
      if (!is.null(tmp$warnings)) {
        
        cond.warnings <- T
        for (warnings.i in 1:length(tmp$warnings)) {
          cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
        }
        
        # when there is only rank-deficient warning, use the local prediction
        if (cond.warnings) {
          prop.sc <- tmp$value$prop.sc
        } else {
          prop.sc <- whl.propsc
        }
        
      } else if (identical(sort(unique(tmp$value$prop.sc)), c(0.1, 0.9))) {
        prop.sc <- whl.propsc
        
      } else {
        prop.sc <- tmp$value$prop.sc
      }
    }
    
    if (parms$adj.mod.loc == "out") {
      est.cond.eff.1 <- whl.est.cond.eff.1
      est.cond.eff.0 <- whl.est.cond.eff.0
      
    } else {
      tmp <- gen.fullrank.g(df            = data.node,
                            adj.form.true = parms$adj.form.true)
      tmp <- withWarnings(est.cond.eff(df        = tmp$df.fullrank,
                                       method    = parms$adj.mthd,
                                       form.true = tmp$adj.form.true.updated,
                                       type.var  = parms$type.var))
      
      # when the outcome model fitting produces warning
      # use outside model fitting results for conditional means
      if (!is.null(tmp$warnings)) {
        
        # when there is only rank-deficient warning, use the local prediction
        cond.warnings <- T
        for (warnings.i in 1:length(tmp$warnings)) {
          cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
        }
        
        if (cond.warnings) {
          est.cond.eff.1 <- tmp$value$pred.A.1
          est.cond.eff.0 <- tmp$value$pred.A.0
        } else {
          est.cond.eff.1 <- whl.est.cond.eff.1
          est.cond.eff.0 <- whl.est.cond.eff.0
        }
        
      } else { # if (!is.null(tmp$warnings)) {
        est.cond.eff.1 <- tmp$value$pred.A.1
        est.cond.eff.0 <- tmp$value$pred.A.0
      }
    }
    
  } #   if (length(unique(data.node.noy$A)) == 1) { else
  
  # Calculating the estimator  
  mu.1 <- mean(A * (Y - est.cond.eff.1) / prop.sc + est.cond.eff.1)
  mu.0 <- mean((1 - A) * (Y - est.cond.eff.0) / (1 - prop.sc) + est.cond.eff.0)
  
  avg.trt.effct <- mu.1 - mu.0
  
  list(label = avg.trt.effct, deviance = rss)
}
