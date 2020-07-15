calc.g.trt.eff <- function(data.node, x, ux, est.cond.eff.0, est.cond.eff.1, 
                           adj.mod.loc, adj.mthd, adj.form.true, type.var, 
                           whl.est.cond.eff.0, whl.est.cond.eff.1){
  
  # calculates ipw mean response difference between treatment and control for each level of x 
  # when x is categorical and the outcome is continuous.
  # data.node: a data frame containing observations at the node.
  # x:         vector of the \texttt{x} values at the node.
  # ux:        unique levels of the categorical covariate x at the node.
  
  trt.eff <- NULL
  
  # Crude mean response difference between treatment and control for each level of x
  for (i in 1:length(ux)){
    
    data.sub.node      <- data.node[x == ux[i], ]
    
    if (adj.mod.loc == "split") { # when the conditional means are estimated within split
      
      if (length(unique(data.sub.node$A)) > 1) {
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
        
        mu.sub <- mean(est.cond.eff.1.sub) - mean(est.cond.eff.0.sub)
      } else { # when there is only treated/untreated unit in that level, use avg.trt.effect for imputation
        
        est.cond.eff.1.sub <- whl.est.cond.eff.1[x == ux[i]]
        est.cond.eff.0.sub <- whl.est.cond.eff.0[x == ux[i]]
        mu.sub <- mean(est.cond.eff.1.sub) - mean(est.cond.eff.0.sub)
      }
      
    } else { # if (is.null(est.cond.eff.0)) { 
      
      est.cond.eff.0.sub <- est.cond.eff.0[x == ux[i]]
      est.cond.eff.1.sub <- est.cond.eff.1[x == ux[i]]
      mu.sub <- mean(est.cond.eff.1.sub) - mean(est.cond.eff.0.sub)
    }
    
    trt.eff     <- c(trt.eff, mu.sub)
    
  } # for loop
  
  return(trt.eff)
  
}



gen.fullrank.g <- function(df, adj.form.true) {
  
  # remove the factor column if the column has only one level
  df.fullrank <- data.frame(Y = df$Y, 
                            A = df$A, 
                            df[, !colnames(df) %in% c("A", "Y")][, sapply(df[, !colnames(df) %in% c("A", "Y")], 
                                                                          function(col) length(unique(col))) > 1])
  # re-factor the factor columns if there is some levels in the factor column missing
  # drop the missing factor
  df.fullrank[] <- lapply(df.fullrank, 
                          function(x) if(is.factor(x)) factor(x) else x)
  
  if (ncol(df.fullrank) < ncol(df) & (!is.null(adj.form.true))) {
    adj.terms <- gsub("Y ~ ", "", adj.form.true)
    adj.terms <- base::strsplit(adj.terms, " + ", fixed = T)[[1]]
    
    # Need to take care of the number of times the covariates appear in each term in the model;
    # if there is one term missing in the interaction, need to delete the interaction term
    adj.terms <- adj.terms[rowSums(sapply(names(df.fullrank), function(covariate) grepl(covariate, adj.terms))) >= 
                             rowSums(sapply(names(df), function(covariate) grepl(covariate, adj.terms)))]
    
    adj.form.true.updated <- paste("Y ~ ", paste(adj.terms, collapse= " + "))
  } else {
    adj.form.true.updated <- adj.form.true
  }
  
  return(list(df.fullrank           = df.fullrank,
              adj.form.true.updated = adj.form.true.updated))
}



# Continuous Outcome
stemp.g <- function(y, wt, x, parms, continuous){
  
  n <- length(y)
  
  # Finding observations in each node
  sub.ind   <- y
  sub.x     <- parms$covariates[sub.ind, ]
  A         <- parms$trt[sub.ind]
  Y         <- parms$response[sub.ind]
  data.node <- data.frame(Y, A, sub.x)
  
  whl.est.cond.eff.1 <- parms$whl.est.cond.eff.1[sub.ind]
  whl.est.cond.eff.0 <- parms$whl.est.cond.eff.0[sub.ind]
  whl.w              <- parms$whl.w[sub.ind, ]
  whl.var.rb         <- parms$whl.var.rb
  
  # created for categorical splits, sort categories
  est.cond.eff.1 <- NULL
  est.cond.eff.0 <- NULL
  
  if (parms$adj.mod.loc == "out") {
    
    est.cond.eff.1 <- whl.est.cond.eff.1
    est.cond.eff.0 <- whl.est.cond.eff.0
    w              <- whl.w
    var.rb         <- whl.var.rb
    
  } else if (parms$adj.mod.loc == "node") {
    
    tmp <- gen.fullrank.g(df            = data.node,
                          adj.form.true = parms$form.true)
    tmp <- withWarnings(est.cond.eff(df        = tmp$df.fullrank,
                                     method    = parms$adj.mthd,
                                     form.true = tmp$adj.form.true.updated,
                                     type.var  = parms$type.var))
    var.rb <- tmp$value$var.rb
    w      <- tmp$value$w
    
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
  
  if (continuous){
    
    # Skip the first 10 and last 10 splits, technically not necessary 
    # as minsplit takes care of that
    # goodness <- NULL
    # direction <- NULL
    goodness <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    if (parms$adj.mod.loc == "split"){ # model in split
      
      for (i in (parms$num.truc.obs:(n-parms$num.truc.obs))) {
        
        # if all observations are treated or not treated in the son nodes,
        # skip this loop
        if ((length(unique(data.node$A[1:i])) == 1) | (length(unique(data.node$A[(i+1):n])) == 1)){
          next
        }
        
        tmp.l <- gen.fullrank.g(df            = data.node[1:i, ],
                                adj.form.true = parms$form.true)
        tmp.r <- gen.fullrank.g(df            = data.node[(i+1):n, ],
                                adj.form.true = parms$form.true)
        
        tmp.l <- withWarnings(est.cond.eff(df        = tmp.l$df.fullrank,
                                           method    = parms$adj.mthd,
                                           form.true = tmp.l$adj.form.true.updated,
                                           type.var  = parms$type.var))
        tmp.r <- withWarnings(est.cond.eff(df        = tmp.r$df.fullrank,
                                           method    = parms$adj.mthd,
                                           form.true = tmp.r$adj.form.true.updated,
                                           type.var  = parms$type.var))
        var.rb.l <- tmp.l$value$var.rb
        var.rb.r <- tmp.r$value$var.rb
        w.l <- tmp.l$value$w
        w.r <- tmp.r$value$w
        
        # when the outcome model fitting produces warning
        # use outside model fitting results for conditional means
        if (!is.null(tmp.l$warnings)) {
          
          # when there is only rank-deficient warning, use the local prediction
          cond.warnings <- T
          for (warnings.i in 1:length(tmp.l$warnings)) {
            cond.warnings <- cond.warnings & (tmp.l$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          if (cond.warnings) {
            est.cond.eff.1.l <- tmp.l$value$pred.A.1
            est.cond.eff.0.l <- tmp.l$value$pred.A.0
          } else {
            est.cond.eff.0.l <- whl.est.cond.eff.0[1:i]
            est.cond.eff.1.l <- whl.est.cond.eff.1[1:i]
          }
          
        } else {
          est.cond.eff.0.l <- tmp.l$value$pred.A.0
          est.cond.eff.1.l <- tmp.l$value$pred.A.1
        }
        
        if (!is.null(tmp.r$warnings)) {
          
          # when there is only rank-deficient warning, use the local prediction
          cond.warnings <- T
          for (warnings.i in 1:length(tmp.r$warnings)) {
            cond.warnings <- cond.warnings & (tmp.r$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          if (cond.warnings) {
            est.cond.eff.0.r <- tmp.r$value$pred.A.0
            est.cond.eff.1.r <- tmp.r$value$pred.A.1
          } else {
            est.cond.eff.0.r <- whl.est.cond.eff.0[(i+1):n]
            est.cond.eff.1.r <- whl.est.cond.eff.1[(i+1):n]
          }
          
        } else {
          est.cond.eff.0.r <- tmp.r$value$pred.A.0
          est.cond.eff.1.r <- tmp.r$value$pred.A.1
        }
        
        mu.1l <- mean(est.cond.eff.1.l)
        mu.0l <- mean(est.cond.eff.0.l)
        mu.1r <- mean(est.cond.eff.1.r)
        mu.0r <- mean(est.cond.eff.0.r)
        
        if (parms$type.var == "cont") {
          
          g.b.l.1 <- data.frame(w.l) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.r.1 <- data.frame(w.r) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.l.0 <- data.frame(w.l) %>% mutate(A = 0) %>% apply(2, mean)
          g.b.r.0 <- data.frame(w.r) %>% mutate(A = 0) %>% apply(2, mean)
          
        } else if (parms$type.var == "bin") {
          
          w.1l <- w.l %>% mutate(A = 1)
          w.0l <- w.l %>% mutate(A = 0)
          w.1r <- w.r %>% mutate(A = 1)
          w.0r <- w.r %>% mutate(A = 0)
          
          g.b.l.1 <- apply(w.1l * est.cond.eff.1.l * (1 - est.cond.eff.1.l), 2, mean)
          g.b.r.1 <- apply(w.1r * est.cond.eff.1.r * (1 - est.cond.eff.1.r), 2, mean)
          g.b.l.0 <- apply(w.0l * est.cond.eff.0.l * (1 - est.cond.eff.0.l), 2, mean)
          g.b.r.0 <- apply(w.0r * est.cond.eff.0.r * (1 - est.cond.eff.0.r), 2, mean)
          
        }
        
        left.p.1 <- as.matrix(data.frame(est.cond.eff.1.l))
        left.p.0 <- as.matrix(data.frame(est.cond.eff.0.l))
        right.p.1 <- as.matrix(data.frame(est.cond.eff.1.r))
        right.p.0 <- as.matrix(data.frame(est.cond.eff.0.r))
        
        # Calculate variance estimators
        var.1l <- g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/i^2 * sum( (left.p.1 - mu.1l)^2 )
        var.1r <- g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/(n-i)^2 * sum( (right.p.1 - mu.1r)^2 )
        var.0l <- g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/i^2 * sum( (left.p.0 - mu.0l)^2 )
        var.0r <- g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/(n-i)^2 * sum( (right.p.0 - mu.0r)^2 )
        
        var.1l <- as.numeric(var.1l)
        var.1r <- as.numeric(var.1r)
        var.0l <- as.numeric(var.0l)
        var.0r <- as.numeric(var.0r)
        
        goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
        direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r))) 
        
      }
      
    } else { # model outside split
      
      # if rank deficient regression exists, need to take care of the dimension of the design matrix
      # if (ncol(var.rb) < ncol(w)) {
      #   w <- w %>%
      #     select(colnames(var.rb))
      # }
      
      for (i in (parms$num.truc.obs:(n-parms$num.truc.obs))) {
        
        mu.1l <- mean(est.cond.eff.1[1:i])
        mu.0l <- mean(est.cond.eff.0[1:i])
        mu.1r <- mean(est.cond.eff.1[(i+1):n])
        mu.0r <- mean(est.cond.eff.0[(i+1):n])
        
        # Need to change if the outcome is not continuous
        # involve the derivative of the linear function on the outcome scale
        if (parms$type.var == "cont") {
          
          g.b.l.1 <- data.frame(w[1:i, ]) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.r.1 <- data.frame(w[(i+1):n, ]) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.l.0 <- data.frame(w[1:i, ]) %>% mutate(A = 0) %>% apply(2, mean)
          g.b.r.0 <- data.frame(w[(i+1):n, ]) %>% mutate(A = 0) %>% apply(2, mean)
          
        } else if (parms$type.var == "bin") {
          
          w.1 <- w %>% mutate(A = 1)
          w.0 <- w %>% mutate(A = 0)
          
          tmp.1 <- w.1 * est.cond.eff.1 * (1 - est.cond.eff.1)
          tmp.0 <- w.0 * est.cond.eff.0 * (1 - est.cond.eff.0)
          
          g.b.l.1 <- apply(tmp.1[1:i, ], 2, mean)
          g.b.r.1 <- apply(tmp.1[(i+1):n, ], 2, mean)
          g.b.l.0 <- apply(tmp.0[1:i, ], 2, mean)
          g.b.r.0 <- apply(tmp.0[(i+1):n, ], 2, mean)
          
        }
        
        left.p.1 <- as.matrix(data.frame(est.cond.eff.1[1:i]))
        left.p.0 <- as.matrix(data.frame(est.cond.eff.0[1:i]))
        right.p.1 <- as.matrix(data.frame(est.cond.eff.1[(i+1):n]))
        right.p.0 <- as.matrix(data.frame(est.cond.eff.0[(i+1):n]))
        
        # Calculate variance estimators
        var.1l <- g.b.l.1 %*% var.rb %*% g.b.l.1 + 1/i^2 * sum( (left.p.1 - mu.1l)^2 )
        var.1r <- g.b.r.1 %*% var.rb %*% g.b.r.1 + 1/(n-i)^2 * sum( (right.p.1 - mu.1r)^2 )
        var.0l <- g.b.l.0 %*% var.rb %*% g.b.l.0 + 1/i^2 * sum( (left.p.0 - mu.0l)^2 )
        var.0r <- g.b.r.0 %*% var.rb %*% g.b.r.0 + 1/(n-i)^2 * sum( (right.p.0 - mu.0r)^2 )
        
        var.1l <- as.numeric(var.1l)
        var.1r <- as.numeric(var.1r)
        var.0l <- as.numeric(var.0l)
        var.0r <- as.numeric(var.0r)
        
        goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
        direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r))) 
        
      }
      
    }
    
    
  } else{
    
    ux <- sort(unique(x))
    
    # if (parms$adj.mod.loc == "split") { # when the conditional mean model is fitted within split, 
    # need to fit conditional mean models within each level
    trt.eff.g <- calc.g.trt.eff(data.node, x, ux,
                                est.cond.eff.0     = est.cond.eff.0,
                                est.cond.eff.1     = est.cond.eff.1,
                                adj.mod.loc        = parms$adj.mod.loc, 
                                adj.mthd           = parms$adj.mthd,
                                adj.form.true      = parms$form.true, 
                                type.var           = parms$type.var,
                                whl.est.cond.eff.0 = whl.est.cond.eff.0,
                                whl.est.cond.eff.1 = whl.est.cond.eff.1)
    
    # } else { 
    #   
    #   trt.eff.g <- calc.g.trt.eff(data.node, x, ux,
    #                               est.cond.eff.0     = est.cond.eff.0,
    #                               est.cond.eff.1     = est.cond.eff.1,
    #                               adj.mod.loc        = parms$adj.mod.loc, 
    #                               adj.mthd           = NULL,
    #                               adj.form.true      = NULL,
    #                               type.var           = NULL,
    #                               whl.est.cond.eff.0 = NULL,
    #                               whl.est.cond.eff.1 = NULL)      # different from IPW, the average will not be weighted to large numbers
    # }
    
    ord.ux <- order(trt.eff.g)
    goodness  <- rep(0, length(ux) - 1)
    
    if (parms$adj.mod.loc == "split") {
      
      for (i in 1:(length(ux) - 1)) {
        
        ind.l <- x %in% ux[ord.ux[1:i]]
        ind.r <- x %in% ux[ord.ux[(i+1):length(ux)]]
        
        n.l <- sum(ind.l)
        n.r <- sum(ind.r)
        
        # if ((n.l == 1) | (n.r == 1)) {
        #   next
        # }
        
        data.node.l <- data.node[ind.l, ]
        data.node.r <- data.node[ind.r, ]
        
        # Skip this split if there is only treated/untreated unit in one of the levels
        if ((length(unique(data.node.l$A)) == 1) | (length(unique(data.node.r$A)) == 1)){
          next
        }
        
        tmp.l <- gen.fullrank.g(df            = data.node.l,
                                adj.form.true = parms$form.true)
        tmp.r <- gen.fullrank.g(df            = data.node.r,
                                adj.form.true = parms$form.true)
        
        tmp.l <- withWarnings(est.cond.eff(df        = tmp.l$df.fullrank,
                                           method    = parms$adj.mthd,
                                           form.true = tmp.l$adj.form.true.updated,
                                           type.var  = parms$type.var))
        tmp.r <- withWarnings(est.cond.eff(df        = tmp.r$df.fullrank,
                                           method    = parms$adj.mthd,
                                           form.true = tmp.r$adj.form.true.updated,
                                           type.var  = parms$type.var))
        var.rb.l <- tmp.l$value$var.rb
        var.rb.r <- tmp.r$value$var.rb
        w.l <- tmp.l$value$w
        w.r <- tmp.r$value$w
        
        # when the outcome model fitting produces warning
        # use outside model fitting results for conditional means
        if (!is.null(tmp.l$warnings)) {
          
          # when there is only rank-deficient warning, use the local prediction
          cond.warnings <- T
          for (warnings.i in 1:length(tmp.l$warnings)) {
            cond.warnings <- cond.warnings & (tmp.l$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          if (cond.warnings) {
            est.cond.eff.1.l <- tmp.l$value$pred.A.1
            est.cond.eff.0.l <- tmp.l$value$pred.A.0
          } else {
            est.cond.eff.0.l <- whl.est.cond.eff.0[ind.l]
            est.cond.eff.1.l <- whl.est.cond.eff.1[ind.l]
          }
          
        } else {
          est.cond.eff.0.l <- tmp.l$value$pred.A.0
          est.cond.eff.1.l <- tmp.l$value$pred.A.1
        }
        
        if (!is.null(tmp.r$warnings)) {
          
          # when there is only rank-deficient warning, use the local prediction
          cond.warnings <- T
          for (warnings.i in 1:length(tmp.r$warnings)) {
            cond.warnings <- cond.warnings & (tmp.r$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          if (cond.warnings) {
            est.cond.eff.0.r <- tmp.r$value$pred.A.0
            est.cond.eff.1.r <- tmp.r$value$pred.A.1
          } else {
            est.cond.eff.0.r <- whl.est.cond.eff.0[ind.r]
            est.cond.eff.1.r <- whl.est.cond.eff.1[ind.r]
          }
          
        } else {
          est.cond.eff.0.r <- tmp.r$value$pred.A.0
          est.cond.eff.1.r <- tmp.r$value$pred.A.1
        }
        
        mu.1l <- mean(est.cond.eff.1.l)
        mu.0l <- mean(est.cond.eff.0.l)
        mu.1r <- mean(est.cond.eff.1.r)
        mu.0r <- mean(est.cond.eff.0.r)
        
        # Need to change if the outcome is not continuous
        # involve the derivative of the linear function on the outcome scale
        if (parms$type.var == "cont") {
          
          g.b.l.1 <- data.frame(w.l) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.r.1 <- data.frame(w.r) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.l.0 <- data.frame(w.l) %>% mutate(A = 0) %>% apply(2, mean)
          g.b.r.0 <- data.frame(w.r) %>% mutate(A = 0) %>% apply(2, mean)
          
        } else if (parms$type.var == "bin") {
          
          w.1l <- w.l %>% mutate(A = 1)
          w.0l <- w.l %>% mutate(A = 0)
          w.1r <- w.r %>% mutate(A = 1)
          w.0r <- w.r %>% mutate(A = 0)
          
          g.b.l.1 <- apply(w.1l * est.cond.eff.1.l * (1 - est.cond.eff.1.l), 2, mean)
          g.b.r.1 <- apply(w.1r * est.cond.eff.1.r * (1 - est.cond.eff.1.r), 2, mean)
          g.b.l.0 <- apply(w.0l * est.cond.eff.0.l * (1 - est.cond.eff.0.l), 2, mean)
          g.b.r.0 <- apply(w.0r * est.cond.eff.0.r * (1 - est.cond.eff.0.r), 2, mean)
        }
        
        left.p.1 <- as.matrix(data.frame(est.cond.eff.1.l))
        left.p.0 <- as.matrix(data.frame(est.cond.eff.0.l))
        right.p.1 <- as.matrix(data.frame(est.cond.eff.1.r))
        right.p.0 <- as.matrix(data.frame(est.cond.eff.0.r))
        
        # Calculate variance estimators
        var.1l <- g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1 / (n.l^2) * sum( (left.p.1 - mu.1l)^2 )
        var.1r <- g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1 / (n.r^2) * sum( (right.p.1 - mu.1r)^2 )
        var.0l <- g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1 / (n.l^2) * sum( (left.p.0 - mu.0l)^2 )
        var.0r <- g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1 / (n.r^2) * sum( (right.p.0 - mu.0r)^2 )
        
        var.1l <- as.numeric(var.1l)
        var.0l <- as.numeric(var.0l)
        var.1r <- as.numeric(var.1r)
        var.0r <- as.numeric(var.0r)
        var    <- var.1l + var.0l + var.1r + var.0r
        
        if (var < 0) {
          next
        }
        
        goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var))^2
        
      }
      
    } else { # if (is.null(est.cond.eff.0)) {} else
      
      for (i in 1:(length(ux) - 1)) {
        
        ind.l <- x %in% ux[ord.ux[1:i]]
        ind.r <- x %in% ux[ord.ux[(i+1):length(ux)]]
        
        n.l <- sum(ind.l)
        n.r <- sum(ind.r)
        
        # if ((n.l == 1) | (n.r == 1)) {
        #   next
        # }
        
        mu.1l <- mean(est.cond.eff.1[ind.l])
        mu.0l <- mean(est.cond.eff.0[ind.l])
        mu.1r <- mean(est.cond.eff.1[ind.r])
        mu.0r <- mean(est.cond.eff.0[ind.r])
        
        # Need to change if the outcome is not continuous
        # involve the derivative of the linear function on the outcome scale
        w.l <- w[ind.l, ]
        w.r <- w[ind.r, ]
        
        if (parms$type.var == "cont") {
          
          g.b.l.1 <- data.frame(w.l) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.r.1 <- data.frame(w.r) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.l.0 <- data.frame(w.l) %>% mutate(A = 0) %>% apply(2, mean)
          g.b.r.0 <- data.frame(w.r) %>% mutate(A = 0) %>% apply(2, mean)
          
        } else if (parms$type.var == "bin") {
          
          w.1l <- w.l %>% mutate(A = 1)
          w.0l <- w.l %>% mutate(A = 0)
          w.1r <- w.r %>% mutate(A = 1)
          w.0r <- w.r %>% mutate(A = 0)
          
          g.b.l.1 <- apply(w.1l * est.cond.eff.1[ind.l] * (1 - est.cond.eff.1[ind.l]), 2, mean)
          g.b.r.1 <- apply(w.1r * est.cond.eff.1[ind.r] * (1 - est.cond.eff.1[ind.r]), 2, mean)
          g.b.l.0 <- apply(w.0l * est.cond.eff.0[ind.l] * (1 - est.cond.eff.0[ind.l]), 2, mean)
          g.b.r.0 <- apply(w.0r * est.cond.eff.0[ind.r] * (1 - est.cond.eff.0[ind.r]), 2, mean)
          
        }
        
        left.p.1 <- as.matrix(data.frame(est.cond.eff.1[ind.l]))
        left.p.0 <- as.matrix(data.frame(est.cond.eff.0[ind.l]))
        right.p.1 <- as.matrix(data.frame(est.cond.eff.1[ind.r]))
        right.p.0 <- as.matrix(data.frame(est.cond.eff.0[ind.r]))
        
        # Calculate variance estimators
        var.1l <- g.b.l.1 %*% var.rb %*% g.b.l.1 + 1 / (n.l^2) * sum( (left.p.1 - mu.1l)^2 )
        var.1r <- g.b.r.1 %*% var.rb %*% g.b.r.1 + 1 / (n.r^2) * sum( (right.p.1 - mu.1r)^2 )
        var.0l <- g.b.l.0 %*% var.rb %*% g.b.l.0 + 1 / (n.l^2) * sum( (left.p.0 - mu.0l)^2 )
        var.0r <- g.b.r.0 %*% var.rb %*% g.b.r.0 + 1 / (n.r^2) * sum( (right.p.0 - mu.0r)^2 )
        
        var.1l <- as.numeric(var.1l)
        var.1r <- as.numeric(var.1r)
        var.0l <- as.numeric(var.0l)
        var.0r <- as.numeric(var.0r)
        var    <- var.1l + var.0l + var.1r + var.0r
        
        if (var < 0) {
          next
        }
        
        goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var))^2
        
      }
      
    } # if (is.null(est.cond.eff.0)) else
    
    
    # Direction is the ordered categories
    direction <- ux[ord.ux] 
    
  }
  
  list(goodness  = goodness,
       direction = direction)
  
}

etemp.g <- function(y, wt, parms) {
  
  n <- length(y)
  
  # Finding observations into each node
  sub.ind   <- y
  sub.x     <- parms$covariates[sub.ind, ]
  A         <- parms$trt[sub.ind]
  Y         <- parms$response[sub.ind]
  data.node <- data.frame(Y, A, sub.x)
  
  whl.est.cond.eff.1 <- parms$whl.est.cond.eff.1[sub.ind]
  whl.est.cond.eff.0 <- parms$whl.est.cond.eff.0[sub.ind]
  
  # Weird how this does not work when rss is used.
  # code still runs without problem
  wmean <- sum(y * wt) / sum(wt)
  rss <- sum(wt * (y - wmean)^2)
  
  if (parms$adj.mod.loc != "out"){
    
    if (length(unique(data.node$A)) == 1) {
      
      est.cond.eff.0 <- whl.est.cond.eff.0
      est.cond.eff.1 <- whl.est.cond.eff.1
      
    } else {
      tmp <- gen.fullrank.g(df            = data.node,
                            adj.form.true = parms$form.true)
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
    
  } else { # "out"
    est.cond.eff.0 <- whl.est.cond.eff.0
    est.cond.eff.1 <- whl.est.cond.eff.1
  }
  
  mu.1 <- mean(est.cond.eff.1)
  mu.0 <- mean(est.cond.eff.0)
  
  # Average treatment effect
  avg.trt.effct <- mu.1 - mu.0
  
  list(label = avg.trt.effct, deviance = rss)
}
