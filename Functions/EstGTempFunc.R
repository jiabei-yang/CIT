calc.g.trt.eff <- function(data.node, x, ux, est.cond.eff.0, est.cond.eff.1, adj.mthd, adj.form.true, type.var = "cont"){
  
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
    
    if (is.null(est.cond.eff.0)) { # when the conditional means are estimated within split
      
      tmp <- gen.fullrank.g(df            = data.sub.node,
                            adj.form.true = adj.form.true)
      
      tmp <- try(est.cond.eff(df        = tmp$df.fullrank,
                              method    = adj.mthd,
                              form.true = tmp$adj.form.true.updated, 
                              type.var  = type.var))
      
      # If all the covariates are deleted in df.fullrank, will include NA in the estimated trt.eff
      if ("try-error" %in% class(tmp)) {
        trt.eff <- c(trt.eff, NA)
        next
      }
      
      est.cond.eff.1.sub <- tmp$pred.A.1
      est.cond.eff.0.sub <- tmp$pred.A.0
      
    }
    
    trt.eff     <- c(trt.eff, 
                     mean(est.cond.eff.1.sub) - mean(est.cond.eff.0.sub))
    
    
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
  
  est.cond.eff.1 <- parms$est.cond.eff.1[sub.ind]
  est.cond.eff.0 <- parms$est.cond.eff.0[sub.ind]
  w              <- parms$w[sub.ind, ]
  
  var.rb       <- parms$var.rb
  adj.mthd     <- parms$adj.mthd
  form.true    <- parms$form.true
  mod.insplt   <- parms$mod.insplt
  
  type.var     <- parms$type.var
  num.truc.obs <- parms$num.truc.obs
  
  if (!is.null(mod.insplt)) {
    if (!mod.insplt) {
      
      tmp <- gen.fullrank.g(df            = data.node,
                            adj.form.true = form.true)
      
      tmp <- est.cond.eff(df        = tmp$df.fullrank,
                          method    = adj.mthd,
                          form.true = tmp$adj.form.true.updated,
                          type.var  = type.var)
      
      est.cond.eff.1 <- tmp$pred.A.1
      est.cond.eff.0 <- tmp$pred.A.0
      var.rb         <- tmp$var.rb
      w              <- tmp$w
      
    }
  }

  if (continuous){
    
    # Skip the first 10 and last 10 splits, technically not necessary 
    # as minsplit takes care of that
    # goodness <- NULL
    # direction <- NULL
    goodness <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    if (is.null(est.cond.eff.0)){ # model in split
      
      for (i in (num.truc.obs:(n-num.truc.obs))) {
        
        # if all observations are treated or not treated in the son nodes,
        # skip this loop
        if ((length(unique(data.node$A[1:i])) == 1) | (length(unique(data.node$A[(i+1):n])) == 1)){
          next
        }
        
        tmp.l <- gen.fullrank.g(df            = data.node[1:i, ],
                                adj.form.true = form.true)
        tmp.r <- gen.fullrank.g(df            = data.node[(i+1):n, ],
                                adj.form.true = form.true)
        
        tmp.l <- est.cond.eff(df        = tmp.l$df.fullrank,
                              method    = adj.mthd,
                              form.true = tmp.l$adj.form.true.updated,
                              type.var  = type.var)
        tmp.r <- est.cond.eff(df        = tmp.r$df.fullrank,
                              method    = adj.mthd,
                              form.true = tmp.r$adj.form.true.updated,
                              type.var  = type.var)
        
        est.cond.eff.1.l <- tmp.l$pred.A.1
        est.cond.eff.0.l <- tmp.l$pred.A.0
        est.cond.eff.1.r <- tmp.r$pred.A.1
        est.cond.eff.0.r <- tmp.r$pred.A.0
        
        var.rb.l <- tmp.l$var.rb
        var.rb.r <- tmp.r$var.rb
        
        w.l <- tmp.l$w
        w.r <- tmp.r$w
        
        mu.1l <- mean(est.cond.eff.1.l)
        mu.0l <- mean(est.cond.eff.0.l)
        mu.1r <- mean(est.cond.eff.1.r)
        mu.0r <- mean(est.cond.eff.0.r)
        
        if (type.var == "cont") {
          
          g.b.l.1 <- data.frame(w.l) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.r.1 <- data.frame(w.r) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.l.0 <- data.frame(w.l) %>% mutate(A = 0) %>% apply(2, mean)
          g.b.r.0 <- data.frame(w.r) %>% mutate(A = 0) %>% apply(2, mean)
          
        } else if (type.var == "bin") {
          
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
      
      for (i in (num.truc.obs:(n-num.truc.obs))) {
        
        mu.1l <- mean(est.cond.eff.1[1:i])
        mu.0l <- mean(est.cond.eff.0[1:i])
        mu.1r <- mean(est.cond.eff.1[(i+1):n])
        mu.0r <- mean(est.cond.eff.0[(i+1):n])
        
        # Need to change if the outcome is not continuous
        # involve the derivative of the linear function on the outcome scale
        if (type.var == "cont") {
          
          g.b.l.1 <- data.frame(w[1:i, ]) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.r.1 <- data.frame(w[(i+1):n, ]) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.l.0 <- data.frame(w[1:i, ]) %>% mutate(A = 0) %>% apply(2, mean)
          g.b.r.0 <- data.frame(w[(i+1):n, ]) %>% mutate(A = 0) %>% apply(2, mean)
          
        } else if (type.var == "bin") {
          
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
    
    if (is.null(est.cond.eff.0)) {
      
      trt.eff.g <- calc.g.trt.eff(data.node, x, ux,
                                  est.cond.eff.0 = NULL,
                                  est.cond.eff.1 = NULL,
                                  adj.mthd       = adj.mthd,
                                  adj.form.true  = form.true, 
                                  type.var       = type.var)
      
    } else { # when the conditional mean model is fitted within split, 
             # need to fit conditional mean models within each level
      
      trt.eff.g <- calc.g.trt.eff(data.node, x, ux,
                                  est.cond.eff.0 = est.cond.eff.0,
                                  est.cond.eff.1 = est.cond.eff.1,
                                  adj.mthd       = NULL,
                                  adj.form.true  = NULL,
                                  type.var       = NULL)
    }
    
    ord.ux <- order(trt.eff.g)
    goodness  <- rep(0, length(ux) - 1)
    
    if (is.null(est.cond.eff.0)) {
      
      for (i in 1:(length(ux) - 1)) {
        
        ind.l <- x %in% ux[ord.ux[1:i]]
        ind.r <- x %in% ux[ord.ux[(i+1):length(ux)]]
        
        n.l <- sum(ind.l)
        n.r <- sum(ind.r)
        
        if ((n.l == 1) | (n.r == 1)) {
          next
        }
        
        data.node.l <- data.node[ind.l, ]
        data.node.r <- data.node[ind.r, ]
        
        if ((length(unique(data.node.l$A)) == 1) | (length(unique(data.node.r$A)) == 1)){
          next
        }
        
        tmp.l <- gen.fullrank.g(df            = data.node.l,
                                adj.form.true = form.true)
        tmp.r <- gen.fullrank.g(df            = data.node.r,
                                adj.form.true = form.true)
        
        tmp.l <- est.cond.eff(df        = tmp.l$df.fullrank,
                              method    = adj.mthd,
                              form.true = tmp.l$adj.form.true.updated,
                              type.var  = type.var)
        tmp.r <- est.cond.eff(df        = tmp.r$df.fullrank,
                              method    = adj.mthd,
                              form.true = tmp.r$adj.form.true.updated,
                              type.var  = type.var)
        
        est.cond.eff.1.l <- tmp.l$pred.A.1
        est.cond.eff.0.l <- tmp.l$pred.A.0
        est.cond.eff.1.r <- tmp.r$pred.A.1
        est.cond.eff.0.r <- tmp.r$pred.A.0
        
        var.rb.l <- tmp.l$var.rb
        var.rb.r <- tmp.r$var.rb
        
        w.l <- tmp.l$w
        w.r <- tmp.r$w
        
        mu.1l <- mean(est.cond.eff.1.l)
        mu.0l <- mean(est.cond.eff.0.l)
        mu.1r <- mean(est.cond.eff.1.r)
        mu.0r <- mean(est.cond.eff.0.r)
        
        # Need to change if the outcome is not continuous
        # involve the derivative of the linear function on the outcome scale
        if (type.var == "cont") {
          
          g.b.l.1 <- data.frame(w.l) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.r.1 <- data.frame(w.r) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.l.0 <- data.frame(w.l) %>% mutate(A = 0) %>% apply(2, mean)
          g.b.r.0 <- data.frame(w.r) %>% mutate(A = 0) %>% apply(2, mean)
          
        } else if (type.var == "bin") {
          
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
        
        goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
        
      }
      
    } else { # if (is.null(est.cond.eff.0)) {} else
      
      for (i in 1:(length(ux) - 1)){
        
        ind.l <- x %in% ux[ord.ux[1:i]]
        ind.r <- x %in% ux[ord.ux[(i+1):length(ux)]]
        
        n.l <- sum(ind.l)
        n.r <- sum(ind.r)
        
        if ((n.l == 1) | (n.r == 1)) {
          next
        }

        mu.1l <- mean(est.cond.eff.1[ind.l])
        mu.0l <- mean(est.cond.eff.0[ind.l])
        mu.1r <- mean(est.cond.eff.1[ind.r])
        mu.0r <- mean(est.cond.eff.0[ind.r])
        
        # Need to change if the outcome is not continuous
        # involve the derivative of the linear function on the outcome scale
        w.l <- w[ind.l, ]
        w.r <- w[ind.r, ]
        
        if (type.var == "cont") {
          
          g.b.l.1 <- data.frame(w.l) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.r.1 <- data.frame(w.r) %>% mutate(A = 1) %>% apply(2, mean)
          g.b.l.0 <- data.frame(w.l) %>% mutate(A = 0) %>% apply(2, mean)
          g.b.r.0 <- data.frame(w.r) %>% mutate(A = 0) %>% apply(2, mean)
          
        } else if (type.var == "bin") {
          
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
        
        goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var.1l + var.0l + var.1r + var.0r))^2
        
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
  
  est.cond.eff.1 <- parms$est.cond.eff.1[sub.ind]
  est.cond.eff.0 <- parms$est.cond.eff.0[sub.ind]
  adj.mthd  <- parms$adj.mthd
  form.true <- parms$form.true
  type.var  <- parms$type.var

  if (is.null(est.cond.eff.0)){
    
    tmp <- gen.fullrank.g(df            = data.node,
                          adj.form.true = form.true)
    tmp <- est.cond.eff(df        = tmp$df.fullrank,
                        method    = adj.mthd,
                        form.true = tmp$adj.form.true.updated,
                        type.var  = type.var)
    
    est.cond.eff.1 = tmp$pred.A.1
    est.cond.eff.0 = tmp$pred.A.0
    
  }
  
  mu.1 <- mean(est.cond.eff.1)
  mu.0 <- mean(est.cond.eff.0)
  
  # Average treatment effect
  avg.trt.effct <- mu.1 - mu.0
  
  # Weird how this does not work when rss is used.
  # code still runs without problem
  wmean <- sum(y * wt) / sum(wt)
  rss <- sum(wt * (y - wmean)^2)
  
  list(label = avg.trt.effct, deviance = rss)
}
