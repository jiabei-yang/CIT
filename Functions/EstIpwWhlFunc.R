stemp.ipw.whl <- function(y, wt, x, parms, continuous){
  
  # y:          vector of response values at the node;
  # wt:         vector of weights for observations at the node;
  # x:          vector of the x values at the node;
  # parms:      trt;
  #             covariates;
  #             response;
  #             ..... response.type: continuous, binary, categorical(?);
  #             ..... family of binary?
  # continuous: logical. Indicator for covariate type of x.
  
  n     <- length(y)
  n.whl <- nrow(parms$covariates)
  
  # sub.ind <- match(y, parms$response)
  sub.ind <- y
  sub.x   <- parms$covariates[sub.ind, ]
  A       <- parms$trt[sub.ind]
  Y       <- parms$response[sub.ind]
  
  data.node     <- data.frame(Y, A, sub.x)
  data.node.noy <- data.node[, -1]
  
  mod.insplt    <- parms$mod.insplt 
  num.truc.obs  <- parms$num.truc.obs
  
  if (!is.null(parms$prop.sc)){                      
    
    # Model outside node or not, assign value to prop.sc, w 
    prop.sc   <- parms$prop.sc[sub.ind]
    w <- parms$w[sub.ind, ]
    
  } else {
    
    if (!mod.insplt){
      
      # print(data.node.noy)
      tmp <- gen.fullrank.ipw(df.noy           = data.node.noy, 
                              propsc.form.true = parms$form.true)
      
      tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                         method    = parms$propsc.mthd,
                         form.true = tmp$propsc.form.true.updated)
      prop.sc <- tmp$prop.sc
      w       <- tmp$w
      
    } # end if (!mod.insplt){
    
  } # end if (!is.null(parms$prop.sc)) else
  
  if (continuous){
    
    goodness <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    if (exists("prop.sc")){
      # if (is.null(mod.insplt) | !mod.insplt){
      
      # Ensure the model is at least fitted on 30 obs
      for (i in (num.truc.obs:(n-num.truc.obs))){
        
        mu.1l <- mean((Y * A/ prop.sc)[1:i])
        mu.0l <- mean((Y * (1 - A)/ (1 - prop.sc))[1:i])
        mu.1r <- mean((Y * A/ prop.sc)[(i+1):n])
        mu.0r <- mean((Y * (1 - A)/ (1 - prop.sc))[(i+1):n])
        
        w.l <- w[1:i, ] 
        w.r <- w[(i+1):n, ]
        
        h.l <- apply((A * Y * (1 - prop.sc) / prop.sc + (1 - A) * Y * prop.sc / (1 - prop.sc))[1:i] * w.l,
                     2,
                     mean)
        h.r <- apply((A * Y * (1 - prop.sc) / prop.sc + (1 - A) * Y * prop.sc / (1 - prop.sc))[(i+1):n] * w.r,
                     2,
                     mean)
        
        e.bb <- as.matrix(t(prop.sc * w)) %*% as.matrix((1 - prop.sc) * w) / n.whl
        
        I.i <- (((x %in% x[1:i]) * Y * A) / (i / n.whl * prop.sc) - ((x %in% x[1:i]) * Y * (1 - A)) / (i / n.whl * (1 - prop.sc))) -
          (((x %in% x[(i+1):n]) * Y * A) / ((n - i) / n.whl * prop.sc) - ((x %in% x[(i+1):n]) * Y * (1 - A)) / ((n - i) / n.whl * (1 - prop.sc))) - 
          ((mu.1l - mu.0l) - (mu.1r - mu.0r)) - 
          as.numeric((A - prop.sc) * (h.l - h.r) %*% solve(e.bb) %*% t(w))
        var <- (mean(I.i^2) - ((n-i) / n.whl * (mu.1l - mu.0l) + i / n.whl * (mu.1r - mu.0r))^2 / (i * (n-i) / n.whl^2)) / n
        
        # If the variance is less than 0, will assume the goodness for the split is 0
        if (var < 0) {
          next
        } else {
          goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var))^2
        }
        
        direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))
      }
      
    } else { # if (exists("prop.sc")) else
      
      for (i in (num.truc.obs:(n-num.truc.obs))){
        
        data.node.noy.l <- data.node.noy[1:i, ]
        data.node.noy.r <- data.node.noy[(i+1):n, ]
        
        
        tmp.l <- gen.fullrank.ipw(df.noy           = data.node.noy.l, 
                                  propsc.form.true = parms$form.true)
        tmp.r <- gen.fullrank.ipw(df.noy           = data.node.noy.r, 
                                  propsc.form.true = parms$form.true)
        
        tmp.l <- est.prop.sc(df.noy    = tmp.l$df.noy.fullrank,
                             method    = parms$propsc.mthd,
                             form.true = tmp.l$propsc.form.true.updated)
        tmp.r <- est.prop.sc(df.noy    = tmp.r$df.noy.fullrank,
                             method    = parms$propsc.mthd,
                             form.true = tmp.r$propsc.form.true.updated)
        
        prop.sc.l <- tmp.l$prop.sc
        prop.sc.r <- tmp.r$prop.sc
        
        w.l <- tmp.l$w
        w.r <- tmp.r$w
        
        mu.1l <- mean((Y * A)[1:i] / prop.sc.l)
        mu.0l <- mean((Y * (1 - A))[1:i] / (1 - prop.sc.l))
        mu.1r <- mean((Y * A)[(i+1):n] / prop.sc.r)
        mu.0r <- mean((Y * (1 - A))[(i+1):n] / (1 - prop.sc.r))
        
        # w.l <- w[1:i, ] 
        # w.r <- w[(i+1):n, ]
        
        h.l <- apply((A * Y)[1:i] * (1 - prop.sc.l) / prop.sc.l + ((1 - A) * Y)[1:i] * prop.sc.l / (1 - prop.sc.l) * w.l,
                     2,
                     mean)
        h.r <- apply((A * Y)[(i+1):n] * (1 - prop.sc.r) / prop.sc.r + ((1 - A) * Y)[(i+1):n] * prop.sc.r / (1 - prop.sc.r) * w.r,
                     2,
                     mean)
        
        prop.sc <- c(prop.sc.l, prop.sc.r)
        e.bb.l <- as.matrix(t(prop.sc.l * w.l)) %*% as.matrix((1 - prop.sc.l) * w.l) / n
        e.bb.r <- as.matrix(t(prop.sc.r * w.r)) %*% as.matrix((1 - prop.sc.r) * w.r) / n
        
        I.i <- (((x %in% x[1:i]) * Y * A) / (i / n * prop.sc) - ((x %in% x[1:i]) * Y * (1 - A)) / (i / n * (1 - prop.sc))) -
          (((x %in% x[(i+1):n]) * Y * A) / ((n - i) / n * prop.sc) - ((x %in% x[(i+1):n]) * Y * (1 - A)) / ((n - i) / n * (1 - prop.sc))) - 
          ((mu.1l - mu.0l) - (mu.1r - mu.0r)) - 
          c(as.numeric((A[1:i] - prop.sc.l) * h.l %*% solve(e.bb.l) %*% t(w.l)), rep(0, n-i)) -
          c(rep(0, i), as.numeric((A[(i+1):n] - prop.sc.r) * h.r %*% solve(e.bb.r) %*% t(w.r)))
        
        # if ("try-error" %in% class(I.i)) {
        #   next
        # }
        
        var <- (mean(I.i^2) - ((n-i) / n * (mu.1l - mu.0l) + i / n * (mu.1r - mu.0r))^2 / (i * (n-i) / n^2)) / n
        
        # If the variance is less than 0, will assume the goodness for the split is 0
        if (var < 0) {
          next
        } else {
          goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var))^2
        }
        
        direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))
      }
    }
    
  } else {
    
    ux <- sort(unique(x))
    
    # Order the levels of the categorical covariates by their treatment effect
    if (exists("prop.sc")) {
      trt.eff.ipw <- calc.ipw.trt.eff(data.node, x, ux, 
                                      prop.sc,
                                      propsc.mthd      = parms$propsc.mthd,
                                      propsc.form.true = parms$form.true)
    } else { # when the propensity score is fitted within split, 
      # need to fit propensity score models within each level
      trt.eff.ipw <- calc.ipw.trt.eff(data.node, x, ux, 
                                      prop.sc          = NULL, 
                                      propsc.mthd      = parms$propsc.mthd,
                                      propsc.form.true = parms$form.true)
    }
    
    ord.ux <- order(trt.eff.ipw)
    goodness  <- rep(0, length(ux) - 1)
    
    if (exists("prop.sc")){
      
      for (i in 1:(length(ux) - 1)){
        
        ind.l <- x %in% ux[ord.ux[1:i]]
        ind.r <- x %in% ux[ord.ux[(i+1):length(ux)]]
        
        mu.1l <- mean((Y * A/ prop.sc)[ind.l])
        mu.0l <- mean((Y * (1 - A)/ (1 - prop.sc))[ind.l])
        mu.1r <- mean((Y * A/ prop.sc)[ind.r])
        mu.0r <- mean((Y * (1 - A)/ (1 - prop.sc))[ind.r])
        
        w.l <- w[ind.l, ] 
        w.r <- w[ind.r, ]
        
        if (sum(ind.l) == 1) {
          h.l <- apply(t((A * Y * (1 - prop.sc) / prop.sc + (1 - A) * Y * prop.sc / (1 - prop.sc))[ind.l] * w.l),
                       2,
                       mean)
        } else {
          h.l <- apply((A * Y * (1 - prop.sc) / prop.sc + (1 - A) * Y * prop.sc / (1 - prop.sc))[ind.l] * w.l,
                       2,
                       mean)
        }
        
        if (sum(ind.r) == 1) {
          h.r <- apply(t((A * Y * (1 - prop.sc) / prop.sc + (1 - A) * Y * prop.sc / (1 - prop.sc))[ind.r] * w.r),
                       2,
                       mean)
        } else {
          h.r <- apply((A * Y * (1 - prop.sc) / prop.sc + (1 - A) * Y * prop.sc / (1 - prop.sc))[ind.r] * w.r,
                       2,
                       mean)
        }
        
        e.bb <- as.matrix(t(prop.sc * w)) %*% as.matrix((1 - prop.sc) * w) / n
        
        n.l <- sum(ind.l)
        n.r <- sum(ind.r)
        
        I.i <- ((ind.l * Y * A) / (n.l / n * prop.sc) - (ind.l * Y * (1 - A)) / (n.l / n * (1 - prop.sc))) -
          ((ind.r * Y * A) / (n.r / n * prop.sc) - (ind.r * Y * (1 - A)) / (n.r / n * (1 - prop.sc))) - 
          ((mu.1l - mu.0l) - (mu.1r - mu.0r)) - 
          as.numeric((A - prop.sc) * (h.l - h.r) %*% solve(e.bb) %*% t(w))
        var <- (mean(I.i^2) - (n.r / n * (mu.1l - mu.0l) + n.l / n * (mu.1r - mu.0r))^2 / (n.l * n.r / n^2)) / n
        
        # If the variance is less than 0, will assume the goodness for the split is 0
        if (var < 0) {
          next
        } else {
          goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var))^2
        }
        
      }
      
      
    } else { # if (exists("prop.sc")){ else
      
      prop.sc  <- rep(0, n)
      var.term <- rep(0, n) 
      
      for (i in 1:(length(ux) - 1)){
        
        ind.l <- x %in% ux[ord.ux[1:i]]
        ind.r <- x %in% ux[ord.ux[(i+1):length(ux)]]
        
        n.l <- sum(ind.l)
        n.r <- sum(ind.r)
        
        if ((n.l == 1) | (n.r == 1)) {
          next
        }
        
        data.node.noy.l <- data.node.noy[ind.l, ]
        data.node.noy.r <- data.node.noy[ind.r, ]
        
        tmp.l <- gen.fullrank.ipw(df.noy           = data.node.noy.l, 
                                  propsc.form.true = parms$form.true)
        
        # if ("try-error" %in% class(tmp.l)) {
        #   next
        # }
        
        tmp.r <- gen.fullrank.ipw(df.noy           = data.node.noy.r, 
                                  propsc.form.true = parms$form.true)
        # if ("try-error" %in% class(tmp.r)) {
        #   next
        # }
        
        tmp.l <- est.prop.sc(df.noy    = tmp.l$df.noy.fullrank,
                             method    = parms$propsc.mthd,
                             form.true = tmp.l$propsc.form.true.updated)
        tmp.r <- est.prop.sc(df.noy    = tmp.r$df.noy.fullrank,
                             method    = parms$propsc.mthd,
                             form.true = tmp.r$propsc.form.true.updated)
        
        prop.sc.l <- tmp.l$prop.sc
        prop.sc.r <- tmp.r$prop.sc
        w.l <- tmp.l$w
        w.r <- tmp.r$w
        
        mu.1l <- mean((Y * A)[ind.l] / prop.sc.l)
        mu.0l <- mean((Y * (1 - A))[ind.l] / (1 - prop.sc.l))
        mu.1r <- mean((Y * A)[ind.r] / prop.sc.r)
        mu.0r <- mean((Y * (1 - A))[ind.r]/ (1 - prop.sc.r))
        
        h.l <- apply(((A * Y)[ind.l] * (1 - prop.sc.l) / prop.sc.l + ((1 - A) * Y)[ind.l] * prop.sc.l / (1 - prop.sc.l)) * w.l,
                     2,
                     mean)
        h.r <- apply(((A * Y)[ind.r] * (1 - prop.sc.r) / prop.sc.r + ((1 - A) * Y)[ind.r] * prop.sc.r / (1 - prop.sc.r)) * w.r,
                     2,
                     mean)
        
        e.bb.l <- as.matrix(t(prop.sc.l * w.l)) %*% as.matrix((1 - prop.sc.l) * w.l) / n
        e.bb.r <- as.matrix(t(prop.sc.r * w.r)) %*% as.matrix((1 - prop.sc.r) * w.r) / n
        
        prop.sc[ind.l] <- prop.sc.l
        prop.sc[ind.r] <- prop.sc.r
        
        term.l <- as.numeric((A[ind.l] - prop.sc.l) * h.l %*% solve(e.bb.l) %*% t(w.l))
        term.r <- as.numeric((A[ind.r] - prop.sc.r) * h.r %*% solve(e.bb.r) %*% t(w.r))
        
        # if (("try-error" %in% class(term.l)) | ("try-error" %in% class(term.r))) {
        #   next
        # }
        
        var.term[ind.l] <- term.l
        var.term[ind.r] <- term.r
        
        I.i <- ((ind.l * Y * A) / (n.l / n * prop.sc) - (ind.l * Y * (1 - A)) / (n.l / n * (1 - prop.sc))) -
          ((ind.r * Y * A) / (n.r / n * prop.sc) - (ind.r * Y * (1 - A)) / (n.r / n * (1 - prop.sc))) - 
          ((mu.1l - mu.0l) - (mu.1r - mu.0r)) - var.term
        var <- (mean(I.i^2) - (n.r / n * (mu.1l - mu.0l) + n.l / n * (mu.1r - mu.0r))^2 / (n.l * n.r / n^2)) / n
        
        # If the variance is less than 0, will assume the goodness for the split is 0
        if (var < 0) {
          next
        } else {
          goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var))^2
        }
        
      }
      
    }
    
    direction <- ux[ord.ux]
    
  }
  
  list(goodness  = goodness,
       direction = direction)
  
}

etemp.ipw.whl <- function(y, wt, parms) {
  
  n <- length(y)
  
  # sub.ind <- match(y, parms$response)
  sub.ind <- y
  sub.x   <- parms$covariates[sub.ind, ]
  A       <- parms$trt[sub.ind]
  Y       <- parms$response[sub.ind]
  
  data.node = data.frame(Y, A, sub.x)
  data.node.noy <- data.node[, -1]
  
  # Need to move this part of code here for return value earlier
  # TO DO: Weird how this does not work when rss is used.
  # not crucial to the performance 
  wmean <- sum(Y*wt)/sum(wt)
  rss <- sum(wt*(Y-wmean)^2)
  
  # # If there is only one observation in this node, avg.trt.eff cannot be calculated so return earlier
  # # Will cause error in mse calculation if use NA as label
  # if (n == 1) {
  #   return(list(label = NA, deviance = rss))
  # }
  
  if (!is.null(parms$prop.sc)){                    # Model outside node or not
    prop.sc   <- parms$prop.sc[sub.ind]
  } else {
    
    tmp <- gen.fullrank.ipw(df.noy           = data.node.noy, 
                            propsc.form.true = parms$form.true)
    
    tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                       method    = parms$propsc.mthd,
                       form.true = tmp$propsc.form.true.updated)
    prop.sc <- tmp$prop.sc
    
  }
  
  mu.1 <- mean(Y * A/ prop.sc)
  mu.0 <- mean(Y * (1 - A)/ (1 - prop.sc))
  avg.trt.effct <- mu.1 - mu.0
  
  list(label = avg.trt.effct, deviance = rss)
}
