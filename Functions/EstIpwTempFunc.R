calc.ipw.trt.eff <- function(data.node, x, ux, prop.sc, propsc.mod.loc, propsc.mthd, propsc.form.true, whl.propsc, avg.trt.effct){
  
  # calculates ipw mean response difference between treatment and control for each level of x 
  # when x is categorical and the outcome is continuous.
  # data.node: a data frame containing observations at the node.
  # x:         vector of the \texttt{x} values at the node.
  # ux:        unique levels of the categorical covariate x at the node.
  
  trt.eff <- NULL
  
  # Crude mean response difference between treatment and control for each level of x
  for (i in 1:length(ux)){
    
    data.sub.node <- data.node[x == ux[i], ]
    
    # when there are both treated and untreated unit in split, use avg.trt.effct
    if (length(unique(data.sub.node$A)) > 1) {
      
      if (propsc.mod.loc != "split") {
        
        prop.sc.sub <- prop.sc[x == ux[i]]
        mu.sub <- mean(data.sub.node$Y * data.sub.node$A / prop.sc.sub) - 
          mean(data.sub.node$Y * (1 - data.sub.node$A) / (1 - prop.sc.sub))
        
      } else { # when the propensity scores are estimated within split
        
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
        
        mu.sub      <- mean(data.sub.node$Y * data.sub.node$A / prop.sc.sub) - 
          mean(data.sub.node$Y * (1 - data.sub.node$A) / (1 - prop.sc.sub))
        
      }
    } else { # if (length(unique(data.sub.node$A)) > 1) {
      # when there is only treated or untreated unit in split, use avg.trt.effct
      # since the model will be misleading
      mu.sub <- avg.trt.effct
    }
    
    trt.eff     <- c(trt.eff, mu.sub)
    
  } # for loop
  
  return(trt.eff)
  
}



gen.fullrank.ipw <- function(df.noy, propsc.form.true) {
  
  # remove the factor column if the column has only one level
  df.noy.fullrank <- data.frame(A = df.noy$A,
                                data.frame(df.noy[, colnames(df.noy) != "A"])[, sapply(df.noy[, colnames(df.noy) != "A"], 
                                                                                       function(col) length(unique(col))) > 1])
  
  # re-factor the factor columns if there is some levels in the factor column missing
  # drop the missing factor
  df.noy.fullrank[] <- lapply(df.noy.fullrank, 
                              function(x) if(is.factor(x)) factor(x) else x)
  
  # print(df.noy)
  
  if (ncol(df.noy.fullrank) < ncol(df.noy) & (!is.null(propsc.form.true))) {
    propsc.terms <- gsub("A ~ ", "", propsc.form.true)
    propsc.terms <- base::strsplit(propsc.terms, " + ", fixed = T)[[1]]
    
    # Need to take care of the other functional forms in the terms
    # propsc.terms.colname <- gsub("exp[(]", "", propsc.terms)
    # propsc.terms.colname <- gsub("[)]", "", propsc.terms.colname)
    
    # find the terms still in df.noy.fullrank
    propsc.terms <- propsc.terms[rowSums(sapply(names(df.noy.fullrank), function(covariate) grepl(covariate, propsc.terms))) >= 
                                   rowSums(sapply(names(df.noy), function(covariate) grepl(covariate, propsc.terms)))]
    
    # when there is nothing left fit an intercept only model
    if (identical(propsc.terms, character(0))) {
      propsc.form.true.updated <- "A ~ 1"
    } else {
      propsc.form.true.updated <- paste0("A ~ ", paste(propsc.terms, collapse= " + "))
    }
    
  } else { # if (ncol(df.noy.fullrank) < ncol(df.noy) & (!is.null(propsc.form.true)))
    # do we need to do something here? ############
    propsc.form.true.updated <- propsc.form.true
  }
  
  return(list(df.noy.fullrank          = df.noy.fullrank,
              propsc.form.true.updated = propsc.form.true.updated))
}



stemp.ipw <- function(y, wt, x, parms, continuous){
  
  # y:          vector of response values at the node;
  # wt:         vector of weights for observations at the node;
  # x:          vector of the x values at the node;
  # parms:      trt;
  #             covariates;
  #             response;
  #             ..... response.type: continuous, binary, categorical(?);
  #             ..... family of binary?
  # continuous: logical. Indicator for covariate type of x.
  
  n <- length(y) 
  
  # sub.ind <- match(y, parms$response)
  sub.ind <- y
  sub.x   <- parms$covariates[sub.ind, ]
  A       <- parms$trt[sub.ind]
  Y       <- parms$response[sub.ind]
  
  data.node     <- data.frame(Y, A, sub.x)
  data.node.noy <- data.node[, -1]
  
  # For categorical covariate, insplit treatment effect calculation
  whl.propsc <- parms$whl.propsc[sub.ind]
  
  # created for categorical splits, sort categories
  prop.sc <- NULL
  
  # Both outside and in node need w, but different ways to calculate propensity scores
  if (parms$propsc.mod.loc != "split") {
    
    # When the propensity score is fitted outside and there are categorical covariates
    # need to take care of the singular inverse
    tmp <- gen.fullrank.ipw(df.noy           = data.node.noy, 
                            propsc.form.true = parms$form.true)
    tmp <- withWarnings(est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                                    method    = parms$propsc.mthd,
                                    form.true = tmp$propsc.form.true.updated))
    # val.w.used will be the same from both following conditions
    w <- tmp$value$w
    
    if (parms$propsc.mod.loc == "out") {
      # Model outside node or not, assign value to prop.sc
      prop.sc <- whl.propsc
      
    } else { # "node"
      # if 1) propensity score model fitting produces warning or 2) only produce thresholded propensity scores
      # use outside model fitting results
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
    
  }
  
  if (continuous){
    
    goodness <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    if (parms$propsc.mod.loc != "split") {
      # if (is.null(mod.insplt) | !mod.insplt){
      
      # Ensure the model is at least fitted on num.truc.obs obs
      for (i in (parms$num.truc.obs:(n-parms$num.truc.obs))){
        
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
        
        e.bb <- as.matrix(t(prop.sc * w)) %*% as.matrix((1 - prop.sc) * w) / n
        
        I.i <- try((((x %in% x[1:i]) * Y * A) / (i / n * prop.sc) - ((x %in% x[1:i]) * Y * (1 - A)) / (i / n * (1 - prop.sc))) -
                     (((x %in% x[(i+1):n]) * Y * A) / ((n - i) / n * prop.sc) - ((x %in% x[(i+1):n]) * Y * (1 - A)) / ((n - i) / n * (1 - prop.sc))) - 
                     ((mu.1l - mu.0l) - (mu.1r - mu.0r)) - 
                     as.numeric((A - prop.sc) * (h.l - h.r) %*% solve(e.bb) %*% t(w)), 
                   silent = T)
        
        # if there is error in solve(e.bb), skip this loop
        # sometimes exp() function makes the matrix unsolvable
        if ("try-error" %in% class(I.i)) {
          next
        }
        
        var <- (mean(I.i^2) - ((n-i) / n * (mu.1l - mu.0l) + i / n * (mu.1r - mu.0r))^2 / (i * (n-i) / n^2)) / n
        
        # If the variance is less than 0, will assume the goodness for the split is 0
        if (var < 0) {
          next
        } else {
          goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var))^2
        }
        
        direction[i] <- c(sign((mu.1l - mu.0l) - (mu.1r - mu.0r)))
      }
      
    } else { # if  (parms$propsc.mod.loc != "split") else, in split
      
      for (i in (parms$num.truc.obs:(n-parms$num.truc.obs))){
        
        data.node.noy.l <- data.node.noy[1:i, ]
        data.node.noy.r <- data.node.noy[(i+1):n, ]
        
        # if all observations are treated or not treated in the son nodes,
        # skip this loop
        if ((length(unique(data.node.noy.l$A)) == 1) | (length(unique(data.node.noy.r$A)) == 1)){
          next
        }
        
        tmp.l <- gen.fullrank.ipw(df.noy           = data.node.noy.l, 
                                  propsc.form.true = parms$form.true)
        tmp.r <- gen.fullrank.ipw(df.noy           = data.node.noy.r, 
                                  propsc.form.true = parms$form.true)
        
        tmp.l <- withWarnings(est.prop.sc(df.noy    = tmp.l$df.noy.fullrank,
                                          method    = parms$propsc.mthd,
                                          form.true = tmp.l$propsc.form.true.updated))
        tmp.r <- withWarnings(est.prop.sc(df.noy    = tmp.r$df.noy.fullrank,
                                          method    = parms$propsc.mthd,
                                          form.true = tmp.r$propsc.form.true.updated))
        w.l <- tmp.l$value$w
        w.r <- tmp.r$value$w
        
        # if 1) propensity score model fitting produces warning or 2) only produce thresholded propensity scores
        # use outside model fitting results
        if (!is.null(tmp.l$warnings)) {
          
          cond.warnings <- T
          for (warnings.i in 1:length(tmp.l$warnings)) {
            cond.warnings <- cond.warnings & (tmp.l$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          # when there is only rank-deficient warning, use the local prediction
          if (cond.warnings) {
            prop.sc.l <- tmp.l$value$prop.sc
          } else {
            # since whl.propsc already subsetted and sorted
            prop.sc.l <- whl.propsc[1:i]
          }

        } else if (identical(sort(unique(tmp.l$value$prop.sc)), c(0.1, 0.9))) {
          prop.sc.l <- whl.propsc[1:i]
        } else {
          prop.sc.l <- tmp.l$value$prop.sc
        }
        
        if (!is.null(tmp.r$warnings)) {
          
          cond.warnings <- T
          for (warnings.i in 1:length(tmp.r$warnings)) {
            cond.warnings <- cond.warnings & (tmp.r$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          # when there is only rank-deficient warning, use the local prediction
          if (cond.warnings) {
            prop.sc.r <- tmp.r$value$prop.sc
          } else {
            # since whl.propsc already subsetted and sorted
            prop.sc.r <- whl.propsc[(i+1):n]
          }
          
        } else if (identical(sort(unique(tmp.r$value$prop.sc)), c(0.1, 0.9))) { 
          prop.sc.r <- whl.propsc[(i+1):n]
        } else {
          prop.sc.r <- tmp.r$value$prop.sc
        }
        
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
        
        I.i <- try((((x %in% x[1:i]) * Y * A) / (i / n * prop.sc) - ((x %in% x[1:i]) * Y * (1 - A)) / (i / n * (1 - prop.sc))) -
                     (((x %in% x[(i+1):n]) * Y * A) / ((n - i) / n * prop.sc) - ((x %in% x[(i+1):n]) * Y * (1 - A)) / ((n - i) / n * (1 - prop.sc))) - 
                     ((mu.1l - mu.0l) - (mu.1r - mu.0r)) - 
                     c(as.numeric((A[1:i] - prop.sc.l) * h.l %*% solve(e.bb.l) %*% t(w.l)), rep(0, n-i)) -
                     c(rep(0, i), as.numeric((A[(i+1):n] - prop.sc.r) * h.r %*% solve(e.bb.r) %*% t(w.r))), 
                   silent = T)
        
        # if there is error in solve(e.bb), skip this loop
        # sometimes exp() function makes the matrix unsolvable
        if ("try-error" %in% class(I.i)) {
          next
        }
        
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
    
  } else { # categorical covariates
    
    ux <- sort(unique(x))
    
    # Order the levels of the categorical covariates by their treatment effect
    # cannot merge this into ifelse because there is issue when assigning value to prop.sc 
    # prop.sc does not exist for in split
    # if (parms$propsc.mod.loc != "split") {
    #   trt.eff.ipw <- calc.ipw.trt.eff(data.node, x, ux, 
    #                                   prop.sc          = prop.sc,
    #                                   propsc.mod.loc   = parms$propsc.mod.loc,
    #                                   propsc.mthd      = NULL,
    #                                   propsc.form.true = NULL,
    #                                   whl.propsc       = NULL,                      # only needed when in split
    #                                   avg.trt.effct    = parms$avg.trt.effct)
    # } else {
    trt.eff.ipw <- calc.ipw.trt.eff(data.node, x, ux, 
                                    prop.sc          = prop.sc,
                                    propsc.mod.loc   = parms$propsc.mod.loc,
                                    propsc.mthd      = parms$propsc.mthd,
                                    propsc.form.true = parms$form.true,
                                    whl.propsc       = whl.propsc,                      # only needed when in split
                                    avg.trt.effct    = parms$avg.trt.effct)
    # }
    
    ord.ux <- order(trt.eff.ipw)
    goodness  <- rep(0, length(ux) - 1)
    
    if (parms$propsc.mod.loc != "split") {
      
      for (i in 1:(length(ux) - 1)) {
        
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
        
        I.i <- try(((ind.l * Y * A) / (n.l / n * prop.sc) - (ind.l * Y * (1 - A)) / (n.l / n * (1 - prop.sc))) -
                     ((ind.r * Y * A) / (n.r / n * prop.sc) - (ind.r * Y * (1 - A)) / (n.r / n * (1 - prop.sc))) - 
                     ((mu.1l - mu.0l) - (mu.1r - mu.0r)) - 
                     as.numeric((A - prop.sc) * (h.l - h.r) %*% solve(e.bb) %*% t(w)), 
                   silent = T)
        
        # if there is error in solve(e.bb), skip this loop
        # sometimes exp() function makes the matrix unsolvable
        if ("try-error" %in% class(I.i)) {
          next
        }
        
        var <- (mean(I.i^2) - (n.r / n * (mu.1l - mu.0l) + n.l / n * (mu.1r - mu.0r))^2 / (n.l * n.r / n^2)) / n
        
        # If the variance is less than 0, will assume the goodness for the split is 0
        if (var < 0) {
          next
        } else {
          goodness[i] <- (((mu.1l - mu.0l) - (mu.1r - mu.0r))/sqrt(var))^2
        }
        
      }
      
      
    } else { # if (parms$propsc.mod.loc != "split") { else
      
      prop.sc  <- rep(0, n)
      var.term <- rep(0, n) 
      
      for (i in 1:(length(ux) - 1)){
        
        ind.l <- x %in% ux[ord.ux[1:i]]
        ind.r <- x %in% ux[ord.ux[(i+1):length(ux)]]
        
        n.l <- sum(ind.l)
        n.r <- sum(ind.r)
        
        data.node.noy.l <- data.node.noy[ind.l, ]
        data.node.noy.r <- data.node.noy[ind.r, ]
        
        # if all observations are treated or not treated in the son nodes,
        # skip this loop
        if ((length(unique(data.node.noy.l$A)) == 1) | (length(unique(data.node.noy.r$A)) == 1)) {
          next
        }
        
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
        
        tmp.l <- withWarnings(est.prop.sc(df.noy    = tmp.l$df.noy.fullrank,
                                          method    = parms$propsc.mthd,
                                          form.true = tmp.l$propsc.form.true.updated))
        tmp.r <- withWarnings(est.prop.sc(df.noy    = tmp.r$df.noy.fullrank,
                                          method    = parms$propsc.mthd,
                                          form.true = tmp.r$propsc.form.true.updated))
        w.l <- tmp.l$value$w
        w.r <- tmp.r$value$w
        
        # if 1) propensity score model fitting produces warning or 2) only produce thresholded propensity scores
        # use outside model fitting results
        if (!is.null(tmp.l$warnings)) {
          
          cond.warnings <- T
          for (warnings.i in 1:length(tmp.l$warnings)) {
            cond.warnings <- cond.warnings & (tmp.l$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          # when there is only rank-deficient warning, use the local prediction
          if (cond.warnings) {
            prop.sc.l <- tmp.l$value$prop.sc
          } else {
            # since whl.propsc already subsetted and sorted
            prop.sc.l <- whl.propsc[ind.l]
          }
          
        } else if (identical(sort(unique(tmp.l$value$prop.sc)), c(0.1, 0.9))) {
          prop.sc.l <- whl.propsc[ind.l]
        } else {
          prop.sc.l <- tmp.l$value$prop.sc
        }
        
        if (!is.null(tmp.r$warnings)) {
          
          cond.warnings <- T
          for (warnings.i in 1:length(tmp.r$warnings)) {
            cond.warnings <- cond.warnings & (tmp.r$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
          }
          
          # when there is only rank-deficient warning, use the local prediction
          if (cond.warnings) {
            prop.sc.r <- tmp.r$value$prop.sc
          } else {
            # since whl.propsc already subsetted and sorted
            prop.sc.r <- whl.propsc[ind.r]
          }
          
        } else if (identical(sort(unique(tmp.r$value$prop.sc)), c(0.1, 0.9))) {
          prop.sc.r <- whl.propsc[ind.r]
        } else {
          prop.sc.r <- tmp.r$value$prop.sc
        }
        
        # Skip this split if there is only treated/untreated unit in one of the levels
        # if (length(unique(data.node.noy.l$A)) == 1) {
        #   mu.l <- avg.trt.effct
        # } else {
        mu.l <- mean((Y * A)[ind.l] / prop.sc.l) - mean((Y * (1 - A))[ind.l] / (1 - prop.sc.l))
        # }
        
        # if (length(unique(data.node.noy.r$A)) == 1) {
        #   mu.r <- avg.trt.effct
        # } else {
        mu.r <- mean((Y * A)[ind.r] / prop.sc.r) - mean((Y * (1 - A))[ind.r]/ (1 - prop.sc.r)) 
        # }
        
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
        
        term.l <- try(as.numeric((A[ind.l] - prop.sc.l) * h.l %*% solve(e.bb.l) %*% t(w.l)), silent = T)
        term.r <- try(as.numeric((A[ind.r] - prop.sc.r) * h.r %*% solve(e.bb.r) %*% t(w.r)), silent = T)
        
        # if there is error in solve(e.bb.l) or solve(e.bb.r), skip this loop
        # sometimes exp() function makes the matrix unsolvable
        if (("try-error" %in% class(term.l)) | ("try-error" %in% class(term.r))) {
          next
        }
        
        var.term[ind.l] <- term.l
        var.term[ind.r] <- term.r
        
        I.i <- ((ind.l * Y * A) / (n.l / n * prop.sc) - (ind.l * Y * (1 - A)) / (n.l / n * (1 - prop.sc))) -
          ((ind.r * Y * A) / (n.r / n * prop.sc) - (ind.r * Y * (1 - A)) / (n.r / n * (1 - prop.sc))) - 
          (mu.l - mu.r) - var.term
        var <- (mean(I.i^2) - (n.r / n * mu.l + n.l / n * mu.r)^2 / (n.l * n.r / n^2)) / n
        
        # If the variance is less than 0, will assume the goodness for the split is 0
        if (var < 0) {
          next
        } else {
          goodness[i] <- ((mu.l - mu.r)/sqrt(var))^2
        }
        
      }
      
    }
    
    direction <- ux[ord.ux]
    
  }
  
  list(goodness  = goodness,
       direction = direction)
  
}



etemp.ipw <- function(y, wt, parms) {
  
  n <- length(y)
  
  # sub.ind <- match(y, parms$response)
  sub.ind <- y
  sub.x   <- parms$covariates[sub.ind, ]
  A       <- parms$trt[sub.ind]
  Y       <- parms$response[sub.ind]
  
  data.node = data.frame(Y, A, sub.x)
  data.node.noy <- data.node[, -1]
  
  whl.propsc <- parms$whl.propsc[sub.ind]
  
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
  
  # when there are both treated and untreated unit in split, use avg.trt.effct
  if (length(unique(data.node.noy$A)) == 1) {
    return(list(label = parms$avg.trt.effct, deviance = rss))
    
  } else {
    
    if (parms$propsc.mod.loc == "out") {                    # Model outside node or not
      prop.sc <- whl.propsc
      
    } else {
      
      tmp <- gen.fullrank.ipw(df.noy           = data.node.noy, 
                              propsc.form.true = parms$form.true)
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
  }
  
  mu.1 <- mean(Y * A/ prop.sc)
  mu.0 <- mean(Y * (1 - A)/ (1 - prop.sc))
  avg.trt.effct <- mu.1 - mu.0
  
  list(label = avg.trt.effct, deviance = rss)
}
