######################################################################################################################
################################################# CV method 1, IPW ###################################################
######################################################################################################################
EstIpw.CvMethod1 = function(data.used, tree.list, lambda.used, val.sample, type.var = "cont",
                            propsc.mod.loc, propsc.mthd = "GLM", propsc.form.true = NULL){
  
  # type.var does not matter for IPW, the parameter can be omitted, but the previous code needs to be modified
  # min.obs.mod: minimum number of observations for the model to be fitted
  
  # if (is.null(val.w)){                          # True model or not, assign value to w
  #   val.w <- cbind(1, val.sample[, 3:ncol(val.sample)])
  # } 
  val.sample.noy <- val.sample[, !colnames(val.sample) %in% c("Y")]
  whl.propsc <- gen.fullrank.ipw(df.noy           = val.sample.noy, 
                                 propsc.form.true = propsc.form.true)
  whl.propsc <- withWarnings(est.prop.sc(df.noy    = whl.propsc$df.noy.fullrank,
                                         method    = propsc.mthd,
                                         form.true = whl.propsc$propsc.form.true.updated))
  # avg.trt.effct <- mean((val.sample$Y * val.sample$A / whl.propsc$prop.sc)) -
  #   mean((val.sample$Y * (1 - val.sample$A)/ (1 - whl.propsc$prop.sc)))
  
  # For in node and in split, whole dataset
  val.w <- whl.propsc$value$w
  
  if (propsc.mod.loc == "out"){
    
    val.sample$prop.sc <- whl.propsc$value$prop.sc
    # val.w <- whl.propsc$w
    
  }
  
  # Storage space for the complexity value associated with each candidate tree
  complex.val = rep(NA, length(tree.list))
  
  # Looping through the candidate trees
  for (m in 1:length(complex.val)){
    # for (m in c(1)){
    # Finding the tree being evaluated
    # print(m)
    tree.used = tree.list[[m]] 
    
    # If only root node there is no internal node
    if(nrow(tree.used$frame) == 1){
      goodness.test = 0
      numb.int = 0
      
    } else {            # If at least one split
      
      is.leaf <- (tree.used$frame$var == "<leaf>")
      goodness.test = 0
      # Calculate the goodness of the tree using the test
      # Finding the test data falling in each terminal node
      numb.int = sum(!is.leaf)
      # Finding all kids on terminal node 
      #kids.i <- order(as.numeric(rownames(tree.used$frame)))
      
      # Calculating split complexity for each internal node of branch
      right.ind  <- 1
      last.right <- NULL
      w.last.right <- NULL
      
      for (h in 1:dim(tree.used$frame)[1]){
        # for (h in 1:385){
        # h <- h + 1
        
        # Finding observations in validation sample falling in that node
        if (tree.used$frame$var[h] == "<leaf>") {
          next
        } else if (tree.used$frame$n[h] == dim(data.used)[1]){
          val.sample.used <- val.sample
          val.w.used      <- val.w
        } else{
          if (!is.null(last.left)){
            val.sample.used <- last.left[[1]]
            val.w.used      <- w.last.left[[1]]
          } else {
            val.sample.used <- as.data.frame(last.right[[right.ind - 1]])
            val.w.used      <- as.data.frame(w.last.right[[right.ind - 1]])
            if (length(last.right) > 2){
              last.right   <- last.right[1:(right.ind - 2)]
              w.last.right <- w.last.right[1:(right.ind - 2)]
            } else if (length(last.right) == 2){
              last.right   <- list(last.right[[1]])
              w.last.right <- list(w.last.right[[1]])
            }else {
              last.right <- NULL
              w.last.right <- NULL
            }
            right.ind <- right.ind - 1
          }
        }
        
        # Predict propensity score inside node (need to be before identifying obs into right/left node since need propsc later)
        if (propsc.mod.loc == "node"){ 
          
          if (length(unique(val.sample.used$A)) > 1) { 
            val.sample.used.noy     <- val.sample.used[, !(colnames(val.sample.used) %in% c("Y", "prop.sc"))]
            node.propsc.fullrank <- gen.fullrank.ipw(df.noy           = val.sample.used.noy, 
                                                     propsc.form.true = propsc.form.true)
            node.propsc <- withWarnings(est.prop.sc(df.noy    = node.propsc.fullrank$df.noy.fullrank,
                                                    method    = propsc.mthd,
                                                    form.true = node.propsc.fullrank$propsc.form.true.updated))
            # val.w.used will be the same from both following conditions
            val.w.used              <- node.propsc$value$w
            
            # if 1) propensity score model fitting produces warning or 2) only produce thresholded propensity scores
            # use outside model fitting results
            if (!is.null(node.propsc$warnings)) {
              
              cond.warnings <- T
              for (warnings.i in 1:length(node.propsc$warnings)) {
                cond.warnings <- cond.warnings & (node.propsc$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
              }
              # when there is only rank-deficient warning, use the local prediction
              if (cond.warnings) {
                val.sample.used$prop.sc <- node.propsc$value$prop.sc
              } else {
                val.sample.used$prop.sc <- whl.propsc$value$prop.sc[rownames(val.sample) %in% rownames(val.sample.used)]
              }
              
            } else if (identical(sort(unique(node.propsc$value$prop.sc)), c(0.1, 0.9))) {
              val.sample.used$prop.sc <- whl.propsc$value$prop.sc[rownames(val.sample) %in% rownames(val.sample.used)]
            } else {
              val.sample.used$prop.sc <- node.propsc$value$prop.sc
              # val.w.used              <- node.propsc$value$w
            }
          }
          #   } # end if (!propsc.mod.insplt){
        } # end if (propsc.mod.loc == "node"){ 
        
        row.ind <- sum((tree.used$frame$var[1:h]) != "<leaf>")
        split.used <- tree.used$splits[row.ind, 4]
        var.used <- tree.used$frame$var[h]
        col.ind <- which(colnames(val.sample.used) == var.used)
        
        # Calculate goodness corresponding to split
        # Finding observations falling in right and left node
        # Need to figure out observations going to the left/right node
        # if ((split.used %% 1) == 0){   
        if (class(data.used[, as.character(var.used)]) == "factor"){   
          # Categorical covariate split
          lvls <- levels(val.sample[, col.ind])
          val.sample.left  <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used, ] == 1], ]
          val.sample.right <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          val.w.used.left  <- val.w.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used, ] == 1], ]
          val.w.used.right <- val.w.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          
        } else{
          # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[row.ind, 2] > 0) {
            val.sample.left  <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
            val.w.used.left  <- val.w.used[val.sample.used[,  col.ind] >= split.used, ]
            val.w.used.right <- val.w.used[val.sample.used[,  col.ind] < split.used, ]
            
          } else {
            val.sample.left <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
            val.w.used.left  <- val.w.used[val.sample.used[,  col.ind] < split.used, ]
            val.w.used.right <- val.w.used[val.sample.used[,  col.ind] >= split.used, ]
          }
          
        }
        
        if (tree.used$frame$var[h+1] != "<leaf>"){
          last.left   <- list(val.sample.left)
          w.last.left <- list(val.w.used.left)
        } else{
          last.left   <- NULL
          w.last.left <- NULL
        }
        
        which.right <- as.numeric(rownames(tree.used$frame)[h+1]) + 1
        if (tree.used$frame$var[as.numeric(rownames(tree.used$frame)) == which.right] != "<leaf>"){
          last.right[[right.ind]]   <- val.sample.right
          w.last.right[[right.ind]] <- val.w.used.right
          right.ind <- right.ind + 1
        } 
        
        # Skip if there is only treated/untreated unit in the node when the model is fitted in node
        # Since if replace with average treatment effect, goodness.test will be 0 anyway
        # if (!propsc.mod.out){
        #   if (!propsc.mod.insplt){
        #     if (length(unique(val.sample.used$A)) <= 1){
        #       next
        #     }
        #   } # end if (!propsc.mod.insplt){
        # } # end if (!propsc.mod.out){
        
        # Skip if the number of observations in any split is smaller than the dimension of covariates to avoid variance computation error
        # if (!is.null(val.w.used)){
        #   if ((nrow(val.sample.left) < ncol(val.w.used)) | (nrow(val.sample.right) < ncol(val.w.used))){
        #     next             
        #   }
        # }
        
        # Skip if there is no observations in the left/right subgroup
        # if (min(nrow(val.sample.left), nrow(val.sample.right)) == 0) {
        #   next
        # }
        
        # Skip if there is no treated/untreated in the left/right subgroup
        if (min(sum(val.sample.left$A == 1), sum(val.sample.left$A == 0), sum(val.sample.right$A == 1), sum(val.sample.right$A == 0)) == 0){
          next
        }
        
        if (propsc.mod.loc != "split") {
          
          mu.1l <- mean(val.sample.left$Y * val.sample.left$A / val.sample.left$prop.sc)
          mu.0l <- mean(val.sample.left$Y * (1 - val.sample.left$A)/ (1 - val.sample.left$prop.sc))
          mu.1r <- mean(val.sample.right$Y * val.sample.right$A / val.sample.right$prop.sc)
          mu.0r <- mean(val.sample.right$Y * (1 - val.sample.right$A)/ (1 - val.sample.right$prop.sc))
          
          mu.l <- mu.1l - mu.0l
          mu.r <- mu.1r - mu.0r
          
          # When the propensity score is fitted outside and there are categorical covariates
          # need to take care of the singular inverse
          if (propsc.mod.loc == "out") {
            val.sample.used.noy     <- val.sample.used[, !(colnames(val.sample.used) %in% c("Y", "prop.sc"))]
            node.propsc.fullrank <- gen.fullrank.ipw(df.noy           = val.sample.used.noy, 
                                                     propsc.form.true = propsc.form.true)
            node.propsc <- withWarnings(est.prop.sc(df.noy    = node.propsc.fullrank$df.noy.fullrank,
                                                    method    = propsc.mthd,
                                                    form.true = node.propsc.fullrank$propsc.form.true.updated))
            # val.w.used will be the same from both following conditions
            val.w.used       <- node.propsc$value$w
            val.w.used.left  <- val.w.used[rownames(val.w.used) %in% rownames(val.sample.left), ]
            val.w.used.right <- val.w.used[rownames(val.w.used) %in% rownames(val.sample.right), ]
          }
          
          h.l <- apply((val.sample.left$A * val.sample.left$Y * (1 - val.sample.left$prop.sc) / val.sample.left$prop.sc + 
                          (1 - val.sample.left$A) * val.sample.left$Y * val.sample.left$prop.sc / (1 - val.sample.left$prop.sc)) * 
                         val.w.used.left,
                       2,
                       mean)
          h.r <- apply((val.sample.right$A * val.sample.right$Y * (1 - val.sample.right$prop.sc) / val.sample.right$prop.sc + 
                          (1 - val.sample.right$A) * val.sample.right$Y * val.sample.right$prop.sc / (1 - val.sample.right$prop.sc)) * 
                         val.w.used.right,
                       2,
                       mean)
          
          n.l    <- nrow(val.sample.left)
          n.r    <- nrow(val.sample.right)
          n.used <- nrow(val.sample.used)
          
          e.bb <- as.matrix(t(val.sample.used$prop.sc * val.w.used)) %*% 
            as.matrix((1 - val.sample.used$prop.sc) * val.w.used) / n.used
          
          I.i <- try((((rownames(val.sample.used) %in% rownames(val.sample.left)) * val.sample.used$Y * val.sample.used$A) / 
                        (n.l / n.used * val.sample.used$prop.sc) - 
                        ((rownames(val.sample.used) %in% rownames(val.sample.left)) * val.sample.used$Y * (1 - val.sample.used$A)) / 
                        (n.l / n.used * (1 - val.sample.used$prop.sc))) -
                       (((rownames(val.sample.used) %in% rownames(val.sample.right)) * val.sample.used$Y * val.sample.used$A) / 
                          (n.r / n.used * val.sample.used$prop.sc) - 
                          ((rownames(val.sample.used) %in% rownames(val.sample.right)) * val.sample.used$Y * (1 - val.sample.used$A)) / 
                          (n.r / n.used * (1 - val.sample.used$prop.sc))) - 
                       (mu.l - mu.r) - 
                       as.numeric((val.sample.used$A - val.sample.used$prop.sc) * (h.l - h.r) %*% solve(e.bb) %*% t(val.w.used)), 
                     silent = T)
          
          # When the number of observations in val.w.used very small, cause problem when solve(e.bb)
          # will skip this loop, so the addition to goodness.test will be 0
          if ("try-error" %in% class(I.i)) {
            next
          }
          
          # var <- (mean(I.i^2) - 
          #           (n.r / n.used * mu.l + n.l / n.used * mu.r)^2 / 
          #           (n.l * n.r / n.used^2)) / n.used
          
        } else { # if (propsc.mod.loc != "split") { else
          
          # if ((nrow(val.sample.left) < min.obs.mod) & (nrow(val.sample.right) < min.obs.mod)) {
          #   next                         # break if the number of observations in the son nodes is fewer than min.obs.mod
          # } else if ((nrow(val.sample.left) == 1) | (nrow(val.sample.right) == 1)) {
          #   next
          # }
          
          val.sample.noy.l  <- val.sample.left[, !(colnames(val.sample.left) %in% c("Y", "prop.sc"))]
          tmp.l.fullrank <- gen.fullrank.ipw(df.noy           = val.sample.noy.l, 
                                             propsc.form.true = propsc.form.true)
          tmp.l <- withWarnings(est.prop.sc(df.noy    = tmp.l.fullrank$df.noy.fullrank,
                                            method    = propsc.mthd,
                                            form.true = tmp.l.fullrank$propsc.form.true.updated))
          val.w.used.left  <- tmp.l$value$w
          
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
              prop.sc.l <- whl.propsc$value$prop.sc[rownames(val.sample) %in% rownames(val.sample.left)]
            }
            
          } else if (identical(sort(unique(tmp.l$value$prop.sc)), c(0.1, 0.9))) {
            prop.sc.l <- whl.propsc$value$prop.sc[rownames(val.sample) %in% rownames(val.sample.left)]
          } else {
            prop.sc.l <- tmp.l$value$prop.sc
          }
          
          val.sample.noy.r <- val.sample.right[, !(colnames(val.sample.right) %in% c("Y", "prop.sc"))]
          tmp.r.fullrank <- gen.fullrank.ipw(df.noy           = val.sample.noy.r, 
                                             propsc.form.true = propsc.form.true)
          tmp.r <- withWarnings(est.prop.sc(df.noy    = tmp.r.fullrank$df.noy.fullrank,
                                            method    = propsc.mthd,
                                            form.true = tmp.r.fullrank$propsc.form.true.updated))
          val.w.used.right <- tmp.r$value$w
          
          # if 1) propensity score model fitting produces warning or 2) only produce thresholded propensity scores
          # use outside model fitting results
          if (!is.null(tmp.r$warnings)) {
        
            cond.warnings <- T
            for (warnings.i in 1:length(tmp.r$warnings)) {
              cond.warnings <- cond.warnings & (tmp.r$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
            }
            # when there is only rank-deficient warning, use the local prediction
            if (cond.warnings) {
              prop.sc.r <- tmp.r$value$prop.sc
            } else {
              prop.sc.r <- whl.propsc$value$prop.sc[rownames(val.sample) %in% rownames(val.sample.right)]
            }
            
          } else if (identical(sort(unique(tmp.r$value$prop.sc)), c(0.1, 0.9))) {
            prop.sc.r <- whl.propsc$value$prop.sc[rownames(val.sample) %in% rownames(val.sample.right)]
          } else {
            prop.sc.r <- tmp.r$value$prop.sc
          }
          
          # Skip if there are fewer observations than the number of columns in w
          # if ((nrow(val.sample.left) < ncol(val.w.used.left)) | (nrow(val.sample.right) < ncol(val.w.used.right))){
          #   next             
          # }
          
          # map prop.sc in the left and right child nodes back to the parent node
          val.sample.used$prop.sc[rownames(val.sample.used) %in% rownames(val.sample.left)] <- prop.sc.l
          val.sample.used$prop.sc[rownames(val.sample.used) %in% rownames(val.sample.right)] <- prop.sc.r
          
          prop.sc <- val.sample.used$prop.sc
          # remove the propensity score column for later if conditions
          val.sample.used <- val.sample.used[, colnames(val.sample.used) != "prop.sc"]
          
          # Numerator in the splitting statistics
          # replace the numerator with average treatment effect from the whole dataset 
          # when there is only treated/untreated unit in the node
          # if (length(unique(val.sample.left$A)) == 1) {
          #   mu.l <- avg.trt.effct
          # } else {
          mu.l <- mean(val.sample.left$Y * val.sample.left$A / prop.sc.l) - 
            mean(val.sample.left$Y * (1 - val.sample.left$A)/ (1 - prop.sc.l))
          # }
          
          # 1st condition is skipped earlier
          # if (length(unique(val.sample.right$A)) == 1) {
          #   mu.r <- avg.trt.effct
          # } else {
          mu.r <- mean(val.sample.right$Y * val.sample.right$A / prop.sc.r) -
            mean(val.sample.right$Y * (1 - val.sample.right$A)/ (1 - prop.sc.r))
          # }
          
          h.l <- apply((val.sample.left$A * val.sample.left$Y * (1 - prop.sc.l) / prop.sc.l + 
                          (1 - val.sample.left$A) * val.sample.left$Y * prop.sc.l / (1 - prop.sc.l)) * 
                         val.w.used.left,
                       2,
                       mean)
          h.r <- apply((val.sample.right$A * val.sample.right$Y * (1 - prop.sc.r) / prop.sc.r + 
                          (1 - val.sample.right$A) * val.sample.right$Y * prop.sc.r / (1 - prop.sc.r)) * 
                         val.w.used.right,
                       2,
                       mean)
          
          n.l    <- nrow(val.sample.left)
          n.r    <- nrow(val.sample.right)
          n.used <- nrow(val.sample.used)
          
          e.bb.l <- as.matrix(t(prop.sc.l * val.w.used.left)) %*%
            as.matrix((1 - prop.sc.l) * val.w.used.left) / n.used
          e.bb.r <- as.matrix(t(prop.sc.r * val.w.used.right)) %*%
            as.matrix((1 - prop.sc.r) * val.w.used.right) / n.used
          
          term.l <- try(as.numeric((val.sample.left$A - prop.sc.l) * h.l %*% solve(e.bb.l) %*% t(val.w.used.left)), silent = T)
          term.r <- try(as.numeric((val.sample.right$A - prop.sc.r) * h.r %*% solve(e.bb.r) %*% t(val.w.used.right)), silent = T)
          
          if (("try-error" %in% class(term.l)) | ("try-error" %in% class(term.r))) {
            next
          }
          
          val.sample.used$var.term[rownames(val.sample.used) %in% rownames(val.sample.left)]  <- term.l
          val.sample.used$var.term[rownames(val.sample.used) %in% rownames(val.sample.right)] <- term.r
          
          var.term <- val.sample.used$var.term
          val.sample.used <- val.sample.used[, colnames(val.sample.used) != "var.term"]
          
          I.i <- (((rownames(val.sample.used) %in% rownames(val.sample.left)) * val.sample.used$Y * val.sample.used$A) / 
                    (n.l / n.used * prop.sc) - 
                    ((rownames(val.sample.used) %in% rownames(val.sample.left)) * val.sample.used$Y * (1 - val.sample.used$A)) / 
                    (n.l / n.used * (1 - prop.sc))) -
            (((rownames(val.sample.used) %in% rownames(val.sample.right)) * val.sample.used$Y * val.sample.used$A) / 
               (n.r / n.used * prop.sc) - 
               ((rownames(val.sample.used) %in% rownames(val.sample.right)) * val.sample.used$Y * (1 - val.sample.used$A)) / 
               (n.r / n.used * (1 - prop.sc))) - 
            (mu.l - mu.r) - var.term
          
          # var <- (mean(I.i^2) - 
          #           (n.r / n.used * mu.l + n.l / n.used * mu.r)^2 / 
          #           (n.l * n.r / n.used^2)) / n.used
        }
        
        var <- (mean(I.i^2) - 
                  (n.r / n.used * mu.l + n.l / n.used * mu.r)^2 / 
                  (n.l * n.r / n.used^2)) / n.used
        
        # mu.l - mu.r
        # ((mu.l - mu.r) / sqrt(var))^2
        # print(var.l + var.r)
        if (var < 0){
          next
        } else{
          goodness.test <- goodness.test + ((mu.l - mu.r) / sqrt(var))^2
        }
        
        # print(var)
      } # End h
    } # End if loop
    # Calculating complexity value
    complex.val[m] = goodness.test - lambda.used * numb.int
  } # End m loop
  
  # Averaging over cross validation sets
  tree.final = tree.list[[which.max(complex.val)]]
  
  return(list(tree.final, complex.val))
}


######################################################################################################################
############################################ CV method 1, G estimation ###############################################
######################################################################################################################
EstG.CvMethod1 = function(data.used, tree.list, lambda.used, val.sample, type.var = "cont",
                          adj.mod.loc, adj.mthd = "GLM", adj.form.true = NULL) {
  
  # min.obs.mod: minimum number of observations for the model to be fitted
  
  # if (is.null(val.w)){                          # True model or not, assign value to w
  #   val.w <- cbind(1, val.sample$A, val.sample[, 3:ncol(val.sample)], val.sample$A * val.sample[, 3:ncol(val.sample)])
  #   colnames(val.w) <- c("(Intercept)", "A", colnames(val.sample)[3:ncol(val.sample)], paste("A:", colnames(val.sample)[3:ncol(val.sample)], sep = ""))
  # }
  
  whl.g <- gen.fullrank.g(df            = val.sample,
                          adj.form.true = adj.form.true)
  whl.g <- withWarnings(est.cond.eff(df        = whl.g$df.fullrank,
                                     method    = adj.mthd,
                                     form.true = whl.g$adj.form.true.updated,
                                     type.var  = type.var))
  # For in node and in split, whole dataset
  val.w <- NULL
  
  if (adj.mod.loc == "out"){
    
    val.sample$est.cond.eff.0 <- whl.g$value$pred.A.0
    val.sample$est.cond.eff.1 <- whl.g$value$pred.A.1
    var.rb                    <- whl.g$value$var.rb
    val.w                     <- whl.g$value$w
    # if (ncol(var.rb) < ncol(val.w)){
    #   val.w <- val.w %>%
    #     select(colnames(var.rb))
    # }
  }
  
  # Storage space for the complexity value associated with each candidate tree
  complex.val = rep(NA, length(tree.list))
  
  # Looping through the candidate trees
  for (m in 1:length(complex.val)){
    # for (m in c(1)){
    # Finding the tree being evaluated
    # print(m)
    tree.used = tree.list[[m]] 
    
    # If only root node there is no internal node
    if(nrow(tree.used$frame) == 1){
      goodness.test = 0
      numb.int = 0
      
    } else {            # If at least one split
      
      is.leaf <- (tree.used$frame$var == "<leaf>")
      goodness.test = 0
      # Calculate the goodness of the tree using the test
      # Finding the test data falling in each terminal node
      numb.int = sum(!is.leaf)
      # Finding all kids on terminal node 
      #kids.i <- order(as.numeric(rownames(tree.used$frame)))
      
      # Calculating split complexity for each internal node of branch
      right.ind  <- 1
      last.right <- NULL
      w.last.right <- NULL
      
      for (h in 1:dim(tree.used$frame)[1]){
        # for (h in 3:129){
        # h <- h + 1
        
        # Finding observations in validation sample falling in that node
        if (tree.used$frame$var[h] == "<leaf>") {
          next
        } else if (tree.used$frame$n[h] == dim(data.used)[1]){
          val.sample.used <- val.sample
          val.w.used      <- val.w
        } else{
          
          if (!is.null(last.left)){
            val.sample.used <- last.left[[1]]
            val.w.used      <- w.last.left[[1]]
          } else {
            val.sample.used <- as.data.frame(last.right[[right.ind - 1]])
            val.w.used      <- as.data.frame(w.last.right[[right.ind - 1]])
            if (length(last.right) > 2){
              last.right   <- last.right[1:(right.ind - 2)]
              w.last.right <- w.last.right[1:(right.ind - 2)]
            } else if (length(last.right) == 2){
              last.right   <- list(last.right[[1]])
              w.last.right <- list(w.last.right[[1]])
            }else {
              last.right <- NULL
              w.last.right <- NULL
            }
            right.ind <- right.ind - 1
          }
          
        }
        
        # Adjustment model inside node (need to be before identifying obs into right/left node since need propsc later)
        if (adj.mod.loc == "node") { 
          # if (!adj.mod.insplt) {
          if (length(unique(val.sample.used$A)) > 1) { 
            
            val.sample.used.nocondeff <- val.sample.used[, !(colnames(val.sample.used) %in% c("est.cond.eff.0", "est.cond.eff.1"))]
            tmp <- gen.fullrank.g(df            = val.sample.used.nocondeff,
                                  adj.form.true = adj.form.true)
            
            tmp <- withWarnings(est.cond.eff(df        = tmp$df.fullrank,
                                             method    = adj.mthd,
                                             form.true = tmp$adj.form.true.updated,
                                             type.var  = type.var))
            var.rb     <- tmp$value$var.rb
            val.w.used <- tmp$value$w
            
            # when the outcome model fitting produces warning
            # use outside model fitting results for conditional means
            if (!is.null(tmp$warnings)) {
              
              # when there is only rank-deficient warning, use the local prediction
              cond.warnings <- T
              for (warnings.i in 1:length(tmp$warnings)) {
                cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
              }
              
              if (cond.warnings) {
                val.sample.used$est.cond.eff.0 <- tmp$value$pred.A.0
                val.sample.used$est.cond.eff.1 <- tmp$value$pred.A.1
              } else {
                val.sample.used$est.cond.eff.0 <- whl.g$value$pred.A.0[rownames(val.sample) %in% rownames(val.sample.used)]
                val.sample.used$est.cond.eff.1 <- whl.g$value$pred.A.1[rownames(val.sample) %in% rownames(val.sample.used)]
              }
              
            } else {
              val.sample.used$est.cond.eff.0 <- tmp$value$pred.A.0
              val.sample.used$est.cond.eff.1 <- tmp$value$pred.A.1
            }
            
          } 
          # } # end if (!propsc.mod.insplt){
        } # end if (adj.mod.loc == "node") { 
        
        row.ind <- sum((tree.used$frame$var[1:h]) != "<leaf>")
        split.used <- tree.used$splits[row.ind, 4]
        var.used <- tree.used$frame$var[h]
        col.ind <- which(colnames(val.sample.used) == var.used)
        
        # Calculate goodness corresponding to split
        # Finding observations falling in right and left node
        # Need to figure out observations going to the left/right node
        # if ((split.used %% 1) == 0){   
        if (class(data.used[, as.character(var.used)]) == "factor"){   
          
          # Categorical covariate split
          lvls <- levels(val.sample[, col.ind])
          val.sample.left  <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used, ] == 1], ]
          val.sample.right <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          val.w.used.left  <- val.w.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used, ] == 1], ]
          val.w.used.right <- val.w.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          
        } else{
          # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[row.ind, 2] > 0) {
            val.sample.left  <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
            val.w.used.left  <- val.w.used[val.sample.used[,  col.ind] >= split.used, ]
            val.w.used.right <- val.w.used[val.sample.used[,  col.ind] < split.used, ]
            
          } else {
            val.sample.left <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
            
            # # Need to take care for w since it is not a dataframe if there is only one row in the node in the validation sample
            # if (class(val.w.used) == "numeric") {
            #   val.w.used.left  <- data.frame(t(val.w.used))[val.sample.used[,  col.ind] < split.used, ]
            #   val.w.used.right <- data.frame(t(val.w.used))[val.sample.used[,  col.ind] >= split.used, ]
            # } else {
            val.w.used.left  <- val.w.used[val.sample.used[,  col.ind] < split.used, ]
            val.w.used.right <- val.w.used[val.sample.used[,  col.ind] >= split.used, ]
            # }
            
          }
          
        }
        
        if (tree.used$frame$var[h+1] != "<leaf>"){
          last.left   <- list(val.sample.left)
          w.last.left <- list(val.w.used.left)
        } else{
          last.left   <- NULL
          w.last.left <- NULL
        }
        
        which.right <- as.numeric(rownames(tree.used$frame)[h+1]) + 1
        if (tree.used$frame$var[as.numeric(rownames(tree.used$frame)) == which.right] != "<leaf>"){
          last.right[[right.ind]]   <- val.sample.right
          w.last.right[[right.ind]] <- val.w.used.right
          right.ind <- right.ind + 1
        } 
        
        # Skip if there is only treated/untreated unit in the node when the model is fitted in node
        # the goodness will be 0 anyway
        # if (!adj.mod.out){
        #   if (!adj.mod.insplt) {
        #     if (length(unique(val.sample.used$A)) <= 1){
        #       next
        #     }
        #   } 
        # else { # if (!adj.mod.insplt){ this else is for models fitted in split
        #   if ((nrow(val.sample.left) < min.obs.mod) & (nrow(val.sample.right) < min.obs.mod)) {
        #     next                         # break if the number of observations in both son nodes is fewer than min.obs.mod
        #   } else if ((nrow(val.sample.left) == 1) | (nrow(val.sample.right) == 1)) {
        #     next
        #   }
        # } # end if (!adj.mod.insplt){
        # } # end if (!adj.mod.out){
        
        # Skip if the number of observations in any split is smaller than the dimension of covariates to avoid variance computation error
        # if (!is.null(val.w.used)) {
        #   if ((nrow(val.sample.left) < ncol(val.w.used)) | (nrow(val.sample.right) < ncol(val.w.used))){
        #     next             
        #   }
        # }
        
        # if (min(nrow(val.sample.left), nrow(val.sample.right)) == 0) {
        #   next
        # }
        
        # Skip if there is no treated/untreated in the left/right subgroup
        if (min(sum(val.sample.left$A == 1), sum(val.sample.left$A == 0), sum(val.sample.right$A == 1), sum(val.sample.right$A == 0)) == 0){
          next
        }
        
        if (adj.mod.loc != "split") {
          
          mu.1l <- mean(val.sample.left$est.cond.eff.1)
          mu.0l <- mean(val.sample.left$est.cond.eff.0)
          mu.1r <- mean(val.sample.right$est.cond.eff.1)
          mu.0r <- mean(val.sample.right$est.cond.eff.0)
          
          # if (ncol(var.rb) < ncol(val.w.used.left)){
          #   val.w.used.left <- val.w.used.left %>%
          #     select(colnames(var.rb))
          # }
          # 
          # if (ncol(var.rb) < ncol(val.w.used.right)){
          #   val.w.used.right <- val.w.used.right %>%
          #     select(colnames(var.rb))
          # }
          
          val.w.1l <- val.w.used.left %>% mutate(A = 1)
          val.w.0l <- val.w.used.left %>% mutate(A = 0)
          val.w.1r <- val.w.used.right %>% mutate(A = 1)
          val.w.0r <- val.w.used.right %>% mutate(A = 0)
          
          if (type.var == "cont") {
            
            g.b.l.1 <- apply(val.w.1l, 2, mean)
            g.b.r.1 <- apply(val.w.1r, 2, mean)
            g.b.l.0 <- apply(val.w.0l, 2, mean)
            g.b.r.0 <- apply(val.w.0r, 2, mean)
            
          } else if (type.var == "bin") {
            
            g.b.l.1 <- apply(val.w.1l * val.sample.left$est.cond.eff.1 * (1 - val.sample.left$est.cond.eff.1), 2, mean)
            g.b.r.1 <- apply(val.w.1r * val.sample.right$est.cond.eff.1 * (1 - val.sample.right$est.cond.eff.1), 2, mean)
            g.b.l.0 <- apply(val.w.0l * val.sample.left$est.cond.eff.0 * (1 - val.sample.left$est.cond.eff.0), 2, mean)
            g.b.r.0 <- apply(val.w.0r * val.sample.right$est.cond.eff.0 * (1 - val.sample.right$est.cond.eff.0), 2, mean)
            
          }
          
          left.p.1 <- as.matrix(data.frame(val.sample.left$est.cond.eff.1))
          left.p.0 <- as.matrix(data.frame(val.sample.left$est.cond.eff.0))
          right.p.1 <- as.matrix(data.frame(val.sample.right$est.cond.eff.1))
          right.p.0 <- as.matrix(data.frame(val.sample.right$est.cond.eff.0))
          
          n.l <- nrow(val.sample.left)
          n.r <- nrow(val.sample.right)
          
          # Calculate variance estimators
          var.1l <- g.b.l.1 %*% var.rb %*% g.b.l.1 + 1/n.l^2 * sum( (left.p.1 - mu.1l)^2 )
          var.1r <- g.b.r.1 %*% var.rb %*% g.b.r.1 + 1/n.r^2 * sum( (right.p.1 - mu.1r)^2 )
          var.0l <- g.b.l.0 %*% var.rb %*% g.b.l.0 + 1/n.l^2 * sum( (left.p.0 - mu.0l)^2 )
          var.0r <- g.b.r.0 %*% var.rb %*% g.b.r.0 + 1/n.r^2 * sum( (right.p.0 - mu.0r)^2 )
          
          # only in node model fitting
          if (adj.mod.loc == "node") {
            rm(var.rb)
          }
          
        } else { #  if (adj.mod.loc != "split") { when adj.mod.insplt = T
          
          # length(unique(val.sample.left$A)) == 1 or 0 is skipped earlier
          # skip if there is only treated/untreated unit in the child node
          # skip since their w will cause problem if only treated/untreated unit
          # if ((length(unique(val.sample.left$A)) == 1) | (length(unique(val.sample.right$A)) == 1)) {
          #   next
          # }
          
          val.sample.left.nocondeff  <- val.sample.left[, !(colnames(val.sample.left) %in% c("est.cond.eff.0", "est.cond.eff.1"))]
          tmp.l <- gen.fullrank.g(df            = val.sample.left.nocondeff,
                                  adj.form.true = adj.form.true)
          tmp.l <- withWarnings(est.cond.eff(df        = tmp.l$df.fullrank,
                                             method    = adj.mthd,
                                             form.true = tmp.l$adj.form.true.updated,
                                             type.var  = type.var))
          var.rb.l        <- tmp.l$value$var.rb
          val.w.used.left <- tmp.l$value$w
          
          # when the outcome model fitting produces warning
          # use outside model fitting results for conditional means
          if (!is.null(tmp.l$warnings)) {
            
            # when there is only rank-deficient warning, use the local prediction
            cond.warnings <- T
            for (warnings.i in 1:length(tmp.l$warnings)) {
              cond.warnings <- cond.warnings & (tmp.l$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
            }
            
            if (cond.warnings) {
              val.sample.left$est.cond.eff.0 <- tmp.l$value$pred.A.0
              val.sample.left$est.cond.eff.1 <- tmp.l$value$pred.A.1
            } else {
              val.sample.left$est.cond.eff.0 <- whl.g$value$pred.A.0[rownames(val.sample) %in% rownames(val.sample.left)]
              val.sample.left$est.cond.eff.1 <- whl.g$value$pred.A.1[rownames(val.sample) %in% rownames(val.sample.left)]
            }
            
          } else {
            val.sample.left$est.cond.eff.0 <- tmp.l$value$pred.A.0
            val.sample.left$est.cond.eff.1 <- tmp.l$value$pred.A.1
          }
          
          
          # For some reason, the var.rb cannot be calculated when the number of observations is small in the child nodes
          # Therefore, skip
          if (all(is.nan(var.rb.l))) {
            next
          }
          
          val.sample.right.nocondeff <- val.sample.right[, !(colnames(val.sample.right) %in% c("est.cond.eff.0", "est.cond.eff.1"))]
          tmp.r <- gen.fullrank.g(df            = val.sample.right.nocondeff,
                                  adj.form.true = adj.form.true)
          tmp.r <- withWarnings(est.cond.eff(df        = tmp.r$df.fullrank,
                                             method    = adj.mthd,
                                             form.true = tmp.r$adj.form.true.updated,
                                             type.var  = type.var))   
          var.rb.r         <- tmp.r$value$var.rb
          val.w.used.right <- tmp.r$value$w
          
          # when the outcome model fitting produces warning
          # use outside model fitting results for conditional means
          if (!is.null(tmp.r$warnings)) {
            
            # when there is only rank-deficient warning, use the local prediction
            cond.warnings <- T
            for (warnings.i in 1:length(tmp.r$warnings)) {
              cond.warnings <- cond.warnings & (tmp.r$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
            }
            
            if (cond.warnings) {
              val.sample.right$est.cond.eff.0 <- tmp.r$value$pred.A.0
              val.sample.right$est.cond.eff.1 <- tmp.r$value$pred.A.1        
            } else {
              val.sample.right$est.cond.eff.0 <- whl.g$value$pred.A.0[rownames(val.sample) %in% rownames(val.sample.right)]
              val.sample.right$est.cond.eff.1 <- whl.g$value$pred.A.1[rownames(val.sample) %in% rownames(val.sample.right)]
            }
            
          } else {
            val.sample.right$est.cond.eff.0 <- tmp.r$value$pred.A.0
            val.sample.right$est.cond.eff.1 <- tmp.r$value$pred.A.1
          }
          
          # For some reason, the var.rb cannot be calculated when the number of observations is small in the child nodes
          # Therefore, skip
          if (all(is.nan(var.rb.r))) {
            next
          }
          
          mu.1l <- mean(val.sample.left$est.cond.eff.1)
          mu.0l <- mean(val.sample.left$est.cond.eff.0)
          mu.1r <- mean(val.sample.right$est.cond.eff.1)
          mu.0r <- mean(val.sample.right$est.cond.eff.0)
          
          # if (ncol(var.rb.l) < ncol(val.w.used.left)) {
          #   val.w.used.left <- val.w.used.left %>%
          #     select(colnames(var.rb.l))
          # }
          # 
          # if (ncol(var.rb.r) < ncol(val.w.used.right)) {
          #   val.w.used.right <- val.w.used.right %>%
          #     select(colnames(var.rb.r))
          # }
          
          val.w.1l <- val.w.used.left %>% mutate(A = 1)
          val.w.0l <- val.w.used.left %>% mutate(A = 0)
          val.w.1r <- val.w.used.right %>% mutate(A = 1)
          val.w.0r <- val.w.used.right %>% mutate(A = 0)
          
          if (type.var == "cont") {
            
            g.b.l.1 <- apply(val.w.1l, 2, mean)
            g.b.r.1 <- apply(val.w.1r, 2, mean)
            g.b.l.0 <- apply(val.w.0l, 2, mean)
            g.b.r.0 <- apply(val.w.0r, 2, mean)
            
          } else if (type.var == "bin") {
            
            g.b.l.1 <- apply(val.w.1l * val.sample.left$est.cond.eff.1 * (1 - val.sample.left$est.cond.eff.1), 2, mean)
            g.b.r.1 <- apply(val.w.1r * val.sample.right$est.cond.eff.1 * (1 - val.sample.right$est.cond.eff.1), 2, mean)
            g.b.l.0 <- apply(val.w.0l * val.sample.left$est.cond.eff.0 * (1 - val.sample.left$est.cond.eff.0), 2, mean)
            g.b.r.0 <- apply(val.w.0r * val.sample.right$est.cond.eff.0 * (1 - val.sample.right$est.cond.eff.0), 2, mean)
            
          }
          
          left.p.1 <- as.matrix(data.frame(val.sample.left$est.cond.eff.1))
          left.p.0 <- as.matrix(data.frame(val.sample.left$est.cond.eff.0))
          right.p.1 <- as.matrix(data.frame(val.sample.right$est.cond.eff.1))
          right.p.0 <- as.matrix(data.frame(val.sample.right$est.cond.eff.0))
          
          n.l <- nrow(val.sample.left)
          n.r <- nrow(val.sample.right)
          
          # Calculate variance estimators
          var.1l <- g.b.l.1 %*% var.rb.l %*% g.b.l.1 + 1/n.l^2 * sum( (left.p.1 - mu.1l)^2 )
          var.1r <- g.b.r.1 %*% var.rb.r %*% g.b.r.1 + 1/n.r^2 * sum( (right.p.1 - mu.1r)^2 )
          var.0l <- g.b.l.0 %*% var.rb.l %*% g.b.l.0 + 1/n.l^2 * sum( (left.p.0 - mu.0l)^2 )
          var.0r <- g.b.r.0 %*% var.rb.r %*% g.b.r.0 + 1/n.r^2 * sum( (right.p.0 - mu.0r)^2 )
          
        }
        
        var.1l <- as.numeric(var.1l)
        var.1r <- as.numeric(var.1r)
        var.0l <- as.numeric(var.0l)
        var.0r <- as.numeric(var.0r)
        var    <- var.1l + var.0l + var.1r + var.0r
        
        if (var < 0) {
          next
        }
        
        # print((mu.1l - mu.0l) - (mu.1r - mu.0r))
        # print(var.l + var.r)
        (((mu.1l - mu.0l) - (mu.1r - mu.0r)) / sqrt(var))^2
        goodness.test <- goodness.test + (((mu.1l - mu.0l) - (mu.1r - mu.0r)) / sqrt(var))^2
        
        
        # print(var)
      } # End h
    } # End if loop
    # Calculating complexity value
    complex.val[m] = goodness.test - lambda.used * numb.int
  } # End m loop
  
  # Averaging over cross validation sets
  tree.final = tree.list[[which.max(complex.val)]]
  
  return(list(tree.final, complex.val))
}



######################################################################################################################
############################################ CV method 1, DR estimation ##############################################
######################################################################################################################
EstDr.CvMethod1 = function(data.used, tree.list, lambda.used, val.sample, type.var = "cont",
                           propsc.mod.loc, propsc.mthd = "GLM", propsc.form.true = NULL,
                           adj.mod.loc, adj.mthd = "GLM", adj.form.true = NULL) {
  # min.obs.mod = NULL
  
  # min.obs.mod: minimum number of observations for the model to be fitted
  # models when there are only treated/untreated unit in the subset
  val.sample.noy <- val.sample[, !colnames(val.sample) %in% c("Y")]
  whl.propsc <- gen.fullrank.ipw(df.noy           = val.sample.noy, 
                                 propsc.form.true = propsc.form.true)
  whl.propsc <- withWarnings(est.prop.sc(df.noy    = whl.propsc$df.noy.fullrank,
                                         method    = propsc.mthd,
                                         form.true = whl.propsc$propsc.form.true.updated))
  
  whl.g <- gen.fullrank.g(df            = val.sample,
                          adj.form.true = adj.form.true)
  whl.g <- withWarnings(est.cond.eff(df        = whl.g$df.fullrank,
                                     method    = adj.mthd,
                                     form.true = whl.g$adj.form.true.updated,
                                     type.var  = type.var))
  
  if (adj.mod.loc == "out") {
    val.sample$est.cond.eff.0 <- whl.g$value$pred.A.0
    val.sample$est.cond.eff.1 <- whl.g$value$pred.A.1
    
  }
  
  # need to add the propensity score back to val.sample 
  # since val.sample should not be affected when calculating cond.eff
  if (propsc.mod.loc == "out") {
    val.sample$prop.sc <- whl.propsc$value$prop.sc
  }
  
  # Calculate avg.trt.effct for when there is only treated/untreated unit in the subset
  # mu.1 <- mean(val.sample$A * (val.sample$Y - whl.g$pred.A.1) / whl.propsc$prop.sc + whl.g$pred.A.1)
  # mu.0 <- mean((1 - val.sample$A) * (val.sample$Y - whl.g$pred.A.0) / (1 - whl.propsc$prop.sc) + whl.g$pred.A.0)
  # 
  # avg.trt.effct <- mu.1 - mu.0
  
  # Storage space for the complexity value associated with each candidate tree
  complex.val = rep(NA, length(tree.list))
  
  # Looping through the candidate trees
  for (m in 1:length(complex.val)){
    # for (m in c(1)){
    # Finding the tree being evaluated
    # print(m)
    tree.used = tree.list[[m]] 
    
    # If only root node there is no internal node
    if(nrow(tree.used$frame) == 1){
      goodness.test = 0
      numb.int = 0
      
    } else {            # If at least one split
      is.leaf <- (tree.used$frame$var == "<leaf>")
      goodness.test = 0
      # tmp.goodness <- NULL
      # Calculate the goodness of the tree using the test
      # Finding the test data falling in each terminal node
      numb.int = sum(!is.leaf)
      # Finding all kids on terminal node 
      #kids.i <- order(as.numeric(rownames(tree.used$frame)))
      
      # Calculating split complexity for each internal node of branch
      right.ind  <- 1
      last.right <- NULL
      # w.last.right <- NULL
      
      for (h in 1:dim(tree.used$frame)[1]){
        # for (h in 3:109){
        # h <- h + 1
        
        # Finding observations in validation sample falling in that node
        if (tree.used$frame$var[h] == "<leaf>") {
          next
        } else if (tree.used$frame$n[h] == dim(data.used)[1]){
          val.sample.used <- val.sample
          # val.w.used      <- val.w
        } else {
          if (!is.null(last.left)){
            val.sample.used <- last.left[[1]]
            # val.w.used      <- w.last.left[[1]]
          } else {
            val.sample.used <- as.data.frame(last.right[[right.ind - 1]])
            # val.w.used      <- as.data.frame(w.last.right[[right.ind - 1]])
            if (length(last.right) > 2){
              last.right   <- last.right[1:(right.ind - 2)]
              # w.last.right <- w.last.right[1:(right.ind - 2)]
            } else if (length(last.right) == 2){
              last.right   <- list(last.right[[1]])
              # w.last.right <- list(w.last.right[[1]])
            }else {
              last.right <- NULL
              # w.last.right <- NULL
            }
            right.ind <- right.ind - 1
          }
        }
        
        # Predict propensity score inside node (need to be before identifying obs into right/left node since need propsc later)
        if (propsc.mod.loc == "node") { 
          # if (!propsc.mod.insplt) {
          # if (nrow(val.sample.used) >= min.obs.mod){ 
          # Fit propensity score model when it is possible to generate the full rank matrix.
          # only fit model when there are both treated and untreated units in the node
          # otherwise, the loop will be skipped later
          if (length(unique(val.sample.used$A)) > 1) { 
            val.sample.used.noy     <- val.sample.used[, !(colnames(val.sample.used) %in% c("Y", 
                                                                                            "prop.sc", 
                                                                                            "est.cond.eff.0", 
                                                                                            "est.cond.eff.1"))]
            tmp <- gen.fullrank.ipw(df.noy           = val.sample.used.noy, 
                                    propsc.form.true = propsc.form.true)
            tmp <- withWarnings(est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                                            method    = propsc.mthd,
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
                val.sample.used$prop.sc <- tmp$value$prop.sc
              } else {
                val.sample.used$prop.sc <- whl.propsc$value$prop.sc[rownames(val.sample) %in% rownames(val.sample.used)]
              }
              
            } else if (identical(sort(unique(tmp$value$prop.sc)), c(0.1, 0.9))) {
              val.sample.used$prop.sc <- whl.propsc$value$prop.sc[rownames(val.sample) %in% rownames(val.sample.used)]
              
            } else {
              val.sample.used$prop.sc <- tmp$value$prop.sc
            }
            
          } 
          # else if (length(unique(val.sample.used$A)) == 1) { # == 0 will be skipped later
          #   val.sample.used$prop.sc <- whl.propsc$prop.sc[rownames(val.sample) %in% rownames(val.sample.used)]
          # }
          
          # } # end if (!propsc.mod.insplt){
        } # end if (propsc.mod.loc == "node") { 
        
        # Adjustment model inside node (need to be before identifying obs into right/left node since need propsc later)
        if (adj.mod.loc == "node") { 
          # if (!adj.mod.insplt) {
          # if (nrow(val.sample.used) >= min.obs.mod) { 
          # Fit outcome model when it is possible to generate the full rank matrix.
          # only fit model when there are both treated and untreated units in the node
          # otherwise, the loop will be skipped later
          if (length(unique(val.sample.used$A)) > 1) { 
            
            val.sample.used.nocondeff <- val.sample.used[, !(colnames(val.sample.used) %in% c("est.cond.eff.0", 
                                                                                              "est.cond.eff.1", 
                                                                                              "prop.sc"))]
            tmp <- gen.fullrank.g(df            = val.sample.used.nocondeff,
                                  adj.form.true = adj.form.true)
            tmp <- withWarnings(est.cond.eff(df        = tmp$df.fullrank,
                                             method    = adj.mthd,
                                             form.true = tmp$adj.form.true.updated,
                                             type.var  = type.var))
            
            # df        = tmp$df.fullrank
            # method    = adj.mthd
            # form.true = tmp$adj.form.true.updated
            # type.var  = type.var
            
            # when the outcome model fitting produces warning
            # use outside model fitting results for conditional means
            if (!is.null(tmp$warnings)) {
              
              # when there is only rank-deficient warning, use the local prediction
              cond.warnings <- T
              for (warnings.i in 1:length(tmp$warnings)) {
                cond.warnings <- cond.warnings & (tmp$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
              }
              
              if (cond.warnings) {
                val.sample.used$est.cond.eff.0 <- tmp$value$pred.A.0
                val.sample.used$est.cond.eff.1 <- tmp$value$pred.A.1
              } else {
                val.sample.used$est.cond.eff.0 <- whl.g$value$pred.A.0[rownames(val.sample) %in% rownames(val.sample.used)]
                val.sample.used$est.cond.eff.1 <- whl.g$value$pred.A.1[rownames(val.sample) %in% rownames(val.sample.used)]
              }

            } else { # if (!is.null(tmp$warnings)) {
              val.sample.used$est.cond.eff.0 <- tmp$value$pred.A.0
              val.sample.used$est.cond.eff.1 <- tmp$value$pred.A.1
            }
            
          } 
          
          # else if (length(unique(val.sample.used$A)) == 1) { # == 0 will be skipped later
          #   val.sample.used$est.cond.eff.0 <- whl.g$pred.A.0[rownames(val.sample) %in% rownames(val.sample.used)]
          #   val.sample.used$est.cond.eff.1 <- whl.g$pred.A.1[rownames(val.sample) %in% rownames(val.sample.used)]
          # }
          # } # end if (!adj.mod.insplt){
        } # end if (!adj.mod.out){ 
        
        row.ind <- sum((tree.used$frame$var[1:h]) != "<leaf>")
        split.used <- tree.used$splits[row.ind, 4]
        var.used <- tree.used$frame$var[h]
        col.ind <- which(colnames(val.sample.used) == var.used)
        
        # Calculate goodness corresponding to split
        # Finding observations falling in right and left node
        # Need to figure out observations going to the left/right node
        # if ((split.used %% 1) == 0){   
        if (class(data.used[, as.character(var.used)]) == "factor"){   
          # Categorical covariate split
          lvls <- levels(val.sample[, col.ind])
          val.sample.left  <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used, ] == 1], ]
          val.sample.right <- val.sample.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          # val.w.used.left  <- val.w.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used, ] == 1], ]
          # val.w.used.right <- val.w.used[val.sample.used[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          
        } else{
          # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[row.ind, 2] > 0) {
            val.sample.left  <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
            # val.w.used.left  <- val.w.used[val.sample.used[,  col.ind] >= split.used, ]
            # val.w.used.right <- val.w.used[val.sample.used[,  col.ind] < split.used, ]
            
          } else {
            val.sample.left <- val.sample.used[val.sample.used[,  col.ind] < split.used, ]
            val.sample.right <- val.sample.used[val.sample.used[,  col.ind] >= split.used, ]
            # val.w.used.left  <- val.w.used[val.sample.used[,  col.ind] < split.used, ]
            # val.w.used.right <- val.w.used[val.sample.used[,  col.ind] >= split.used, ]
          }
          
        }
        
        if (tree.used$frame$var[h+1] != "<leaf>"){
          last.left   <- list(val.sample.left)
          # w.last.left <- list(val.w.used.left)
        } else{
          last.left   <- NULL
          # w.last.left <- NULL
        }
        
        which.right <- as.numeric(rownames(tree.used$frame)[h+1]) + 1
        if (tree.used$frame$var[as.numeric(rownames(tree.used$frame)) == which.right] != "<leaf>"){
          last.right[[right.ind]]   <- val.sample.right
          # w.last.right[[right.ind]] <- val.w.used.right
          right.ind <- right.ind + 1
        } 
        
        # Skip if the number of observations in node is fewer than prespecified number
        # Since there is no model fitted within node uptill here
        # if (!propsc.mod.out){
        #   if (!propsc.mod.insplt){
        #     if (nrow(val.sample.used) < min.obs.mod){
        #       next
        #     } 
        #   } else {# end if (!propsc.mod.insplt){ 
        #     
        #     if ((nrow(val.sample.left) < min.obs.mod) & (nrow(val.sample.right) < min.obs.mod)) {
        #       next                         # break if the number of observations in both son nodes is fewer than min.obs.mod
        #     } else if ((nrow(val.sample.left) == 1) | (nrow(val.sample.right) == 1)) {
        #       next
        #     }
        #     
        #   }
        # } # end if (!propsc.mod.out){
        
        # Skip if the number of observations in node is fewer than prespecified number
        # Since there is no model fitted within node uptill here 
        # but we still need to figure out what observations are in the right and left son
        # if (!adj.mod.out){
        #   if (!adj.mod.insplt){
        #     if (nrow(val.sample.used) < min.obs.mod){
        #       next
        #     }
        #   } else { # end if (!adj.mod.insplt){
        #     
        #     if ((nrow(val.sample.left) < min.obs.mod) & (nrow(val.sample.right) < min.obs.mod)) {
        #       next                         # break if the number of observations in both son nodes is fewer than min.obs.mod
        #     } else if ((nrow(val.sample.left) == 1) | (nrow(val.sample.right) == 1)) {
        #       next
        #     }
        #     
        #   }
        # } # end if (!adj.mod.out){
        
        # Skip if there is no obervation in left/right subgroup 
        # if (min(nrow(val.sample.left), nrow(val.sample.right)) == 0) {
        #   next
        # }
        
        # Skip if there is no treated/untreated in the left/right subgroup
        if (min(sum(val.sample.left$A == 1), sum(val.sample.left$A == 0), sum(val.sample.right$A == 1), sum(val.sample.right$A == 0)) == 0) {
          next
        }
        
        if (adj.mod.loc == "split") {
          
          # if (length(unique(val.sample.left$A)) > 1) {
          # will have both treated and untreated units in child node since skipped earlier
          # left
          val.sample.used.nocondeff.l <- val.sample.left[, !(colnames(val.sample.left) %in% c("est.cond.eff.0", 
                                                                                              "est.cond.eff.1", 
                                                                                              "prop.sc"))]
          tmp <- gen.fullrank.g(df            = val.sample.used.nocondeff.l,
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
              val.sample.left$est.cond.eff.0 <- tmp$value$pred.A.0
              val.sample.left$est.cond.eff.1 <- tmp$value$pred.A.1
            } else {
              val.sample.left$est.cond.eff.0 <- whl.g$value$pred.A.0[rownames(val.sample) %in% rownames(val.sample.left)]
              val.sample.left$est.cond.eff.1 <- whl.g$value$pred.A.1[rownames(val.sample) %in% rownames(val.sample.left)]
            }

          } else {
            val.sample.left$est.cond.eff.0 <- tmp$value$pred.A.0
            val.sample.left$est.cond.eff.1 <- tmp$value$pred.A.1
          }
          
          # } else { # length(unique(val.sample.left$A)) == 1, nrow equal 0 is skipped earlier
          #   val.sample.left$est.cond.eff.0 <- whl.g$pred.A.0[rownames(val.sample) %in% rownames(val.sample.left)]
          #   val.sample.left$est.cond.eff.1 <- whl.g$pred.A.1[rownames(val.sample) %in% rownames(val.sample.left)]
          # }
          
          # if (length(unique(val.sample.right$A)) > 1) {
          # will have both treated and untreated units in child node since skipped earlier
          # right
          val.sample.used.nocondeff.r <- val.sample.right[, !(colnames(val.sample.right) %in% c("est.cond.eff.0", 
                                                                                                "est.cond.eff.1", 
                                                                                                "prop.sc"))]
          tmp <- gen.fullrank.g(df            = val.sample.used.nocondeff.r,
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
              val.sample.right$est.cond.eff.0 <- tmp$value$pred.A.0
              val.sample.right$est.cond.eff.1 <- tmp$value$pred.A.1        
            } else {
              val.sample.right$est.cond.eff.0 <- whl.g$value$pred.A.0[rownames(val.sample) %in% rownames(val.sample.right)]
              val.sample.right$est.cond.eff.1 <- whl.g$value$pred.A.1[rownames(val.sample) %in% rownames(val.sample.right)]
            }

          } else { # if (!is.null(tmp$warnings)) {
            val.sample.right$est.cond.eff.0 <- tmp$value$pred.A.0
            val.sample.right$est.cond.eff.1 <- tmp$value$pred.A.1            
          }
          
          # } else { # nrow(val.sample.right) == 1, nrow equal 0 is skipped earlier
          #   val.sample.right$est.cond.eff.0 <- whl.g$pred.A.0[rownames(val.sample) %in% rownames(val.sample.right)]
          #   val.sample.right$est.cond.eff.1 <- whl.g$pred.A.1[rownames(val.sample) %in% rownames(val.sample.right)]
          # }
          
          
          val.sample.used$est.cond.eff.0[rownames(val.sample.used) %in% rownames(val.sample.left)]  <- val.sample.left$est.cond.eff.0
          val.sample.used$est.cond.eff.0[rownames(val.sample.used) %in% rownames(val.sample.right)] <- val.sample.right$est.cond.eff.0
          
          val.sample.used$est.cond.eff.1[rownames(val.sample.used) %in% rownames(val.sample.left)]  <- val.sample.left$est.cond.eff.1
          val.sample.used$est.cond.eff.1[rownames(val.sample.used) %in% rownames(val.sample.right)] <- val.sample.right$est.cond.eff.1
          
        }
        
        if (propsc.mod.loc == "split") { # propensity score model in split
          
          # Propensity score model
          # if (length(unique(val.sample.left$A)) > 1) {
          # will have both treated and untreated units in child node since skipped earlier
          # left 
          val.sample.noy.l  <- val.sample.left[, !(colnames(val.sample.left) %in% c("Y", 
                                                                                    "prop.sc", 
                                                                                    "est.cond.eff.0", 
                                                                                    "est.cond.eff.1"))]
          tmp.l <- gen.fullrank.ipw(df.noy           = val.sample.noy.l, 
                                    propsc.form.true = propsc.form.true)
          tmp.l <- withWarnings(est.prop.sc(df.noy    = tmp.l$df.noy.fullrank,
                                            method    = propsc.mthd,
                                            form.true = tmp.l$propsc.form.true.updated))
          
          # if 1) propensity score model fitting produces warning or 2) only produce thresholded propensity scores
          # use outside model fitting results
          if (!is.null(tmp.l$warnings)) {
            
            cond.warnings <- T
            for (warnings.i in 1:length(tmp.l$warnings)) {
              cond.warnings <- cond.warnings & (tmp.l$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
            }
            # when there is only rank-deficient warning, use the local prediction
            if (cond.warnings) {
              val.sample.left$prop.sc <- tmp.l$value$prop.sc
            } else {
              val.sample.left$prop.sc <- whl.propsc$value$prop.sc[rownames(val.sample) %in% rownames(val.sample.left)]
            }
            
          } else if (identical(sort(unique(tmp.l$value$prop.sc)), c(0.1, 0.9))) {
            val.sample.left$prop.sc <- whl.propsc$value$prop.sc[rownames(val.sample) %in% rownames(val.sample.left)]
            
          } else {
            val.sample.left$prop.sc <- tmp.l$value$prop.sc
          }
          
          # } else { # length(unique(val.sample.left$A)) == 1, nrow equal 0 is skipped earlier
          #   val.sample.left$prop.sc <- whl.propsc$prop.sc[rownames(val.sample) %in% rownames(val.sample.left)]
          # }
          
          # if (length(unique(val.sample.right$A)) > 1) {
          # will have both treated and untreated units in child node since skipped earlier
          # right 
          val.sample.noy.r <- val.sample.right[, !(colnames(val.sample.right) %in% c("Y", 
                                                                                     "prop.sc",
                                                                                     "est.cond.eff.0", 
                                                                                     "est.cond.eff.1"))]
          tmp.r <- gen.fullrank.ipw(df.noy           = val.sample.noy.r, 
                                    propsc.form.true = propsc.form.true)
          tmp.r <- withWarnings(est.prop.sc(df.noy    = tmp.r$df.noy.fullrank,
                                            method    = propsc.mthd,
                                            form.true = tmp.r$propsc.form.true.updated))
          
          # if 1) propensity score model fitting produces warning or 2) only produce thresholded propensity scores
          # use outside model fitting results
          if (!is.null(tmp.r$warnings)) {
            
            cond.warnings <- T
            for (warnings.i in 1:length(tmp.r$warnings)) {
              cond.warnings <- cond.warnings & (tmp.r$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
            }
            
            # when there is only rank-deficient warning, use the local prediction
            if (cond.warnings) {
              val.sample.right$prop.sc <- tmp.r$value$prop.sc
            } else {
              val.sample.right$prop.sc <- whl.propsc$value$prop.sc[rownames(val.sample) %in% rownames(val.sample.right)]
            }
            
          } else if (identical(sort(unique(tmp.r$value$prop.sc)), c(0.1, 0.9))) {
            val.sample.right$prop.sc <- whl.propsc$value$prop.sc[rownames(val.sample) %in% rownames(val.sample.right)]
          } else {
            val.sample.right$prop.sc <- tmp.r$value$prop.sc
          }
          
          # } else { # length(unique(val.sample.right$A)) == 1, nrow equal 0 is skipped earlier
          #   val.sample.right$prop.sc <- whl.propsc$prop.sc[rownames(val.sample) %in% rownames(val.sample.right)]
          # }
          
          val.sample.used$prop.sc[rownames(val.sample.used) %in% rownames(val.sample.left)] <- val.sample.left$prop.sc
          val.sample.used$prop.sc[rownames(val.sample.used) %in% rownames(val.sample.right)] <- val.sample.right$prop.sc
          
        } 
        
        mu.1l <- mean(val.sample.left$A * (val.sample.left$Y - val.sample.left$est.cond.eff.1) / 
                        val.sample.left$prop.sc + val.sample.left$est.cond.eff.1)
        mu.0l <- mean((1 - val.sample.left$A) * (val.sample.left$Y - val.sample.left$est.cond.eff.0) / 
                        (1 - val.sample.left$prop.sc) + val.sample.left$est.cond.eff.0)
        mu.1r <- mean(val.sample.right$A * (val.sample.right$Y - val.sample.right$est.cond.eff.1) / 
                        val.sample.right$prop.sc + val.sample.right$est.cond.eff.1)
        mu.0r <- mean((1 - val.sample.right$A) * (val.sample.right$Y - val.sample.right$est.cond.eff.0) / 
                        (1 - val.sample.right$prop.sc) + val.sample.right$est.cond.eff.0)
        
        t.dr <- (mu.1l - mu.0l) - (mu.1r - mu.0r)
        
        n.l <- nrow(val.sample.left)
        n.r <- nrow(val.sample.right)
        p.l <- n.l / (n.l + n.r)
        p.r <- n.r / (n.l + n.r)
        
        # Calculate variance estimators
        I_i <- (rownames(val.sample.used) %in% rownames(val.sample.left)) / p.l * 
          ((val.sample.used$A * (val.sample.used$Y - val.sample.used$est.cond.eff.1) / 
              val.sample.used$prop.sc + val.sample.used$est.cond.eff.1) -
             ((1 - val.sample.used$A) * (val.sample.used$Y - val.sample.used$est.cond.eff.0) / 
                (1 - val.sample.used$prop.sc) + val.sample.used$est.cond.eff.0)) - 
          (rownames(val.sample.used) %in% rownames(val.sample.right)) / p.r * 
          ((val.sample.used$A * (val.sample.used$Y - val.sample.used$est.cond.eff.1) / 
              val.sample.used$prop.sc + val.sample.used$est.cond.eff.1) - 
             ((1 - val.sample.used$A) * (val.sample.used$Y - val.sample.used$est.cond.eff.0) / 
                (1 - val.sample.used$prop.sc) + val.sample.used$est.cond.eff.0)) - 
          t.dr
        
        
        # print((mu.1l - mu.0l) - (mu.1r - mu.0r))
        # print(var.l + var.r)
        
        # variance smaller if both are glm models
        if ((propsc.mthd == "GLM") & (adj.mthd == "GLM")) {
          var <- (mean(I_i^2) - (p.r * (mu.1l - mu.0l) + p.l * (mu.1r - mu.0r))^2 / (p.l * p.r)) / (n.l + n.r)
          
          # Skip this loop if there is negative estimated variance
          if (var < 0) {
            next
          } 
        } else {
          var <- mean(I_i^2) / (n.l + n.r)
        }
        
        # tmp.goodness <- c(tmp.goodness, (t.dr / sqrt(var))^2)
        # (t.dr / sqrt(var))^2
        goodness.test <- goodness.test + (t.dr / sqrt(var))^2
        
        if (propsc.mod.loc != "out") {
          val.sample.used  <- val.sample.used[, !(colnames(val.sample.used) %in% c("prop.sc"))]
          val.sample.left  <- val.sample.left[, !(colnames(val.sample.left) %in% c("prop.sc"))]
          val.sample.right <- val.sample.right[, !(colnames(val.sample.right) %in% c("prop.sc"))]
        }
        
        if (adj.mod.loc != "out") {
          val.sample.used  <- val.sample.used[, !(colnames(val.sample.used) %in% c("est.cond.eff.0", "est.cond.eff.1"))]
          val.sample.left  <- val.sample.left[, !(colnames(val.sample.left) %in% c("est.cond.eff.0", "est.cond.eff.1"))]
          val.sample.right <- val.sample.right[, !(colnames(val.sample.right) %in% c("est.cond.eff.0", "est.cond.eff.1"))]
        }
        
        # print(var)
      } # End h
    } # End if loop
    # Calculating complexity value
    complex.val[m] = goodness.test - lambda.used * numb.int
  } # End m loop
  
  # Averaging over cross validation sets
  tree.final = tree.list[[which.max(complex.val)]]
  
  return(list(tree.final, complex.val))
}

