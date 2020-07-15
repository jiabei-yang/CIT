######################################################################################################################
############################################ CV method 1, DR estimation ##############################################
######################################################################################################################
EstDr.CvMethod1.MltLamda = function(data.used, tree.list, lambda.used, val.sample, type.var = "cont",
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
  complex.val = matrix(NA, ncol = length(tree.list), nrow = length(lambda.used))
  
  # Looping through the candidate trees
  for (m in 1:ncol(complex.val)){
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
    
    for (lambda.i in 1:length(lambda.used)){
      complex.val[lambda.i, m] = goodness.test - lambda.used[lambda.i] * numb.int
    }
    
  } # End m loop
  
  tree.final <- list()
  # Averaging over cross validation sets
  for (lambda.i in 1:length(lambda.used)){
    tree.final[[lambda.i]] = tree.list[[which.max(complex.val[lambda.i, ])]]
  }
  
  return(list(tree.final, complex.val))
}
