######################################################################################################################
################################################# CV method 1, IPW ###################################################
######################################################################################################################
EstIpw.CvMethod1 = function(data.used, tree.list, lambda.used, val.sample, type.var = "cont",
                            propsc.mod.out = T, propsc.mthd = "GLM", propsc.form.true = NULL, val.w = NULL,
                            propsc.mod.insplt = NULL, min.obs.mod = NULL){
  
  # type.var does not matter for IPW, the parameter can be omitted, but the previous code needs to be modified
  # min.obs.mod: minimum number of observations for the model to be fitted
  
  # if (is.null(val.w)){                          # True model or not, assign value to w
  #   val.w <- cbind(1, val.sample[, 3:ncol(val.sample)])
  # }
  
  if (propsc.mod.out){
    
    val.sample.noy <- val.sample[, !colnames(val.sample) %in% c("Y")]
    
    tmp <- gen.fullrank.ipw(df.noy           = val.sample.noy, 
                            propsc.form.true = propsc.form.true)
    tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                       method    = propsc.mthd,
                       form.true = tmp$propsc.form.true.updated)
    
    val.sample$prop.sc <- tmp$prop.sc
    val.w <- tmp$w
    
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
        # for (h in 1:10){
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
        if (!propsc.mod.out){ 
          if (!propsc.mod.insplt){
            if (nrow(val.sample.used) >= min.obs.mod){ 
              val.sample.used.noy     <- val.sample.used[, !(colnames(val.sample.used) %in% c("Y", "prop.sc"))]
              
              tmp <- gen.fullrank.ipw(df.noy           = val.sample.used.noy, 
                                      propsc.form.true = propsc.form.true)
              tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                                 method    = propsc.mthd,
                                 form.true = tmp$propsc.form.true.updated)
              val.sample.used$prop.sc <- tmp$prop.sc
              val.w.used              <- tmp$w
              
            }
          } # end if (!propsc.mod.insplt){
        } # end if (!propsc.mod.out){ 
        
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
        
        # Skip if the number of observations in node is fewer than prespecified number
        # Since there is no model fitted within node uptill here
        if (!propsc.mod.out){
          if (!propsc.mod.insplt){
            if (nrow(val.sample.used) < min.obs.mod){
              next
            }
          } # end if (!propsc.mod.insplt){
        } # end if (!propsc.mod.out){
        
        # Skip if the number of observations in any split is smaller than the dimension of covariates to avoid variance computation error
        # if (!is.null(val.w.used)){
        #   if ((nrow(val.sample.left) < ncol(val.w.used)) | (nrow(val.sample.right) < ncol(val.w.used))){
        #     next             
        #   }
        # }
        
        # Skip if there is no treated/untreated in the left/right subgroup
        if (min(sum(val.sample.left$A == 1), sum(val.sample.left$A == 0), sum(val.sample.right$A == 1), sum(val.sample.right$A == 0)) == 0){
          next
        }
        
        if (!is.null(val.sample.used$prop.sc)) {
          
          mu.1l <- mean(val.sample.left$Y * val.sample.left$A / val.sample.left$prop.sc)
          mu.0l <- mean(val.sample.left$Y * (1 - val.sample.left$A)/ (1 - val.sample.left$prop.sc))
          mu.1r <- mean(val.sample.right$Y * val.sample.right$A / val.sample.right$prop.sc)
          mu.0r <- mean(val.sample.right$Y * (1 - val.sample.right$A)/ (1 - val.sample.right$prop.sc))
          
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
                       ((mu.1l - mu.0l) - (mu.1r - mu.0r)) - 
                       as.numeric((val.sample.used$A - val.sample.used$prop.sc) * (h.l - h.r) %*% solve(e.bb) %*% t(val.w.used)))
          
          # When the number of observations in val.w.used very sumall, call problem when solve(e.bb)
          # will skip this loop, so the addition to goodness.test will be 0
          if ("try-error" %in% class(I.i)) {
            next
          }
          
          var <- (mean(I.i^2) - 
                    (n.r / n.used * (mu.1l - mu.0l) + n.l / n.used * (mu.1r - mu.0r))^2 / 
                    (n.l * n.r / n.used^2)) / n.used
          
        } else { # only when propsc.mod.insplt is TRUE
          
          if ((nrow(val.sample.left) < min.obs.mod) & (nrow(val.sample.right) < min.obs.mod)) {
            next                         # break if the number of observations in the son nodes is fewer than min.obs.mod
          } else if ((nrow(val.sample.left) == 1) | (nrow(val.sample.right) == 1)) {
            next
          }
          
          val.sample.noy.l  <- val.sample.left[, !(colnames(val.sample.left) %in% c("Y", "prop.sc"))]
          val.sample.noy.r <- val.sample.right[, !(colnames(val.sample.right) %in% c("Y", "prop.sc"))]
          
          tmp.l <- gen.fullrank.ipw(df.noy           = val.sample.noy.l, 
                                    propsc.form.true = propsc.form.true)
          tmp.r <- gen.fullrank.ipw(df.noy           = val.sample.noy.r, 
                                    propsc.form.true = propsc.form.true)
          
          tmp.l <- est.prop.sc(df.noy    = tmp.l$df.noy.fullrank,
                               method    = propsc.mthd,
                               form.true = tmp.l$propsc.form.true.updated)
          tmp.r <- est.prop.sc(df.noy    = tmp.r$df.noy.fullrank,
                               method    = propsc.mthd,
                               form.true = tmp.r$propsc.form.true.updated)
          
          prop.sc.l <- tmp.l$prop.sc
          prop.sc.r <- tmp.r$prop.sc
          val.w.used.left  <- tmp.l$w
          val.w.used.right <- tmp.r$w
          
          # Skip if there are fewer observations than the number of columns in w
          # if ((nrow(val.sample.left) < ncol(val.w.used.left)) | (nrow(val.sample.right) < ncol(val.w.used.right))){
          #   next             
          # }
          
          val.sample.used$prop.sc[rownames(val.sample.used) %in% rownames(val.sample.left)] <- prop.sc.l
          val.sample.used$prop.sc[rownames(val.sample.used) %in% rownames(val.sample.right)] <- prop.sc.r
          
          prop.sc <- val.sample.used$prop.sc
          val.sample.used <- val.sample.used[, colnames(val.sample.used) != "prop.sc"]
          
          # Numerator in the splitting statistics
          mu.1l <- mean(val.sample.left$Y * val.sample.left$A / prop.sc.l)
          mu.0l <- mean(val.sample.left$Y * (1 - val.sample.left$A)/ (1 - prop.sc.l))
          mu.1r <- mean(val.sample.right$Y * val.sample.right$A / prop.sc.r)
          mu.0r <- mean(val.sample.right$Y * (1 - val.sample.right$A)/ (1 - prop.sc.r))
          
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
          
          term.l <- as.numeric((val.sample.left$A - prop.sc.l) * h.l %*% solve(e.bb.l) %*% t(val.w.used.left))
          term.r <- as.numeric((val.sample.right$A - prop.sc.r) * h.r %*% solve(e.bb.r) %*% t(val.w.used.right))
          
          # if (("try-error" %in% class(term.l)) | ("try-error" %in% class(term.r))) {
          #   next
          # }
          
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
            ((mu.1l - mu.0l) - (mu.1r - mu.0r)) - var.term
          
          var <- (mean(I.i^2) - 
                    (n.r / n.used * (mu.1l - mu.0l) + n.l / n.used * (mu.1r - mu.0r))^2 / 
                    (n.l * n.r / n.used^2)) / n.used
        }
        
        
        # print((mu.1l - mu.0l) - (mu.1r - mu.0r))
        # print(var.l + var.r)
        if (var < 0){
          next
        } else{
          goodness.test <- goodness.test + (((mu.1l - mu.0l) - (mu.1r - mu.0r)) / sqrt(var))^2
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
                          adj.mod.out = T, adj.mthd = "GLM", adj.form.true = NULL, val.w = NULL,
                          adj.mod.insplt = NULL, min.obs.mod = NULL){
  
  # min.obs.mod: minimum number of observations for the model to be fitted
  
  # if (is.null(val.w)){                          # True model or not, assign value to w
  #   val.w <- cbind(1, val.sample$A, val.sample[, 3:ncol(val.sample)], val.sample$A * val.sample[, 3:ncol(val.sample)])
  #   colnames(val.w) <- c("(Intercept)", "A", colnames(val.sample)[3:ncol(val.sample)], paste("A:", colnames(val.sample)[3:ncol(val.sample)], sep = ""))
  # }
  
  if (adj.mod.out){
    
    tmp <- est.cond.eff(df        = val.sample,
                        method    = adj.mthd,
                        form.true = adj.form.true,
                        type.var  = type.var)
    
    val.sample$est.cond.eff.0 <- tmp$pred.A.0
    val.sample$est.cond.eff.1 <- tmp$pred.A.1
    var.rb                    <- tmp$var.rb
    val.w                     <- tmp$w
    
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
        # for (h in 1:10){
        
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
        if (!adj.mod.out) { 
          if (!adj.mod.insplt) {
            if (nrow(val.sample.used) >= min.obs.mod) { 
              
              val.sample.used.nocondeff <- val.sample.used[, !(colnames(val.sample.used) %in% c("est.cond.eff.0", "est.cond.eff.1"))]
              tmp <- gen.fullrank.g(df            = val.sample.used.nocondeff,
                                    adj.form.true = adj.form.true)
              
              tmp <- est.cond.eff(df        = tmp$df.fullrank,
                                  method    = adj.mthd,
                                  form.true = tmp$adj.form.true.updated,
                                  type.var  = type.var)
              
              val.sample.used$est.cond.eff.0 <- tmp$pred.A.0
              val.sample.used$est.cond.eff.1 <- tmp$pred.A.1
              var.rb                         <- tmp$var.rb
              val.w.used                     <- tmp$w
              
            }
          } # end if (!propsc.mod.insplt){
        } # end if (!propsc.mod.out){ 
        
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
        
        # Skip if the number of observations in node is fewer than prespecified number
        # Since there is no model fitted within node uptill here 
        # but we still need to figure out what observations are in the right and left son
        if (!adj.mod.out){
          if (!adj.mod.insplt) {
            if (nrow(val.sample.used) < min.obs.mod){
              next
            }
          } else { # if (!adj.mod.insplt){ this else is for models fitted in split
            if ((nrow(val.sample.left) < min.obs.mod) & (nrow(val.sample.right) < min.obs.mod)) {
              next                         # break if the number of observations in both son nodes is fewer than min.obs.mod
            } else if ((nrow(val.sample.left) == 1) | (nrow(val.sample.right) == 1)) {
              next
            }
          } # end if (!adj.mod.insplt){
        } # end if (!adj.mod.out){
        
        # Skip if the number of observations in any split is smaller than the dimension of covariates to avoid variance computation error
        # if (!is.null(val.w.used)) {
        #   if ((nrow(val.sample.left) < ncol(val.w.used)) | (nrow(val.sample.right) < ncol(val.w.used))){
        #     next             
        #   }
        # }
        
        # Skip if there is no treated/untreated in the left/right subgroup
        if (min(sum(val.sample.left$A == 1), sum(val.sample.left$A == 0), sum(val.sample.right$A == 1), sum(val.sample.right$A == 0)) == 0){
          next
        }
        
        if (exists("var.rb")) {
          
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
          
          var.1l <- as.numeric(var.1l)
          var.1r <- as.numeric(var.1r)
          var.0l <- as.numeric(var.0l)
          var.0r <- as.numeric(var.0r)
          
          if (!is.null(adj.mod.insplt)){
            rm(var.rb)
          }
          
        } else { # if(exists("var.rb")) when adj.mod.insplt = T
          
          val.sample.left.nocondeff  <- val.sample.left[, !(colnames(val.sample.left) %in% c("est.cond.eff.0", "est.cond.eff.1"))]
          val.sample.right.nocondeff <- val.sample.right[, !(colnames(val.sample.right) %in% c("est.cond.eff.0", "est.cond.eff.1"))]
          
          tmp.l <- gen.fullrank.g(df            = val.sample.left.nocondeff,
                                  adj.form.true = adj.form.true)
          
          tmp.r <- gen.fullrank.g(df            = val.sample.right.nocondeff,
                                  adj.form.true = adj.form.true)
          
          tmp.l <- est.cond.eff(df        = tmp.l$df.fullrank,
                                method    = adj.mthd,
                                form.true = tmp.l$adj.form.true.updated,
                                type.var  = type.var)
          tmp.r <- est.cond.eff(df        = tmp.r$df.fullrank,
                                method    = adj.mthd,
                                form.true = tmp.r$adj.form.true.updated,
                                type.var  = type.var)          
          
          val.sample.left$est.cond.eff.0 <- tmp.l$pred.A.0
          val.sample.left$est.cond.eff.1 <- tmp.l$pred.A.1
          var.rb.l                       <- tmp.l$var.rb
          val.w.used.left                <- tmp.l$w
          
          val.sample.right$est.cond.eff.0 <- tmp.r$pred.A.0
          val.sample.right$est.cond.eff.1 <- tmp.r$pred.A.1
          var.rb.r                        <- tmp.r$var.rb
          val.w.used.right                <- tmp.r$w
          
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
          
          var.1l <- as.numeric(var.1l)
          var.1r <- as.numeric(var.1r)
          var.0l <- as.numeric(var.0l)
          var.0r <- as.numeric(var.0r)
          
        }
        
        # print((mu.1l - mu.0l) - (mu.1r - mu.0r))
        # print(var.l + var.r)
        goodness.test <- goodness.test + (((mu.1l - mu.0l) - (mu.1r - mu.0r)) / sqrt(var.1l + var.0l + var.1r + var.0r))^2
        
        
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
                           propsc.mod.out = T, propsc.mthd = "GLM", 
                           propsc.form.true = NULL, propsc.mod.insplt = NULL, 
                           adj.mod.out = T, adj.mthd = "GLM", adj.form.true = NULL, 
                           adj.mod.insplt = NULL, min.obs.mod = NULL){
  
  # min.obs.mod: minimum number of observations for the model to be fitted
  if (propsc.mod.out) {
    
    val.sample.noy <- val.sample[, !colnames(val.sample) %in% c("Y")]
    tmp.propsc <- gen.fullrank.ipw(df.noy           = val.sample.noy, 
                                   propsc.form.true = propsc.form.true)
    tmp.propsc <- est.prop.sc(df.noy    = tmp.propsc$df.noy.fullrank,
                              method    = propsc.mthd,
                              form.true = tmp.propsc$propsc.form.true.updated)
    
  }
  
  if (adj.mod.out) {
    
    tmp <- gen.fullrank.g(df            = val.sample,
                          adj.form.true = adj.form.true)
    tmp <- est.cond.eff(df        = tmp$df.fullrank,
                        method    = adj.mthd,
                        form.true = tmp$adj.form.true.updated,
                        type.var  = type.var)
    
    val.sample$est.cond.eff.0 <- tmp$pred.A.0
    val.sample$est.cond.eff.1 <- tmp$pred.A.1
    
  }
  
  # need to add the propensity score back to val.sample 
  # since val.sample should not be affected when calculating cond.eff
  if (propsc.mod.out) {
    val.sample$prop.sc <- tmp.propsc$prop.sc
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
      # w.last.right <- NULL
      
      for (h in 1:dim(tree.used$frame)[1]){
        # for (h in 1:25){
        # Finding observations in validation sample falling in that node
        if (tree.used$frame$var[h] == "<leaf>"){
          next
        } else if (tree.used$frame$n[h] == dim(data.used)[1]){
          val.sample.used <- val.sample
          # val.w.used      <- val.w
        } else{
          
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
        if (!propsc.mod.out) { 
          if (!propsc.mod.insplt) {
            if (nrow(val.sample.used) >= min.obs.mod){ 
              
              val.sample.used.noy     <- val.sample.used[, !(colnames(val.sample.used) %in% c("Y", 
                                                                                              "prop.sc", 
                                                                                              "est.cond.eff.0", 
                                                                                              "est.cond.eff.1"))]
              tmp <- gen.fullrank.ipw(df.noy           = val.sample.used.noy, 
                                      propsc.form.true = propsc.form.true)
              val.sample.used$prop.sc <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                                                     method    = propsc.mthd,
                                                     form.true = tmp$propsc.form.true.updated)$prop.sc
            }
          } # end if (!propsc.mod.insplt){
        } # end if (!propsc.mod.out){ 
        
        # Adjustment model inside node (need to be before identifying obs into right/left node since need propsc later)
        if (!adj.mod.out) { 
          if (!adj.mod.insplt) {
            if (nrow(val.sample.used) >= min.obs.mod) { 
              
              val.sample.used.nocondeff <- val.sample.used[, !(colnames(val.sample.used) %in% c("est.cond.eff.0", 
                                                                                                "est.cond.eff.1", 
                                                                                                "prop.sc"))]
              tmp <- gen.fullrank.g(df            = val.sample.used.nocondeff,
                                    adj.form.true = adj.form.true)
              tmp <- est.cond.eff(df        = tmp$df.fullrank,
                                  method    = adj.mthd,
                                  form.true = tmp$adj.form.true.updated,
                                  type.var  = type.var)
              
              val.sample.used$est.cond.eff.0 <- tmp$pred.A.0
              val.sample.used$est.cond.eff.1 <- tmp$pred.A.1
              
            }
          } # end if (!propsc.mod.insplt){
        } # end if (!propsc.mod.out){ 
        
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
        if (!propsc.mod.out){
          if (!propsc.mod.insplt){
            if (nrow(val.sample.used) < min.obs.mod){
              next
            } 
          } else {# end if (!propsc.mod.insplt){ 
            
            if ((nrow(val.sample.left) < min.obs.mod) & (nrow(val.sample.right) < min.obs.mod)) {
              next                         # break if the number of observations in both son nodes is fewer than min.obs.mod
            } else if ((nrow(val.sample.left) == 1) | (nrow(val.sample.right) == 1)) {
              next
            }
            
          }
        } # end if (!propsc.mod.out){
        
        # Skip if the number of observations in node is fewer than prespecified number
        # Since there is no model fitted within node uptill here 
        # but we still need to figure out what observations are in the right and left son
        if (!adj.mod.out){
          if (!adj.mod.insplt){
            if (nrow(val.sample.used) < min.obs.mod){
              next
            }
          } else { # end if (!adj.mod.insplt){
            
            if ((nrow(val.sample.left) < min.obs.mod) & (nrow(val.sample.right) < min.obs.mod)) {
              next                         # break if the number of observations in both son nodes is fewer than min.obs.mod
            } else if ((nrow(val.sample.left) == 1) | (nrow(val.sample.right) == 1)) {
              next
            }
            
          }
        } # end if (!adj.mod.out){
        
        # Skip if there is no treated/untreated in the left/right subgroup
        if (min(sum(val.sample.left$A == 1), sum(val.sample.left$A == 0), sum(val.sample.right$A == 1), sum(val.sample.right$A == 0)) == 0){
          next
        }
        
        if (is.null(val.sample.left$est.cond.eff.1)) {
          
          val.sample.used.nocondeff.l <- val.sample.left[, !(colnames(val.sample.left) %in% c("est.cond.eff.0", 
                                                                                              "est.cond.eff.1", 
                                                                                              "prop.sc"))]
          tmp <- gen.fullrank.g(df            = val.sample.used.nocondeff.l,
                                adj.form.true = adj.form.true)
          tmp <- est.cond.eff(df        = tmp$df.fullrank,
                              method    = adj.mthd,
                              form.true = tmp$adj.form.true.updated,
                              type.var  = type.var)
          
          val.sample.left$est.cond.eff.0 <- tmp$pred.A.0
          val.sample.left$est.cond.eff.1 <- tmp$pred.A.1
          
          val.sample.used.nocondeff.r <- val.sample.right[, !(colnames(val.sample.right) %in% c("est.cond.eff.0", 
                                                                                                "est.cond.eff.1", 
                                                                                                "prop.sc"))]
          tmp <- gen.fullrank.g(df            = val.sample.used.nocondeff.r,
                                adj.form.true = adj.form.true)
          tmp <- est.cond.eff(df        = tmp$df.fullrank,
                              method    = adj.mthd,
                              form.true = tmp$adj.form.true.updated,
                              type.var  = type.var)
          
          val.sample.right$est.cond.eff.0 <- tmp$pred.A.0
          val.sample.right$est.cond.eff.1 <- tmp$pred.A.1
          
          val.sample.used$est.cond.eff.0[rownames(val.sample.used) %in% rownames(val.sample.left)]  <- val.sample.left$est.cond.eff.0
          val.sample.used$est.cond.eff.0[rownames(val.sample.used) %in% rownames(val.sample.right)] <- val.sample.right$est.cond.eff.0

          val.sample.used$est.cond.eff.1[rownames(val.sample.used) %in% rownames(val.sample.left)]  <- val.sample.left$est.cond.eff.1
          val.sample.used$est.cond.eff.1[rownames(val.sample.used) %in% rownames(val.sample.right)] <- val.sample.right$est.cond.eff.1
                    
        }
        
        if (is.null(val.sample.left$prop.sc)) { # Both models in split
          
          # Propensity score model
          val.sample.noy.l  <- val.sample.left[, !(colnames(val.sample.left) %in% c("Y", 
                                                                                    "prop.sc", 
                                                                                    "est.cond.eff.0", 
                                                                                    "est.cond.eff.1"))]
          val.sample.noy.r <- val.sample.right[, !(colnames(val.sample.right) %in% c("Y", 
                                                                                     "prop.sc",
                                                                                     "est.cond.eff.0", 
                                                                                     "est.cond.eff.1"))]
          
          tmp.l <- gen.fullrank.ipw(df.noy           = val.sample.noy.l, 
                                    propsc.form.true = propsc.form.true)
          tmp.r <- gen.fullrank.ipw(df.noy           = val.sample.noy.r, 
                                    propsc.form.true = propsc.form.true)
          
          tmp.l <- est.prop.sc(df.noy    = tmp.l$df.noy.fullrank,
                               method    = propsc.mthd,
                               form.true = tmp.l$propsc.form.true.updated)
          tmp.r <- est.prop.sc(df.noy    = tmp.r$df.noy.fullrank,
                               method    = propsc.mthd,
                               form.true = tmp.r$propsc.form.true.updated)
          
          prop.sc.l <- tmp.l$prop.sc
          prop.sc.r <- tmp.r$prop.sc
          
          val.sample.left$prop.sc <- prop.sc.l
          val.sample.right$prop.sc <- prop.sc.r
          val.sample.used$prop.sc[rownames(val.sample.used) %in% rownames(val.sample.left)] <- prop.sc.l
          val.sample.used$prop.sc[rownames(val.sample.used) %in% rownames(val.sample.right)] <- prop.sc.r
          
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
        var <- (mean(I_i^2) - (p.r * (mu.1l - mu.0l) + p.l * (mu.1r - mu.0r))^2 / (p.l * p.r)) / (n.l + n.r)
        
        # Skip this loop if there is negative estimated variance
        if (var < 0) {
          next
        }
        
        goodness.test <- goodness.test + (t.dr / sqrt(var))^2
        
        if (!propsc.mod.out) {
          val.sample.used  <- val.sample.used[, !(colnames(val.sample.used) %in% c("prop.sc"))]
          val.sample.left  <- val.sample.left[, !(colnames(val.sample.left) %in% c("prop.sc"))]
          val.sample.right <- val.sample.right[, !(colnames(val.sample.right) %in% c("prop.sc"))]
        }
        
        if (!adj.mod.out) {
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



