######################################################################################################################
################################################# CV method 2, IPW ###################################################
######################################################################################################################
EstIpw.CvMethod2 = function(data.used, tree.list, type.var = "cont", seed = NULL, n.cv = 5,
                            propsc.mod.out = T, propsc.mthd = "GLM", propsc.form.true = NULL, min.obs.mod = NULL){
  
  if (!is.null(seed)){
    set.seed(seed)
  }
  
  cross.val.ind = sample(cut(seq(1,nrow(data.used)), breaks=n.cv, labels=FALSE), 
                         nrow(data.used))
  cv.err = matrix(NA, ncol = n.cv, nrow = length(tree.list))
  
  for (l in 1:n.cv){
    
    test.data = data.used[which(cross.val.ind == l), ]
    train.data = data.used[which(cross.val.ind != l), ]
    
    # Fit the random forest procedure on the test data to find the gold standard
    train.data.noy <- train.data[, colnames(train.data) != "Y"]
    train.data.noy <- train.data.noy %>%
      mutate(prop.sc = est.prop.sc(train.data.noy, method = "GLM")$prop.sc) %>%
      mutate(wt      = ifelse(A == 1, prop.sc, 1 - prop.sc))
    
    rand.for.train = rfsrc(Y~., data = train.data, case.wt = 1 / train.data.noy$wt)
    test.A.1 = test.data
    test.A.1$A = 1
    pred.A.1 = predict(rand.for.train, newdata = test.A.1)
    test.A.0 = test.data
    test.A.0$A = 0
    pred.A.0 = predict(rand.for.train, newdata = test.A.0)
    treat.eff.test = pred.A.1$predicted - pred.A.0$predicted
    
    # Calculate prop.sc if prop.sc model fitted outside
    train.data.noy <- train.data[, colnames(train.data) != "Y"]
    prop.sc.train <- est.prop.sc(df.noy    = train.data.noy, 
                                 method    = propsc.mthd,
                                 form.true = propsc.form.true)$prop.sc

    # Calculating g(h) for the cross.validated tree
    # Looping through the candidate trees
    for(m in 1:length(tree.list)){

      tree.used = tree.list[[m]]
      
      if (nrow(tree.used$frame) >3) {
        
        pred.train = predict(tree.used, newdata = train.data)
        pred.tree = predict(tree.used, newdata = test.data)
        for(k in 1:length(unique(pred.train))){
          # Calculating the treatment effect estimator for each node using only the 
          # training data
          data.node = train.data[which(pred.train == unique(pred.train)[k]), ]
          test.index = which(pred.tree == unique(pred.train)[k])
          
          if (!propsc.mod.out) {
            
            if (nrow(data.node) >= min.obs.mod){  # When prop.mod.out = T, min.obs.mod = NULL
              data.node.noy <- data.node[, colnames(data.node) != "Y"]
              
              tmp <- gen.fullrank.ipw(df.noy           = data.node.noy, 
                                      propsc.form.true = propsc.form.true)
              
              tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                                 method    = propsc.mthd,
                                 form.true = tmp$propsc.form.true.updated)
              
              prop.sc.node <- tmp$prop.sc
              
            } else {
              prop.sc.node <- prop.sc.train[which(pred.train == unique(pred.train)[k])]
            }
            
          } else { # model prop.sc.out / the number of observations in the node small
            prop.sc.node <- prop.sc.train[which(pred.train == unique(pred.train)[k])]
          }
          
          mu.1 <- mean(data.node$Y * data.node$A/ prop.sc.node)
          mu.0 <- mean(data.node$Y * (1 - data.node$A)/ (1 - prop.sc.node))
          
          pred.tree[test.index] <- mu.1 - mu.0
          
        }
        
      } else if (nrow(tree.used$frame) == 3){
        
        # If there is one or zero splits there is a weird memory error so need to do manually
        pred.tree = rep(NA, nrow(test.data))
        split.used = tree.used$splits[, 4]
        var.used = tree.used$frame$var[1]
        col.ind <- which(colnames(train.data) == var.used)
        
        # Need to figure out observations going to the left/right node
        if ((split.used %% 1) == 0){   
          
          # Categorical covariate split
          lvls <- levels(train.data[, col.ind])
          data.node.l <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1], ]
          data.node.r <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          
          if (propsc.mod.out){ # When prop.mod.out = T, min.obs.mod = NULL
            prop.sc.node.l <- prop.sc.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]]
            prop.sc.node.r <- prop.sc.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]]
          } else { 
            if (nrow(data.node.l) < min.obs.mod) {
              prop.sc.node.l <- prop.sc.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]]
            }
            
            if (nrow(data.node.r) < min.obs.mod) {
              prop.sc.node.r <- prop.sc.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]]
            }
          }
          
        } else{
          # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            data.node.l <- train.data[train.data[,  col.ind] >= split.used, ]
            data.node.r <- train.data[train.data[,  col.ind] < split.used, ]
            
            if (propsc.mod.out){
              prop.sc.node.l <- prop.sc.train[train.data[, col.ind] >= split.used]
              prop.sc.node.r <- prop.sc.train[train.data[, col.ind] < split.used]
            } else {
              if (nrow(data.node.l) < min.obs.mod) {
                prop.sc.node.l <- prop.sc.train[train.data[, col.ind] >= split.used]
              }
              
              if (nrow(data.node.r) < min.obs.mod){
                prop.sc.node.r <- prop.sc.train[train.data[, col.ind] < split.used]
              } 
            }
            
          } else {
            
            data.node.l <- train.data[train.data[,  col.ind] < split.used, ]
            data.node.r <- train.data[train.data[,  col.ind] >= split.used, ]
            
            if (propsc.mod.out){
              prop.sc.node.l <- prop.sc.train[train.data[, col.ind] < split.used]
              prop.sc.node.r <- prop.sc.train[train.data[, col.ind] >= split.used]
            } else {
              if (nrow(data.node.l) < min.obs.mod) {
                prop.sc.node.l <- prop.sc.train[train.data[, col.ind] < split.used]
              }
              
              if (nrow(data.node.r) < min.obs.mod) {
                prop.sc.node.r <- prop.sc.train[train.data[, col.ind] >= split.used]
              }
            }
            
          }
        }
        
        if (!propsc.mod.out){
          if (nrow(data.node.l) >= min.obs.mod) {
            
            data.node.l.noy <- data.node.l[, colnames(data.node.l) != "Y"]
            tmp <- gen.fullrank.ipw(df.noy           = data.node.l.noy, 
                                    propsc.form.true = propsc.form.true)
            
            tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                               method    = propsc.mthd,
                               form.true = tmp$propsc.form.true.updated)
            prop.sc.node.l <- tmp$prop.sc
            
          }
          
          if (nrow(data.node.r) >= min.obs.mod) {
            
            data.node.r.noy <- data.node.r[, colnames(data.node.r) != "Y"]
            tmp <- gen.fullrank.ipw(df.noy           = data.node.r.noy, 
                                    propsc.form.true = propsc.form.true)
            
            tmp <- est.prop.sc(df.noy    = tmp$df.noy.fullrank,
                               method    = propsc.mthd,
                               form.true = tmp$propsc.form.true.updated)
            prop.sc.node.r <- tmp$prop.sc
  
          }
        }
        
        ### Left node
        # Calculating unadjusted estimator  
        mu.1l <- mean(data.node.l$Y * data.node.l$A/ prop.sc.node.l)
        mu.0l <- mean(data.node.l$Y * (1 - data.node.l$A)/ (1 - prop.sc.node.l))
        
        ### Right node
        # Calculating unadjusted estimator  
        mu.1r <- mean(data.node.r$Y * data.node.r$A/ prop.sc.node.r)
        mu.0r <- mean(data.node.r$Y * (1 - data.node.r$A)/ (1 - prop.sc.node.r))
        
        col.ind <- which(colnames(test.data) == var.used)
        if ((split.used %% 1) == 0){   
          
          # Test set prediction on categorical covariate split
          lvls <- levels(train.data[, col.ind])
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]] <- mu.1l - mu.0l
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]] <- mu.1r - mu.0r
          
        } else{
          # Test set prediction on continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] < split.used] <- mu.1r - mu.0r
          } else {
            pred.tree[test.data[, col.ind] < split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1r - mu.0r
          }
          
        }
        
      } else {
        
        mu.1 <- mean(train.data$Y * train.data$A/ prop.sc.train)
        mu.0 <- mean(train.data$Y * (1 - train.data$A)/ (1 - prop.sc.train))   
        pred.tree = mu.1 - mu.0
      }
      
      cv.err[m, l] <- mean((pred.tree - treat.eff.test)^2)
      
    }

  }
  
  # Averaging over cross validation sets
  cv.err.fin = apply(cv.err, 1, mean)
  tree.final = tree.list[[which(cv.err.fin == min(cv.err.fin))[length(which(cv.err.fin == min(cv.err.fin)))]]]
  
  return(list(tree.final, cv.err.fin))
  
}



######################################################################################################################
############################################ CV method 2, G estimation ###############################################
######################################################################################################################
EstG.CvMethod2 = function(data.used, tree.list, type.var = "cont", seed = NULL, n.cv = 5,
                          adj.mod.out = T, adj.mthd = "GLM", adj.form.true = NULL, min.obs.mod = NULL){
  
  if (!is.null(seed)){
    set.seed(seed)
  }
  
  cross.val.ind = sample(cut(seq(1,nrow(data.used)), breaks=n.cv, labels=FALSE), 
                         nrow(data.used))
  cv.err = matrix(NA, ncol = n.cv, nrow = length(tree.list))
  
  for (l in 1:n.cv){
    
    test.data = data.used[which(cross.val.ind == l), ]
    train.data = data.used[which(cross.val.ind != l), ]
    
    # Fit the random forest procedure on the test data to find the gold standard
    train.data.noy <- train.data %>%
      dplyr::select(-Y)
    train.data.noy <- train.data.noy %>%
      mutate(prop.sc = est.prop.sc(train.data.noy, method = "GLM")$prop.sc) %>%
      mutate(wt      = ifelse(A == 1, prop.sc, 1 - prop.sc))
    
    rand.for.train = rfsrc(Y~., data = train.data, case.wt = 1 / train.data.noy$wt)
    test.A.1 = test.data
    test.A.1$A = 1
    pred.A.1 = predict(rand.for.train, newdata = test.A.1)
    test.A.0 = test.data
    test.A.0$A = 0
    pred.A.0 = predict(rand.for.train, newdata = test.A.0)
    treat.eff.test = pred.A.1$predicted - pred.A.0$predicted
    
    # train
    tmp <- est.cond.eff(df        = train.data,
                        method    = adj.mthd,
                        form.true = adj.form.true,
                        type.var  = type.var)
    
    est.cond.eff.1.train <- tmp$pred.A.1
    est.cond.eff.0.train <- tmp$pred.A.0
    
    # Calculating g(h) for the cross.validated tree
    # Looping through the candidate trees
    for(m in 1:length(tree.list)){
      
      tree.used = tree.list[[m]]
      
      if (nrow(tree.used$frame) >3) {
        
        pred.train <- predict(tree.used, newdata = train.data)
        pred.tree <- predict(tree.used, newdata = test.data)
        
        for(k in 1:length(unique(pred.train))){
          # Calculating the treatment effect estimator for each node using only the 
          # training data
          data.node = train.data[which(pred.train == unique(pred.train)[k]), ]
          test.index = which(pred.tree == unique(pred.train)[k])
          
          if (!adj.mod.out){
            
            if (nrow(data.node) >= min.obs.mod){  # When prop.mod.out = T, min.obs.mod = NULL
              
              tmp <- gen.fullrank.g(df            = data.node,
                                    adj.form.true = adj.form.true)
              tmp <- est.cond.eff(df        = tmp$df.fullrank,
                                  method    = adj.mthd,
                                  form.true = tmp$adj.form.true.updated,
                                  type.var  = type.var)
              
              est.cond.eff.1.node <- tmp$pred.A.1
              est.cond.eff.0.node <- tmp$pred.A.0
            } else {
              est.cond.eff.1.node <- est.cond.eff.1.train[which(pred.train == unique(pred.train)[k])]
              est.cond.eff.0.node <- est.cond.eff.0.train[which(pred.train == unique(pred.train)[k])]
            }
            
          } else { # model prop.sc.out / the number of observations in the node small
            est.cond.eff.1.node <- est.cond.eff.1.train[which(pred.train == unique(pred.train)[k])]
            est.cond.eff.0.node <- est.cond.eff.0.train[which(pred.train == unique(pred.train)[k])]
          }
          
          mu.1 <- mean(est.cond.eff.1.node)
          mu.0 <- mean(est.cond.eff.0.node)
          
          pred.tree[test.index] <- mu.1 - mu.0
          
        }
        
      } else if (nrow(tree.used$frame) == 3){
        # If there is one or zero splits there is a weird memory error so need to do manually
        pred.tree = rep(NA, nrow(test.data))
        split.used = tree.used$splits[, 4]
        var.used = tree.used$frame$var[1]
        col.ind <- which(colnames(train.data) == var.used)
        
        # Need to figure out observations going to the left/right node
        if ((split.used %% 1) == 0){   
          
          # Categorical covariate split
          lvls <- levels(train.data[, col.ind])
          data.node.l <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1], ]
          data.node.r <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          
          if (adj.mod.out){ # When prop.mod.out = T, min.obs.mod = NULL
            
            est.cond.eff.1.node.l <- est.cond.eff.1.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]]
            est.cond.eff.0.node.l <- est.cond.eff.0.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]]
            est.cond.eff.1.node.r <- est.cond.eff.1.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]]
            est.cond.eff.0.node.r <- est.cond.eff.0.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]]
            
          } else { 
            if (nrow(data.node.l) < min.obs.mod) {
              est.cond.eff.1.node.l <- est.cond.eff.1.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]]
              est.cond.eff.0.node.l <- est.cond.eff.0.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]]
            }
            
            if  (nrow(data.node.r) < min.obs.mod) {
              est.cond.eff.1.node.r <- est.cond.eff.1.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]]
              est.cond.eff.0.node.r <- est.cond.eff.0.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]]
            }
          }
          
        } else {
          # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            data.node.l <- train.data[train.data[,  col.ind] >= split.used, ]
            data.node.r <- train.data[train.data[,  col.ind] < split.used, ]
            
            if (adj.mod.out){
              est.cond.eff.1.node.l <- est.cond.eff.1.train[train.data[,  col.ind] >= split.used]
              est.cond.eff.0.node.l <- est.cond.eff.0.train[train.data[,  col.ind] >= split.used]
              est.cond.eff.1.node.r <- est.cond.eff.1.train[train.data[,  col.ind] < split.used]
              est.cond.eff.0.node.r <- est.cond.eff.0.train[train.data[,  col.ind] < split.used]

            } else {
              if (nrow(data.node.l) < min.obs.mod) {
                est.cond.eff.1.node.l <- est.cond.eff.1.train[train.data[,  col.ind] >= split.used]
                est.cond.eff.0.node.l <- est.cond.eff.0.train[train.data[,  col.ind] >= split.used]
              }
              
              if (nrow(data.node.r) < min.obs.mod){
                est.cond.eff.1.node.r <- est.cond.eff.1.train[train.data[,  col.ind] < split.used]
                est.cond.eff.0.node.r <- est.cond.eff.0.train[train.data[,  col.ind] < split.used]
              } 
            }
            
          } else {
            
            data.node.l <- train.data[train.data[,  col.ind] < split.used, ]
            data.node.r <- train.data[train.data[,  col.ind] >= split.used, ]
            
            if (adj.mod.out){
              est.cond.eff.1.node.l <- est.cond.eff.1.train[train.data[,  col.ind] < split.used]
              est.cond.eff.0.node.l <- est.cond.eff.0.train[train.data[,  col.ind] < split.used]
              est.cond.eff.1.node.r <- est.cond.eff.1.train[train.data[,  col.ind] >= split.used]
              est.cond.eff.0.node.r <- est.cond.eff.0.train[train.data[,  col.ind] >= split.used]

            } else {
              if (nrow(data.node.l) < min.obs.mod) {
                est.cond.eff.1.node.l <- est.cond.eff.1.train[train.data[,  col.ind] < split.used]
                est.cond.eff.0.node.l <- est.cond.eff.0.train[train.data[,  col.ind] < split.used]
              }
              
              if (nrow(data.node.r) < min.obs.mod) {
                est.cond.eff.1.node.r <- est.cond.eff.1.train[train.data[,  col.ind] >= split.used]
                est.cond.eff.0.node.r <- est.cond.eff.0.train[train.data[,  col.ind] >= split.used]
              }
            }
            
          }
        }
        
        if (!adj.mod.out){
          if (nrow(data.node.l) >= min.obs.mod) {
            
            tmp <- gen.fullrank.g(df            = data.node.l,
                                  adj.form.true = adj.form.true)
            tmp <- est.cond.eff(df        = tmp$df.fullrank,
                                method    = adj.mthd,
                                form.true = tmp$adj.form.true.updated,
                                type.var  = type.var)
            
            est.cond.eff.1.node.l <- tmp$pred.A.1
            est.cond.eff.0.node.l <- tmp$pred.A.0
          }
          
          if (nrow(data.node.r) >= min.obs.mod) {
            
            tmp <- gen.fullrank.g(df            = data.node.r,
                                  adj.form.true = adj.form.true)
            tmp <- est.cond.eff(df        = tmp$df.fullrank,
                                method    = adj.mthd,
                                form.true = tmp$adj.form.true.updated,
                                type.var  = type.var)
            
            est.cond.eff.1.node.r <- tmp$pred.A.1
            est.cond.eff.0.node.r <- tmp$pred.A.0
          }
        }
        
        ### Left node
        # Calculating unadjusted estimator  
        mu.1l <- mean(est.cond.eff.1.node.l)
        mu.0l <- mean(est.cond.eff.0.node.l)
        
        ### Right node
        # Calculating unadjusted estimator  
        mu.1r <- mean(est.cond.eff.1.node.r)
        mu.0r <- mean(est.cond.eff.0.node.r)
        
        col.ind <- which(colnames(test.data) == var.used)
        if ((split.used %% 1) == 0){   
          
          # Test set prediction on categorical covariate split
          lvls <- levels(train.data[, col.ind])
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]] <- mu.1l - mu.0l
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]] <- mu.1r - mu.0r
          
        } else{
          # Test set prediction on continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] < split.used] <- mu.1r - mu.0r
          } else {
            pred.tree[test.data[, col.ind] < split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1r - mu.0r
          }
          
        }
        
      } else { # nrow(tree.used$frame) == 1
        
        mu.1 <- mean(est.cond.eff.1.train)
        mu.0 <- mean(est.cond.eff.0.train)   
        pred.tree <- mu.1 - mu.0
      }
      
      cv.err[m, l] <- mean((pred.tree - treat.eff.test)^2)
    }

  }

  # Averaging over cross validation sets
  cv.err.fin = apply(cv.err, 1, mean)
  tree.final = tree.list[[which(cv.err.fin == min(cv.err.fin))[length(which(cv.err.fin == min(cv.err.fin)))]]]
  
  return(list(tree.final, cv.err.fin))
    
}


######################################################################################################################
########################################### CV method 2, DR estimation ###############################################
######################################################################################################################
EstDr.CvMethod2 = function(data.used, tree.list, type.var = "cont", seed = NULL, n.cv = 5,
                           propsc.mod.out = T, propsc.mthd = "GLM", propsc.form.true = NULL, min.obs.propsc = NULL,
                           adj.mod.out = T, adj.mthd = "GLM", adj.form.true = NULL, min.obs.adj = NULL){
  
  if (!is.null(seed)){
    set.seed(seed)
  }
  
  cross.val.ind = sample(cut(seq(1,nrow(data.used)), breaks=n.cv, labels=FALSE), 
                         nrow(data.used))
  cv.err = matrix(NA, ncol = n.cv, nrow = length(tree.list))
  
  for (l in 1:n.cv){
    
    test.data = data.used[which(cross.val.ind == l), ]
    train.data = data.used[which(cross.val.ind != l), ]
    
    # Fit the random forest procedure on the test data to find the gold standard
    train.data.noy <- train.data %>%
      dplyr::select(-Y)
    train.data.noy <- train.data.noy %>%
      mutate(prop.sc = est.prop.sc(train.data.noy, method = "GLM")$prop.sc) %>%
      mutate(wt      = ifelse(A == 1, prop.sc, 1 - prop.sc))
    
    rand.for.train = rfsrc(Y~., data = train.data, case.wt = 1 / train.data.noy$wt)
    test.A.1 = test.data
    test.A.1$A = 1
    pred.A.1 = predict.rfsrc(rand.for.train, newdata = test.A.1)
    test.A.0 = test.data
    test.A.0$A = 0
    pred.A.0 = predict(rand.for.train, newdata = test.A.0)
    treat.eff.test = pred.A.1$predicted - pred.A.0$predicted
    
    # Calculate the conditional mean if adj.mod fitted outside
    tmp <- gen.fullrank.g(df            = train.data,
                          adj.form.true = adj.form.true)
    tmp <- est.cond.eff(df        = tmp$df.fullrank,
                        method    = adj.mthd,
                        form.true = tmp$adj.form.true.updated, 
                        type.var  = type.var)
    
    est.cond.eff.1.train <- tmp$pred.A.1
    est.cond.eff.0.train <- tmp$pred.A.0
    
    # Calculate prop.sc if prop.sc model fitted outside
    train.data.noy <- train.data %>%
      dplyr::select(-Y)
    tmp <- gen.fullrank.ipw(df.noy           = train.data.noy, 
                            propsc.form.true = propsc.form.true)
    prop.sc.train <- est.prop.sc(df.noy    = tmp$df.noy.fullrank, 
                                 method    = propsc.mthd,
                                 form.true = tmp$propsc.form.true.updated)$prop.sc
    
    # Calculating g(h) for the cross.validated tree
    # Looping through the candidate trees
    for(m in 1:length(tree.list)){
      
      tree.used = tree.list[[m]]
      
      if (nrow(tree.used$frame) >3) {
        pred.train <- predict(tree.used, newdata = train.data)
        pred.tree <- predict(tree.used, newdata = test.data)
        
        for(k in 1:length(unique(pred.train))){
          # Calculating the treatment effect estimator for each node using only the 
          # training data
          data.node = train.data[which(pred.train == unique(pred.train)[k]), ]
          test.index = which(pred.tree == unique(pred.train)[k])
          
          # Fit adj.mod if not fit outside
          if (!adj.mod.out){
            
            if (nrow(data.node) >= min.obs.adj){  # When prop.mod.out = T, min.obs.mod = NULL
              
              tmp <- gen.fullrank.g(df            = data.node,
                                    adj.form.true = adj.form.true)
              tmp <- est.cond.eff(df        = tmp$df.fullrank,
                                  method    = adj.mthd,
                                  form.true = tmp$adj.form.true.updated,
                                  type.var  = type.var)
              
              est.cond.eff.1.node <- tmp$pred.A.1
              est.cond.eff.0.node <- tmp$pred.A.0
    
            } else { # If the number of observations fewer than min.obs.adj, use the result fitted outside to replace
              est.cond.eff.1.node <- est.cond.eff.1.train[which(pred.train == unique(pred.train)[k])]
              est.cond.eff.0.node <- est.cond.eff.0.train[which(pred.train == unique(pred.train)[k])]
            }
            
          } else { # model prop.sc.out / the number of observations in the node small
            est.cond.eff.1.node <- est.cond.eff.1.train[which(pred.train == unique(pred.train)[k])]
            est.cond.eff.0.node <- est.cond.eff.0.train[which(pred.train == unique(pred.train)[k])]
          }
          
          # Fit propsc.mod if not fit outside
          if (!propsc.mod.out){
            
            if (nrow(data.node) >= min.obs.propsc){  # When prop.mod.out = T, min.obs.mod = NULL
              data.node.noy <- data.node %>%
                dplyr::select(-Y)
              tmp <- gen.fullrank.ipw(df.noy           = data.node.noy, 
                                      propsc.form.true = propsc.form.true)
              prop.sc.node <- est.prop.sc(df.noy    = tmp$df.noy.fullrank, 
                                          method    = propsc.mthd,
                                          form.true = tmp$propsc.form.true.updated)$prop.sc
            } else {
              prop.sc.node <- prop.sc.train[which(pred.train == unique(pred.train)[k])]
            }
            
          } else { # model prop.sc.out / the number of observations in the node small
            prop.sc.node <- prop.sc.train[which(pred.train == unique(pred.train)[k])]
          }
  
          mu.1 <- mean(data.node$A * (data.node$Y - est.cond.eff.1.node) / prop.sc.node + est.cond.eff.1.node)
          mu.0 <- mean((1 - data.node$A) * (data.node$Y - est.cond.eff.0.node)/ (1 - prop.sc.node) + est.cond.eff.0.node)

          pred.tree[test.index] <- mu.1 - mu.0
          
        }
        
      } else if (nrow(tree.used$frame) == 3) {
        # If there is one or zero splits there is a weird memory error so need to do manually
        pred.tree = rep(NA, nrow(test.data))
        split.used = tree.used$splits[, 4]
        var.used = tree.used$frame$var[1]
        col.ind <- which(colnames(train.data) == var.used)
        
        # Need to figure out observations going to the left/right node
        if (class(data.used[, as.character(var.used)]) == "factor"){   
          
          # Categorical covariate split
          lvls <- levels(train.data[, col.ind])
          data.node.l <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1], ]
          data.node.r <- train.data[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3], ]
          
          if (adj.mod.out){ # When prop.mod.out = T, min.obs.mod = NULL
            
            est.cond.eff.1.node.l <- est.cond.eff.1.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]]
            est.cond.eff.0.node.l <- est.cond.eff.0.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]]
            est.cond.eff.1.node.r <- est.cond.eff.1.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]]
            est.cond.eff.0.node.r <- est.cond.eff.0.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]]
            
          } else { 
            if (nrow(data.node.l) < min.obs.adj) {
              est.cond.eff.1.node.l <- est.cond.eff.1.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]]
              est.cond.eff.0.node.l <- est.cond.eff.0.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]]
            }
            
            if  (nrow(data.node.r) < min.obs.adj) {
              est.cond.eff.1.node.r <- est.cond.eff.1.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]]
              est.cond.eff.0.node.r <- est.cond.eff.0.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]]
            }
          }
          
          if (propsc.mod.out){ # When prop.mod.out = T, min.obs.mod = NULL
            prop.sc.node.l <- prop.sc.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]]
            prop.sc.node.r <- prop.sc.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]]
          } else { 
            if (nrow(data.node.l) < min.obs.propsc) {
              prop.sc.node.l <- prop.sc.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]]
            }
            
            if  (nrow(data.node.r) < min.obs.propsc) {
              prop.sc.node.r <- prop.sc.train[train.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]]
            }
          }
          
        } else { # Continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            data.node.l <- train.data[train.data[,  col.ind] >= split.used, ]
            data.node.r <- train.data[train.data[,  col.ind] < split.used, ]
            
            if (adj.mod.out){
              est.cond.eff.1.node.l <- est.cond.eff.1.train[train.data[,  col.ind] >= split.used]
              est.cond.eff.0.node.l <- est.cond.eff.0.train[train.data[,  col.ind] >= split.used]
              est.cond.eff.1.node.r <- est.cond.eff.1.train[train.data[,  col.ind] < split.used]
              est.cond.eff.0.node.r <- est.cond.eff.0.train[train.data[,  col.ind] < split.used]
              
            } else {
              if (nrow(data.node.l) < min.obs.adj) {
                est.cond.eff.1.node.l <- est.cond.eff.1.train[train.data[,  col.ind] >= split.used]
                est.cond.eff.0.node.l <- est.cond.eff.0.train[train.data[,  col.ind] >= split.used]
              }
              
              if (nrow(data.node.r) < min.obs.adj){
                est.cond.eff.1.node.r <- est.cond.eff.1.train[train.data[,  col.ind] < split.used]
                est.cond.eff.0.node.r <- est.cond.eff.0.train[train.data[,  col.ind] < split.used]
              } 
            }
            
            if (propsc.mod.out){
              prop.sc.node.l <- prop.sc.train[train.data[, col.ind] >= split.used]
              prop.sc.node.r <- prop.sc.train[train.data[, col.ind] < split.used]
            } else {
              if (nrow(data.node.l) < min.obs.propsc) {
                prop.sc.node.l <- prop.sc.train[train.data[, col.ind] >= split.used]
              }
              
              if (nrow(data.node.r) < min.obs.propsc){
                prop.sc.node.r <- prop.sc.train[train.data[, col.ind] < split.used]
              } 
            }
            
          } else {
            
            data.node.l <- train.data[train.data[,  col.ind] < split.used, ]
            data.node.r <- train.data[train.data[,  col.ind] >= split.used, ]
            
            if (adj.mod.out){
              est.cond.eff.1.node.l <- est.cond.eff.1.train[train.data[,  col.ind] < split.used]
              est.cond.eff.0.node.l <- est.cond.eff.0.train[train.data[,  col.ind] < split.used]
              est.cond.eff.1.node.r <- est.cond.eff.1.train[train.data[,  col.ind] >= split.used]
              est.cond.eff.0.node.r <- est.cond.eff.0.train[train.data[,  col.ind] >= split.used]
              
            } else {
              if (nrow(data.node.l) < min.obs.adj) {
                est.cond.eff.1.node.l <- est.cond.eff.1.train[train.data[,  col.ind] < split.used]
                est.cond.eff.0.node.l <- est.cond.eff.0.train[train.data[,  col.ind] < split.used]
              }
              
              if (nrow(data.node.r) < min.obs.adj) {
                est.cond.eff.1.node.r <- est.cond.eff.1.train[train.data[,  col.ind] >= split.used]
                est.cond.eff.0.node.r <- est.cond.eff.0.train[train.data[,  col.ind] >= split.used]
              }
            }
            
            if (propsc.mod.out){
              prop.sc.node.l <- prop.sc.train[train.data[, col.ind] < split.used]
              prop.sc.node.r <- prop.sc.train[train.data[, col.ind] >= split.used]
            } else {
              if (nrow(data.node.l) < min.obs.propsc) {
                prop.sc.node.l <- prop.sc.train[train.data[, col.ind] < split.used]
              }
              
              if (nrow(data.node.r) < min.obs.propsc) {
                prop.sc.node.r <- prop.sc.train[train.data[, col.ind] >= split.used]
              }
            }
            
          }  # else [if (tree.used$splits[2] > 0)]
        } # else [if ((split.used %% 1) == 0)]
        
        if (!adj.mod.out){
          if (nrow(data.node.l) >= min.obs.adj) {
            tmp <- gen.fullrank.g(df            = data.node.l,
                                  adj.form.true = adj.form.true)
            tmp <- est.cond.eff(df        = tmp$df.fullrank,
                                method    = adj.mthd,
                                form.true = tmp$adj.form.true.updated,
                                type.var  = type.var)
            
            est.cond.eff.1.node.l <- tmp$pred.A.1
            est.cond.eff.0.node.l <- tmp$pred.A.0
          }
          
          if (nrow(data.node.r) >= min.obs.adj) {
            tmp <- gen.fullrank.g(df            = data.node.r,
                                  adj.form.true = adj.form.true)
            tmp <- est.cond.eff(df        = tmp$df.fullrank,
                                method    = adj.mthd,
                                form.true = tmp$adj.form.true.updated,
                                type.var  = type.var)
            
            est.cond.eff.1.node.r <- tmp$pred.A.1
            est.cond.eff.0.node.r <- tmp$pred.A.0
          }
        }
        
        if (!propsc.mod.out){
          if (nrow(data.node.l) >= min.obs.propsc) {
            data.node.l.noy <- data.node.l %>%
              dplyr::select(-Y)
            tmp <- gen.fullrank.ipw(df.noy           = data.node.l.noy, 
                                    propsc.form.true = propsc.form.true)
            prop.sc.node.l <- est.prop.sc(df.noy    = tmp$df.noy.fullrank, 
                                          method    = propsc.mthd,
                                          form.true = tmp$propsc.form.true.updated)$prop.sc
          }
          
          if (nrow(data.node.r) >= min.obs.propsc) {
            data.node.r.noy <- data.node.r %>%
              dplyr::select(-Y)
            tmp <- gen.fullrank.ipw(df.noy           = data.node.r.noy, 
                                    propsc.form.true = propsc.form.true)
            prop.sc.node.r <- est.prop.sc(df.noy    = tmp$df.noy.fullrank, 
                                          method    = propsc.mthd,
                                          form.true = tmp$propsc.form.true.updated)$prop.sc
          }
        }
        
        ### Left node
        # Calculating unadjusted estimator  
        mu.1l <- mean(data.node.l$A * (data.node.l$Y - est.cond.eff.1.node.l) / prop.sc.node.l + est.cond.eff.1.node.l)
        mu.0l <- mean((1 - data.node.l$A) * (data.node.l$Y - est.cond.eff.0.node.l)/ (1 - prop.sc.node.l) + est.cond.eff.0.node.l)

        ### Right node
        # Calculating unadjusted estimator  
        mu.1r <- mean(data.node.r$A * (data.node.r$Y - est.cond.eff.1.node.r) / prop.sc.node.r + est.cond.eff.1.node.r)
        mu.0r <- mean((1 - data.node.r$A) * (data.node.r$Y - est.cond.eff.0.node.r)/ (1 - prop.sc.node.r) + est.cond.eff.0.node.r)
        
        col.ind <- which(colnames(test.data) == var.used)
        if (class(data.used[, as.character(var.used)]) == "factor"){   
          
          # Test set prediction on categorical covariate split
          lvls <- levels(train.data[, col.ind])
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 1]] <- mu.1l - mu.0l
          pred.tree[test.data[, col.ind] %in% lvls[tree.used$csplit[split.used,] == 3]] <- mu.1r - mu.0r
          
        } else {
          # Test set prediction on continuous covariate split
          # Need to take care of left or right
          if (tree.used$splits[2] > 0) {
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] < split.used] <- mu.1r - mu.0r
          } else {
            pred.tree[test.data[, col.ind] < split.used] <- mu.1l - mu.0l
            pred.tree[test.data[, col.ind] >= split.used] <- mu.1r - mu.0r
          }
          
        }
        
      } else { # else if (nrow(tree.used$frame) == 3) {

        mu.1 <- mean(train.data$A * (train.data$Y - est.cond.eff.1.train) / prop.sc.train + est.cond.eff.1.train)
        mu.0 <- mean((1 - train.data$A) * (train.data$Y - est.cond.eff.0.train)/ (1 - prop.sc.train) + est.cond.eff.0.train)

        pred.tree <- mu.1 - mu.0
      }
      
      # If there is only 1 observation in the terminal node of the tree, NA will be produced in calculation
      cv.err[m, l] <- mean((pred.tree - treat.eff.test)^2, na.rm = T)
    }
    
  }
  
  # Averaging over cross validation sets
  cv.err.fin = apply(cv.err, 1, mean)
  tree.final = tree.list[[which(cv.err.fin == min(cv.err.fin))[length(which(cv.err.fin == min(cv.err.fin)))]]]
  
  return(list(tree.final, cv.err.fin))
  
}
