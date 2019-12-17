eval.measures.eff = function(final.tree, test.data, true.trt.eff, noise.var, 
                             corr.split, where.split, dir.split, split.cate = NULL, CT = F){
  
  # split.cate used for categorical splits. Specify which categories should be in the same node for categorical variables
  
  # Calculate the size of tree 
  size.tree <- sum(final.tree$frame$var == "<leaf>")
  
  # Calculate the number of splits made
  num.splt <- sum(final.tree$frame$var != "<leaf>")
  
  # Calcualte the number/proportion of times the tree splits on each variables 
  # Number of Noise variables
  numb.noise <- sum(final.tree$frame$var %in% noise.var)
  if (num.splt == 0){
    prop.noise.all <- 0                      # NEED TO TAKE CARE OF THIS WHEN NO SPLITS SHOULD BE MADE
  } else {
    prop.noise.all <- numb.noise / num.splt
  }
  
  # Number of correct trees
  if (is.null(corr.split)){ 
    exact.corr <- sum(size.tree == 1)
    
  } else if (size.tree != (length(corr.split) + 1)) { # First confirm the size of the tree is correct
    exact.corr <- 0
    
  } else {
    
    for (i in 1:length(where.split)) {
      
      # CURRENTLY ONLY MAKING 1 SPLIT WORKS
      cond <- T
      
      for (j in 1:length(where.split[[i]])) {
        
        cond <- cond & (final.tree$frame$var[where.split[[i]][j]] == corr.split[j])
        
        # This condition is used for when there are more than 1 split
        # make sure the direction of the first split is correct so that the second split can be made correctly
        if (cond){
          if (!is.null(dir.split[[i]][j])){
            cond <- cond & (final.tree$splits[j, 2] == dir.split[[i]][j])
          }
          
          if (!cond){
            break
          }
          
        } else {
          break
        }
        
        # This condition is used for categorical splits
        if (!is.null(final.tree$csplit)) {
          
          # which is the current split
          if (CT) { # in CausalTree, the split.stat doesn't exist; need to use direct way to specify the indices
            split.used <- final.tree$splits[1, 4]
          } else {
            row.ind <- which(final.tree$splits[, 3] == final.tree$frame$split.stat[where.split[[i]][j]])
            split.used <- final.tree$splits[row.ind, 4]
          }
          
          if (split.used %% 1 == 0) {
            
            tmp.csplit <- final.tree$csplit[split.used,]
            cond <- cond & (length(unique(tmp.csplit[split.cate[[i]]])) == 1)
            
            rest.csplit <- tmp.csplit[-split.cate[[i]]]
            rest.csplit <- rest.csplit[rest.csplit != 2]
            cond <- cond & (length(unique(rest.csplit)) == 1)
            
            if (!cond) {
              break
            }
          } # do nothing if the current split is not categorical split
          
        } # end categorical split evaluation
        
      } # for loop j
      
      if (cond){
        break
      }
      
    } # for loop i
      
    # make sure there is no noise variable if splits are made correctly
    if (cond){
      for (i in 1:length(noise.var)){
        cond <- cond & (sum(final.tree$frame$var == noise.var[i]) == 0)
      }
    }
    
    exact.corr <- sum(cond)
    
  }
  
  # Number/proportion of correct splits made on the first few splits
  num.corr.splt.frst <- sum(final.tree$frame[as.character(1:(2^(length(corr.split)) - 1)), ]$var %in% corr.split)
  prop.corr.frst     <- num.corr.splt.frst / (2^(length(corr.split)) - 1)
  # Proportion of correct splits over the number of splits made
  if (num.splt == 0){
    prop.corr.frst.tree <- 0                      # NEED TO TAKE CARE OF THIS WHEN NO SPLITS SHOULD BE MADE
  } else {
    prop.corr.frst.tree <- num.corr.splt.frst / num.splt
  }
  
  # Number/proportion of splits made on the signal variables
  num.corr.splt.all <- sum(final.tree$frame$var %in% corr.split)
  if (num.splt == 0){
    prop.corr.all <- 0                      # NEED TO TAKE CARE OF THIS WHEN NO SPLITS SHOULD BE MADE
  } else {
    prop.corr.all <- num.corr.splt.all / num.splt
  }
  
  # Calcualte the mean square error
  # Calculate the prediction on the test data
  if (nrow(final.tree$frame) > 3) {
    pred.tree = predict(final.tree, newdata = test.data)
  } else if (nrow(final.tree$frame) == 3){
    # If there is one or zero splits there is a weird memory error so need to do manually
    pred.tree = rep(NA, nrow(test.data))
    split.used = final.tree$splits[, 4]
    var.used = final.tree$frame$var[1]
    col.ind <- which(colnames(test.data) == var.used)
    
    if (length(split.used) > 1) {
      split.used <- split.used[1]
      direction <- final.tree$splits[1, 2]
    } else {
      direction <- final.tree$splits[2]
    }
    
    # Need to figure out observations going to the left/right node
    if ((split.used %% 1) == 0){   
      
      # Categorical covariate split
      lvls <- levels(test.data[, col.ind])
      pred.tree[test.data[, col.ind] %in% lvls[final.tree$csplit[split.used,] == 1]] <- final.tree$frame$yval[2]
      pred.tree[test.data[, col.ind] %in% lvls[final.tree$csplit[split.used,] == 3]] <- final.tree$frame$yval[3]
      
    } else{
      # Continuous covariate split
      # Need to take care of left or right
      if (direction > 0) {
        pred.tree[test.data[,  col.ind] >= split.used] <- final.tree$frame$yval[2]
        pred.tree[test.data[,  col.ind] < split.used]  <- final.tree$frame$yval[3]
      } else {
        pred.tree[test.data[,  col.ind] < split.used]  <- final.tree$frame$yval[2]
        pred.tree[test.data[,  col.ind] >= split.used] <- final.tree$frame$yval[3]
      }
      
    }
  } else {
    pred.tree = rep(final.tree$frame$yval, nrow(test.data))
  } 
  # mse will be NA if there is only 1 observation in the terminal node
  mse = mean((pred.tree - true.trt.eff)^2, na.rm = T)
  
  # PPS
  pps = 1
  for(i in 1:(nrow(test.data)-1)){
    for(j in (i+1):nrow(test.data)){
      a=b=0
      if(true.trt.eff[i] == true.trt.eff[j]){a = 1}
      
      # Treat observations in the 1-observation terminal node be in the same node
      if(is.na(pred.tree[i]) | is.na(pred.tree[j])) {
        b = 1
      } else if (pred.tree[i] == pred.tree[j]) {
        b = 1
      }
      pps = pps - abs(a-b)/choose(nrow(test.data), 2)
    }
  }
  
  return(list(mse                 = mse, 
              exact.corr          = exact.corr, 
              size.tree           = size.tree,
              num.splt            = num.splt,
              numb.noise          = numb.noise,
              pps                 = pps,
              prop.noise.all      = prop.noise.all,
              num.corr.splt.frst  = num.corr.splt.frst,
              prop.corr.frst      = prop.corr.frst, 
              prop.corr.frst.tree = prop.corr.frst.tree,
              num.corr.splt.all   = num.corr.splt.all,
              prop.corr.all       = prop.corr.all))
}

# Evaluate whether the first split is made correctly if categorical split
eval.cate.corr.frst.splt <- function(large.tree, corr.split, split.cate) {
  
  cond <- as.character(large.tree$frame$var[1]) == corr.split
  
  # Since we only look at the first split, the index is deterministic and equals 1
  tmp.csplit <- large.tree$csplit[1, ]
  cond <- cond & (length(unique(tmp.csplit[split.cate[[1]]])) == 1)
  
  rest.csplit <- tmp.csplit[-split.cate[[1]]]
  rest.csplit <- rest.csplit[rest.csplit != 2]
  cond <- cond & (length(unique(rest.csplit)) == 1)
  
  return(cond)
  
}
