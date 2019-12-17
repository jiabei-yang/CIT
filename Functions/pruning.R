create.sequence = function(data.used, est.used, type.var = "cont",        # For ipw, type.var does not matter
                           propsc.mod.out = T, propsc.mthd = "GLM", propsc.form.true = NULL, propsc.mod.insplt = NULL, # params for ipw, dr
                           adj.mod.out = T, adj.mthd = "GLM", adj.form.true = NULL, adj.mod.insplt = NULL,             # params for g, dr
                           w = NULL,                                                 
                           num.truc.obs = 30, min.node = 30){
  
  # This function calculates the sequence of candidate trees.
  # Input: data.used the dataset with treatment as first column, and outcome as second column
  # est.used: c("IPW", "G", "DR")
  # need.cond.exp = a logical statement if an outside estimator for the conditional expectation
  # is needed
  # type.var = "cont" for continuous outcomes, "bin" for binary outcomes
  # Output: A list with two elements tree.list which is the sequence of candidate trees 
  # and lambda.list which is the sequence of penalization parameters
  
  # Need different formulas for different outcomes
  # if(type.var == "cont"){
  #   form.used <- as.formula(paste("Y ~ ", paste(names(data.used)[!(names(data.used) %in% c("Y", "A"))], collapse= "+")))
  # } else if(type.var == "bin"){
  
  data.used.covariates  <- data.used[!colnames(data.used) %in% c("Y", "A")]
  form.used             <- as.formula(paste("rownumb ~ ", paste(colnames(data.used.covariates), collapse= "+")))
  data.used.wiz.rownumb <- cbind(rownumb = 1:dim(data.used)[1], data.used)
  data.used.noy         <- data.used[, !colnames(data.used) %in% c("Y")]
  
  # }
  
  # Creating ulist.used and parameter vector used
  if (est.used == "IPW"){
    ulist.used <- list(eval = etemp.ipw, split = stemp.ipw, init = itemp)
    # ulist.used <- list(eval = etemp.ipw.whl, split = stemp.ipw.whl, init = itemp)
        
    if (propsc.mod.out){
      
      tmp <- est.prop.sc(df.noy    = data.used.noy,
                         method    = propsc.mthd,
                         form.true = propsc.form.true)
      parms.used <- list(trt          = data.used$A,
                         covariates   = data.used.covariates,
                         response     = data.used$Y,
                         prop.sc      = tmp$prop.sc,
                         propsc.mthd  = NULL,
                         form.true    = NULL,
                         w            = tmp$w,            # w is the design matrix for the true propensity score model 
                         mod.insplt   = NULL,
                         num.truc.obs = num.truc.obs) 
      
      # df.noy    = data.used.noy
      # method    = propsc.mthd
      # form.true = propsc.form.true
      
    } else{ # if (propsc.mod.out)
      
      parms.used <- list(trt          = data.used$A,
                         covariates   = data.used.covariates,
                         response     = data.used$Y,
                         prop.sc      = NULL,
                         propsc.mthd  = propsc.mthd,
                         form.true    = propsc.form.true,
                         w            = NULL,          # w is the design matrix for the true propensity score model 
                         mod.insplt   = propsc.mod.insplt,
                         num.truc.obs = num.truc.obs)
    }
    
  } else if (est.used == "G"){
    ulist.used <- list(eval = etemp.g, split = stemp.g, init = itemp)
    
    if (adj.mod.out){
      
      tmp <- est.cond.eff(df        = data.used,
                          method    = adj.mthd,
                          form.true = adj.form.true, 
                          type.var  = type.var)
      
      # df        = data.used
      # method    = adj.mthd
      # form.true = adj.form.true
      
      # if (ncol(tmp$var.rb) < ncol(w)) {
      #   w <- w %>%
      #     select(colnames(tmp$var.rb))
      # }
      
      parms.used <- list(trt            = data.used$A,
                         covariates     = data.used.covariates,
                         response       = data.used$Y,
                         est.cond.eff.1 = tmp$pred.A.1,
                         est.cond.eff.0 = tmp$pred.A.0,
                         var.rb         = tmp$var.rb,
                         adj.mthd       = NULL,
                         form.true      = NULL,
                         w              = tmp$w,            # w is the design matrix for the covariate adjusted model 
                         mod.insplt     = NULL,
                         type.var       = type.var,
                         num.truc.obs   = num.truc.obs) 
      
    } else{ # if (adj.mod.out)
      parms.used <- list(trt            = data.used$A,
                         covariates     = data.used.covariates,
                         response       = data.used$Y,
                         est.cond.eff.1 = NULL,
                         est.cond.eff.0 = NULL,
                         var.rb         = NULL,
                         adj.mthd       = adj.mthd,
                         form.true      = adj.form.true,
                         w              = NULL,            # w is the design matrix for the covariate adjusted model 
                         mod.insplt     = adj.mod.insplt,
                         type.var       = type.var,
                         num.truc.obs   = num.truc.obs) 
    }

    # test split on a continuous covariate
    # x     <- sort(data.used$X1)
    # y     <- order(data.used.wiz.rownumb$X1)
    # wt    <- rep(1, length(y))
    # parms <- parms.used
    
    # test split on a categorical covariate
    # x     <- data.used$X4
    # y     <- data.used.wiz.rownumb$rownumb
    # wt    <- rep(1, length(y))
    # parms <- parms.used
    
  } else {  # "DR"
    ulist.used <- list(eval = etemp.dr, split = stemp.dr, init = itemp)

    # Only implemented both propensity score model and adjustment model 
    # in g-formula fit outside node
    if (propsc.mod.out) {
      
      tmp.propsc <- est.prop.sc(df.noy    = data.used.noy,
                                method    = propsc.mthd,
                                form.true = propsc.form.true)
      
      if (adj.mod.out) {
        tmp.g <- est.cond.eff(df        = data.used,
                              method    = adj.mthd,
                              form.true = adj.form.true,
                              type.var  = type.var)
        
        # Whether have type.var in parms.used depends on adj.mod.out;
        # if out, type.var = NULL; if not, type.var = type.var
        parms.used <- list(trt               = data.used$A,
                           covariates        = data.used.covariates,
                           response          = data.used$Y,
                           prop.sc           = tmp.propsc$prop.sc,
                           propsc.mthd       = NULL,
                           propsc.form.true  = NULL,
                           propsc.mod.insplt = NULL,
                           est.cond.eff.1    = tmp.g$pred.A.1,
                           est.cond.eff.0    = tmp.g$pred.A.0,
                           adj.mthd          = NULL,
                           adj.form.true     = NULL,
                           adj.mod.insplt    = NULL,
                           type.var          = NULL,
                           num.truc.obs      = num.truc.obs)
      } else {
        
        parms.used <- list(trt               = data.used$A,
                           covariates        = data.used.covariates,
                           response          = data.used$Y,
                           prop.sc           = tmp.propsc$prop.sc,
                           propsc.mthd       = NULL,
                           propsc.form.true  = NULL,
                           propsc.mod.insplt = NULL,
                           est.cond.eff.1    = NULL,
                           est.cond.eff.0    = NULL,
                           adj.mthd          = adj.mthd,
                           adj.form.true     = adj.form.true,
                           adj.mod.insplt    = adj.mod.insplt,
                           type.var          = type.var,
                           num.truc.obs      = num.truc.obs)
      } 
      
    } else { # if (propsc.mod.out)
      
      if (adj.mod.out) {
        tmp.g <- est.cond.eff(df        = data.used,
                              method    = adj.mthd,
                              form.true = adj.form.true, 
                              type.var  = type.var)
        
        parms.used <- list(trt               = data.used$A,
                           covariates        = data.used.covariates,
                           response          = data.used$Y,
                           prop.sc           = NULL,
                           propsc.mthd       = propsc.mthd,
                           propsc.form.true  = propsc.form.true,
                           propsc.mod.insplt = propsc.mod.insplt,
                           est.cond.eff.1    = tmp.g$pred.A.1,
                           est.cond.eff.0    = tmp.g$pred.A.0,
                           adj.mthd          = NULL,
                           adj.form.true     = NULL,
                           adj.mod.insplt    = NULL,
                           type.var          = NULL,
                           num.truc.obs      = num.truc.obs)
      } else {
        
        parms.used <- list(trt               = data.used$A,
                           covariates        = data.used.covariates,
                           response          = data.used$Y,
                           prop.sc           = NULL,
                           propsc.mthd       = propsc.mthd,
                           propsc.form.true  = propsc.form.true,
                           propsc.mod.insplt = propsc.mod.insplt,
                           est.cond.eff.1    = NULL,
                           est.cond.eff.0    = NULL,
                           adj.mthd          = adj.mthd,
                           adj.form.true     = adj.form.true,
                           adj.mod.insplt    = adj.mod.insplt,
                           type.var          = type.var,
                           num.truc.obs      = num.truc.obs)
      } 

    } # if (propsc.mod.out) else

  } # else "DR"
  
  # Tree is implemented differently for a binary and a continuous outcome
  # if(type.var == "cont"){
  #   # Fit a large tree using the user written splitting functions
  #   a <- rpart(form.used, 
  #              data     = data.used,
  #              method   = ulist.used,
  #              parms    = parms.used,
  #              control  = rpart.control(cp = 0, minsplit = num.truc.obs * 2, minbucket = min.node, maxsurrogate = 0, maxcompete = 0))
  # } else {
  
    # Fit a large tree using the user written splitting functions
  a <- rpart(form.used, data = data.used.wiz.rownumb,
             method   = ulist.used,
             parms    = parms.used,
             control  = rpart.control(cp = 0, minsplit = num.truc.obs * 2, minbucket = min.node, maxsurrogate = 0, maxcompete = 0))
  # }
  
  if (dim(a$frame)[1] == 1){ # Deal with the root only tree
    
    tree.list = list(a)
    lambda.list = list(Inf)
    return(list(tree.list = tree.list, lambda.list = lambda.list))
    
  } else {
    
    # Finding which variables are leaf nodes
    is.leaf <- (a$frame$var == "<leaf>")
    
    # A function adapted from the partykit package that identifies the rows of the frame 
    # which correspond to child nodes of row i in frame matrix
    rpart.kids <- function(i, is.leaf) {
      if (is.leaf[i]) return(NULL)
      else return(c(i + 1L, 
                    which((cumsum(!is.leaf[-(1L:i)]) + 1L) == cumsum(is.leaf[-(1L:i)]))[1L] + 1L + i))
    }
    
    # Finding goodness of the split
    a$frame$split.stat = 0
    a$frame$split.stat[!is.leaf] = a$splits[, 3]
    
    # Calculating the g(h) parameter for each non-terminal node
    g.h = rep(0, nrow(a$frame))
    for(i in 1:nrow(a$frame)){
      if(is.leaf[i]){g.h[i] = Inf} else{
        # Find all kids of node i
        kids.i = i
        stop.loop = FALSE
        while(stop.loop == FALSE){
          kids.old = kids.i
          for(j in 1:length(kids.i)){
            kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf)))
          }
          if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
          # Calculating g.h for node i
          g.h[i] = sum(a$frame$split.stat[kids.i])/sum(a$frame$split.stat[kids.i] != 0)
        }
      }
    }
    
    # Adding g.h to frame
    a$frame$g.h = g.h
    
    # Start pruning
    # First tree is the large tree
    tree.list = list(a)
    lambda.list = list(0)
    stop.prune = FALSE
    k = 1
    
    while(stop.prune == FALSE){
      tree.used = tree.list[[k]]
      # Calculating the g(h) parameter for each non-terminal node
      tree.used$frame$g.h = rep(0, nrow(tree.used$frame))
      is.leaf.prune <- (tree.used$frame$var == "<leaf>")
      # Setting splitting statistics for new terminal nodes to 0
      tree.used$frame$split.stat[is.leaf.prune] = 0
      
      # Calculating the g(h) function for each non-terminal node
      for(i in 1:nrow(tree.used$frame)){
        if(is.leaf.prune[i]){tree.used$frame$g.h[i] = Inf} else{
          # Find all kids of node i
          kids.i = i
          stop.loop = FALSE
          while(stop.loop == FALSE){
            kids.old = kids.i
            for(j in 1:length(kids.i)){
              kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf.prune)))
            }
            if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
            tree.used$frame$g.h[i] = sum(tree.used$frame$split.stat[kids.i])/sum(tree.used$frame$split.stat[kids.i] != 0)
          }
        }
      }
      
      # Finding the value which minimizes g(h) (among internal nodes)
      to.prune = which.min(tree.used$frame$g.h)
      # Finding the minimum g.h value
      g.h.min = min(tree.used$frame$g.h)
      
      # Find all kids of node to.prune
      kids.i = to.prune
      stop.loop = FALSE
      while(stop.loop == FALSE){
        kids.old = kids.i
        for(j in 1:length(kids.i)){
          kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf.prune)))
        }
        if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
      }
      
      
      # Finding number of splits to prune
      split.to.prune = length(kids.i[which(!is.leaf.prune[kids.i])])
      
      # Creating the new splits and frames for new tree
      splits.new = tree.used$splits[-c(sum(!is.leaf.prune[1:to.prune]):(sum(!is.leaf.prune[1:to.prune]) + split.to.prune - 1)), ]
      frame.new = tree.used$frame[-setdiff(kids.i, to.prune), ]
      
      # Changing all nodes that were internal nodes and are now terminal node to terminal nodes
      frame.new$var[to.prune] =  "<leaf>"
      
      tree.new = tree.used
      tree.new$frame = frame.new
      if(class(splits.new) == "matrix"){tree.new$splits = splits.new}
      if(class(splits.new) == "numeric"){
        tree.new$splits = matrix(splits.new, nrow = 1)
        colnames(tree.new$splits) = colnames(tree.used$splits)
      }
      
      # Changing the terminal node for $where in rpart object
      tree.new$where = tree.used$where
      tree.new$where[tree.new$where %in% kids.i] = to.prune
      tree.new$where[tree.new$where > max(kids.i)] = tree.new$where[tree.new$where > max(kids.i)] - length(kids.i) + 1
      tree.new$where = as.integer(tree.new$where)
      
      k = k+1
      # Add tree and lambda to the list
      tree.list[[k]] <- tree.new
      lambda.list[[k]] <- g.h.min
      if(sum(tree.new$frame$var == "<leaf>") == 1){stop.prune = TRUE}  
    }
    return(list(tree.list = tree.list, lambda.list = lambda.list))
  }
  
}
