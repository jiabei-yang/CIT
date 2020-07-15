create.sequence = function(data.used, est.used, type.var = "cont",        # For ipw, type.var does not matter
                           propsc.mod.loc, propsc.mthd = "GLM", propsc.form.true = NULL, # params for ipw, dr
                           adj.mod.loc, adj.mthd = "GLM", adj.form.true = NULL,          # params for g, dr
                           num.truc.obs = 30, min.node = 20){
  
  # This function calculates the sequence of candidate trees.
  # Input: data.used the dataset with treatment as first column, and outcome as second column
  # est.used: c("IPW", "G", "DR")
  # need.cond.exp = a logical statement if an outside estimator for the conditional expectation
  # is needed
  # type.var = "cont" for continuous outcomes, "bin" for binary outcomes
  # Output: A list with two elements tree.list which is the sequence of candidate trees 
  # and lambda.list which is the sequence of penalization parameters
  # propsc.mod.loc: c("out", "node", "split")
  # adj.mod.loc: c("out", "node", "split")
  
  # All fitting used rownumb as the outcome regardless of type.var
  data.used.covariates  <- data.used[!colnames(data.used) %in% c("Y", "A")]
  form.used             <- as.formula(paste("rownumb ~ ", paste(colnames(data.used.covariates), collapse= "+")))
  data.used.wiz.rownumb <- cbind(rownumb = 1:dim(data.used)[1], data.used)
  data.used.noy         <- data.used[, !colnames(data.used) %in% c("Y")]
  
  # Creating ulist.used and parameter vector used
  if (est.used == "IPW"){
    
    ulist.used <- list(eval = etemp.ipw, split = stemp.ipw, init = itemp)
    # ulist.used <- list(eval = etemp.ipw.whl, split = stemp.ipw.whl, init = itemp)
    
    # Calculate the avg.trt.eff for the whole dataset
    # for those terminal node when there are only treated/untreated units in the node
    whl.propsc <- gen.fullrank.ipw(df.noy           = data.used.noy, 
                                   propsc.form.true = propsc.form.true)
    whl.propsc <- withWarnings(est.prop.sc(df.noy    = whl.propsc$df.noy.fullrank,
                                           method    = propsc.mthd,
                                           form.true = whl.propsc$propsc.form.true.updated))
    mu.1 <- mean((data.used$Y * data.used$A/ whl.propsc$value$prop.sc))
    mu.0 <- mean((data.used$Y * (1 - data.used$A)/ (1 - whl.propsc$value$prop.sc)))
    avg.trt.effct <- mu.1 - mu.0
    
    parms.used <- list(trt            = data.used$A,
                       covariates     = data.used.covariates,
                       response       = data.used$Y,
                       whl.propsc     = whl.propsc$value$prop.sc, 
                       propsc.mod.loc = propsc.mod.loc,
                       propsc.mthd    = propsc.mthd,                     # need these 2 parameters since need to find w if singular within node
                       form.true      = propsc.form.true,
                       # w             = whl.propsc$value$w,            # w is the design matrix for the true propensity score model 
                       avg.trt.effct  = avg.trt.effct,                   # avg.trt.effct is needed for calculating order of categories when there is only treated/untreated unit
                       num.truc.obs   = num.truc.obs) 
    
    # if (propsc.mod.out) {
    #   
    #   # tmp <- est.prop.sc(df.noy    = data.used.noy,
    #   #                    method    = propsc.mthd,
    #   #                    form.true = propsc.form.true)
    #   parms.used <- list(trt           = data.used$A,
    #                      covariates    = data.used.covariates,
    #                      response      = data.used$Y,
    #                      prop.sc       = whl.propsc$value$prop.sc,
    #                      propsc.mthd   = propsc.mthd,                     # need these 2 parameters since need to find w if singular within node
    #                      form.true     = propsc.form.true,
    #                      # w             = whl.propsc$value$w,            # w is the design matrix for the true propensity score model 
    #                      mod.insplt    = NULL,
    #                      whl.propsc    = NULL,                            # do not need whl.propsc since we have prop.sc already
    #                      avg.trt.effct = avg.trt.effct,                   # avg.trt.effct is needed for calculating order of categories when there is only treated/untreated unit
    #                      num.truc.obs  = num.truc.obs) 
    #   
    #   # df.noy    = data.used.noy
    #   # method    = propsc.mthd
    #   # form.true = propsc.form.true
    #   
    # } else { # if (propsc.mod.out)
    #   
    #   parms.used <- list(trt           = data.used$A,
    #                      covariates    = data.used.covariates,
    #                      response      = data.used$Y,
    #                      prop.sc       = NULL,
    #                      propsc.mthd   = propsc.mthd,
    #                      form.true     = propsc.form.true,
    #                      # w             = NULL,          # w is the design matrix for the true propensity score model 
    #                      mod.insplt    = propsc.mod.insplt,
    #                      whl.propsc    = whl.propsc$value$prop.sc,
    #                      avg.trt.effct = avg.trt.effct,
    #                      num.truc.obs  = num.truc.obs)
    # }
    
  } else if (est.used == "G"){
    ulist.used <- list(eval = etemp.g, split = stemp.g, init = itemp)
    
    # when there is warning: use est.cond.eff.1/0 for imputing predicted outcome
    # when there is only treated/untreated unit use average treatment effect
    whl.g <- gen.fullrank.g(df            = data.used,
                            adj.form.true = adj.form.true)
    whl.g <- withWarnings(est.cond.eff(df        = whl.g$df.fullrank,
                                       method    = adj.mthd,
                                       form.true = whl.g$adj.form.true.updated,
                                       type.var  = type.var))
    avg.trt.effct <- mean(whl.g$value$pred.A.1) - mean(whl.g$value$pred.A.0)
    
    parms.used <- list(trt                = data.used$A,
                       covariates         = data.used.covariates,
                       response           = data.used$Y,
                       whl.est.cond.eff.1 = whl.g$value$pred.A.1,
                       whl.est.cond.eff.0 = whl.g$value$pred.A.0,
                       whl.var.rb         = whl.g$value$var.rb,
                       adj.mod.loc        = adj.mod.loc,
                       adj.mthd           = adj.mthd,
                       form.true          = adj.form.true,
                       whl.w              = whl.g$value$w,            # w is the design matrix for the covariate adjusted model 
                       type.var           = type.var,
                       avg.trt.effct      = avg.trt.effct,
                       num.truc.obs       = num.truc.obs) 
    
    # if (adj.mod.out) {
    #   
    #   # tmp <- est.cond.eff(df        = data.used,
    #   #                     method    = adj.mthd,
    #   #                     form.true = adj.form.true, 
    #   #                     type.var  = type.var)
    #   # df        = data.used
    #   # method    = adj.mthd
    #   # form.true = adj.form.true
    #   # type.var  = type.var
    #   
    #   # if (ncol(tmp$var.rb) < ncol(w)) {
    #   #   w <- w %>%
    #   #     select(colnames(tmp$var.rb))
    #   # }
    #   
    #   parms.used <- list(trt            = data.used$A,
    #                      covariates     = data.used.covariates,
    #                      response       = data.used$Y,
    #                      est.cond.eff.1 = whl.g$value$pred.A.1,
    #                      est.cond.eff.0 = whl.g$value$pred.A.0,
    #                      var.rb         = whl.g$value$var.rb,
    #                      adj.mthd       = NULL,
    #                      form.true      = NULL,
    #                      w              = whl.g$value$w,            # w is the design matrix for the covariate adjusted model 
    #                      mod.insplt     = NULL,
    #                      type.var       = type.var,
    #                      avg.trt.effct  = avg.trt.effct,
    #                      num.truc.obs   = num.truc.obs) 
    #   
    # } else{ # if (adj.mod.out)
    #   parms.used <- list(trt            = data.used$A,
    #                      covariates     = data.used.covariates,
    #                      response       = data.used$Y,
    #                      est.cond.eff.1 = NULL,
    #                      est.cond.eff.0 = NULL,
    #                      var.rb         = NULL,
    #                      adj.mthd       = adj.mthd,
    #                      form.true      = adj.form.true,
    #                      w              = NULL,            # w is the design matrix for the covariate adjusted model 
    #                      mod.insplt     = adj.mod.insplt,
    #                      type.var       = type.var,
    #                      avg.trt.effct  = avg.trt.effct,
    #                      num.truc.obs   = num.truc.obs) 
    # }
    
    # test split on a continuous covariate
    # x     <- sort(data.used$X1)
    # y     <- order(data.used$X1)
    # wt    <- rep(1, length(y))
    # parms <- parms.used
    
    # test split on a continuous covariate, in a branch
    # x     <- sort(data.used$X1[data.used$X4 %in% c("B", "D")])
    # y     <- order(data.used$X1)[which(data.used$X4[order(data.used$X1)] %in% c("B", "D"))]
    # wt    <- rep(1, length(y))
    # parms <- parms.used
    
    # x     <- sort(data.used$surv2md1[data.used$imputed.cat2 %in% c("NA")])
    # y     <- order(data.used$surv2md1)[which(data.used$imputed.cat2[order(data.used$surv2md1)] %in% c("NA"))]
    # wt    <- rep(1, length(y))
    # parms <- parms.used
    
    # test split on a categorical covariate
    # x     <- data.used$X4
    # y     <- data.used.wiz.rownumb$rownumb
    # wt    <- rep(1, length(y))
    # parms <- parms.used
    
    # test split on a categorical covariate, in a deep branch
    # x     <- data.used$X4[data.used$X4 %in% c("B", "D")]
    # y     <- data.used.wiz.rownumb$rownumb[data.used$X4 %in% c("B", "D")]
    # wt    <- rep(1, length(y))
    # parms <- parms.used
    
    x     <- data.used$imputed.cat2[data.used$imputed.cat2 %in% c("NA", "Colon Cancer", "MOSF w/Malignancy")]
    y     <- data.used.wiz.rownumb$rownumb[data.used$imputed.cat2 %in% c("NA", "Colon Cancer", "MOSF w/Malignancy")]
    wt    <- rep(1, length(y))
    parms <- parms.used
    
    # x     <- data.used$sex
    # y     <- data.used.wiz.rownumb$rownumb
    # wt    <- rep(1, length(y))
    # parms <- parms.used
    
  } else {  # "DR"
    ulist.used <- list(eval = etemp.dr, split = stemp.dr, init = itemp)
    
    # Calculate the avg.trt.eff for the whole dataset
    # for those terminal node when there are only treated/untreated units in the node
    whl.propsc <- gen.fullrank.ipw(df.noy           = data.used.noy, 
                                   propsc.form.true = propsc.form.true)
    whl.propsc <- withWarnings(est.prop.sc(df.noy    = whl.propsc$df.noy.fullrank,
                                           method    = propsc.mthd,
                                           form.true = whl.propsc$propsc.form.true.updated))
    
    whl.g <- gen.fullrank.g(df            = data.used,
                            adj.form.true = adj.form.true)
    whl.g <- withWarnings(est.cond.eff(df        = whl.g$df.fullrank,
                                       method    = adj.mthd,
                                       form.true = whl.g$adj.form.true.updated, 
                                       type.var  = type.var))
    
    mu.1 <- mean(data.used$A * (data.used$Y - whl.g$value$pred.A.1) / whl.propsc$value$prop.sc + whl.g$value$pred.A.1)
    mu.0 <- mean((1 - data.used$A) * (data.used$Y - whl.g$value$pred.A.0) / (1 - whl.propsc$value$prop.sc) + whl.g$value$pred.A.0)
    avg.trt.effct <- mu.1 - mu.0
    
    parms.used <- list(trt                = data.used$A,
                       covariates         = data.used.covariates,
                       response           = data.used$Y,
                       whl.propsc         = whl.propsc$value$prop.sc,
                       propsc.mod.loc     = propsc.mod.loc,
                       propsc.mthd        = propsc.mthd,
                       propsc.form.true   = propsc.form.true,
                       whl.est.cond.eff.1 = whl.g$value$pred.A.1,
                       whl.est.cond.eff.0 = whl.g$value$pred.A.0,
                       adj.mod.loc        = adj.mod.loc,
                       adj.mthd           = adj.mthd,
                       adj.form.true      = adj.form.true,
                       type.var           = type.var,
                       avg.trt.effct      = avg.trt.effct,
                       num.truc.obs       = num.truc.obs)
    
    # Only implemented both propensity score model and adjustment model 
    # in g-formula fit outside node
    # if (propsc.mod.out) {
    #   
    #   tmp.propsc <- est.prop.sc(df.noy    = data.used.noy,
    #                             method    = propsc.mthd,
    #                             form.true = propsc.form.true)
    #   
    #   if (adj.mod.out) {
    #     tmp.g <- est.cond.eff(df        = data.used,
    #                           method    = adj.mthd,
    #                           form.true = adj.form.true,
    #                           type.var  = type.var)
    #     
    #     # Whether have type.var in parms.used depends on adj.mod.out;
    #     # if out, type.var = NULL; if not, type.var = type.var
    #     parms.used <- list(trt               = data.used$A,
    #                        covariates        = data.used.covariates,
    #                        response          = data.used$Y,
    #                        prop.sc           = tmp.propsc$prop.sc,
    #                        propsc.mthd       = NULL,
    #                        propsc.form.true  = NULL,
    #                        propsc.mod.insplt = NULL,
    #                        est.cond.eff.1    = tmp.g$pred.A.1,
    #                        est.cond.eff.0    = tmp.g$pred.A.0,
    #                        adj.mthd          = NULL,
    #                        adj.form.true     = NULL,
    #                        adj.mod.insplt    = NULL,
    #                        type.var          = NULL,
    #                        avg.trt.effct     = NULL,
    #                        num.truc.obs      = num.truc.obs)
    #   } else {
    #     
    #     parms.used <- list(trt               = data.used$A,
    #                        covariates        = data.used.covariates,
    #                        response          = data.used$Y,
    #                        prop.sc           = tmp.propsc$prop.sc,
    #                        propsc.mthd       = NULL,
    #                        propsc.form.true  = NULL,
    #                        propsc.mod.insplt = NULL,
    #                        est.cond.eff.1    = NULL,
    #                        est.cond.eff.0    = NULL,
    #                        adj.mthd          = adj.mthd,
    #                        adj.form.true     = adj.form.true,
    #                        adj.mod.insplt    = adj.mod.insplt,
    #                        type.var          = type.var,
    #                        avg.trt.effct     = avg.trt.effct,
    #                        num.truc.obs      = num.truc.obs)
    #   } 
    #   
    # } else { # if (propsc.mod.out)
    #   
    #   if (adj.mod.out) {
    #     tmp.g <- est.cond.eff(df        = data.used,
    #                           method    = adj.mthd,
    #                           form.true = adj.form.true, 
    #                           type.var  = type.var)
    #     
    #     parms.used <- list(trt               = data.used$A,
    #                        covariates        = data.used.covariates,
    #                        response          = data.used$Y,
    #                        prop.sc           = NULL,
    #                        propsc.mthd       = propsc.mthd,
    #                        propsc.form.true  = propsc.form.true,
    #                        propsc.mod.insplt = propsc.mod.insplt,
    #                        est.cond.eff.1    = tmp.g$pred.A.1,
    #                        est.cond.eff.0    = tmp.g$pred.A.0,
    #                        adj.mthd          = NULL,
    #                        adj.form.true     = NULL,
    #                        adj.mod.insplt    = NULL,
    #                        type.var          = NULL,
    #                        avg.trt.effct     = avg.trt.effct,
    #                        num.truc.obs      = num.truc.obs)
    #   } else {
    #     
    #     parms.used <- list(trt               = data.used$A,
    #                        covariates        = data.used.covariates,
    #                        response          = data.used$Y,
    #                        prop.sc           = NULL,
    #                        propsc.mthd       = propsc.mthd,
    #                        propsc.form.true  = propsc.form.true,
    #                        propsc.mod.insplt = propsc.mod.insplt,
    #                        est.cond.eff.1    = NULL,
    #                        est.cond.eff.0    = NULL,
    #                        adj.mthd          = adj.mthd,
    #                        adj.form.true     = adj.form.true,
    #                        adj.mod.insplt    = adj.mod.insplt,
    #                        type.var          = type.var,
    #                        avg.trt.effct     = avg.trt.effct,
    #                        num.truc.obs      = num.truc.obs)
    #   } 
    #   
    # } # if (propsc.mod.out) else
    
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
