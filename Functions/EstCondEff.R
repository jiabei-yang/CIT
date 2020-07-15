# Catch both the results and the warning message from expr
withWarnings <- function(expr) {
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(expr, warning = wHandler)
  list(value = val, warnings = myWarnings)
} 

est.cond.eff <- function(df, method = "RF", form.true = NULL, type.var = "cont"){
  
  # function used when the outcome is continuous.
  # df:     a data frame with only the ordered columns A, Y, X, where A, Y are the first 2 columns 
  #         and other covariates start from the 3rd column.
  # method: one of "RF", "GLM", or "GAM". The default is "RF".
  
  # Need to move into is.na(form.true) since is.na() cause problem when the parmaeter is formula
  # if (!is.null(form.true)){
  #   form.true <- as.formula(form.true)
  # }
  
  if (method == "RF"){
    
    if (is.null(form.true)){
      fit <- rfsrc(Y ~ ., data = df)
    } else {
      form.true <- as.formula(form.true)
      fit <- rfsrc(form.true, data = df)
    }
    
    df.A.1   <- df
    df.A.1$A <- 1
    pred.A.1 <- predict(fit, newdata = df.A.1)$predicted
    
    df.A.0   <- df
    df.A.0$A <- 0
    pred.A.0 = predict(fit, newdata = df.A.0)$predicted
    
  } else{
    
    if (method == "GLM"){
      
      if (is.null(form.true)) {
        
        form.glm <- "Y ~ A"
        for (i in 3:dim(df)[2]){
          form.glm <- paste(form.glm, colnames(df)[i], sep = " + ")
          form.glm <- paste(form.glm, paste("A:", colnames(df)[i], sep = ""), sep = " + ")
        }
        form.glm <- as.formula(form.glm)
        
        if (type.var == "cont") {
          fit <- glm(form.glm, data = df, family = "gaussian")
        } else if (type.var == "bin") {
          
          # fit <- ressWarnings(glm(form.glm, data = df, family = binomial(link = "logit")))
          fit <- glm(form.glm, data = df, family = binomial(link = "logit"))
          
          # Warning message if there are fitted probability 0 / 1
          # tmp.fit <- tryCatch(glm(form.glm, data = df, family = binomial(link = "logit")), warning = function(w) w)
          # if(is(tmp.fit, "warning")) {
          #   print(df)
          #   print(summary(fit))
          #   print(predict(fit, type = "response"))
          # }
          
        } # Currently only binary and continuous outcome will work
        w   <- model.matrix(form.glm, data = df)[, !is.na(fit$coefficients)] 
        
      } else {     # if (is.null(form.true)) else
        
        form.true <- as.formula(form.true)
        
        if (type.var == "cont") {
          fit <- glm(form.true, data = df, family = "gaussian")
        } else if (type.var == "bin") {
          fit <- glm(form.true, data = df, family = binomial(link = "logit"))
        } # Currently only binary and continuous outcome will work
        w   <- model.matrix(form.true, data = df)[, !is.na(fit$coefficients)]
      }
      
    } else { # "GAM"
      
      if (is.null(form.true)){
        form.gam <- "Y ~ A"
        
        for (i in 3:dim(df)[2]){
          if (class(df[, i]) == "numeric") {
            form.gam <- paste(form.gam, paste("s(", colnames(df)[i], ")", sep = ""), sep = "+")
          } else {
            form.gam <- paste(form.gam, colnames(df)[i], sep = "+")
          }
        }
        form.gam <- as.formula(form.gam)
        
        if (type.var == "cont") {
          fit <- gam(form.gam, family = gaussian(link = identity), data = df)
        } else if (type.var == "bin") {
          fit <- gam(form.gam, family = binomial(link = "logit"), data = df)
        }
        
      } else {
        
        form.true <- as.formula(form.true)
        if (type.var == "cont") {
          fit <- gam(form.true, family = gaussian(link = identity), data = df)
        } else if (type.var == "bin") {
          fit <- gam(form.true, family = binomial(link = "logit"), data = df)
        }
        
      }
      
    } # end else "GAM"
    
    # var.rb = try(plm::vcovHC(fit, type = "HC"))
    var.rb = plm::vcovHC(fit, type = "HC")
    
    # if ("try-error" %in% class(var.rb)) {
    #   print(summary(fit))
    #   print(fit)
    #   print(df)
    #   return(NULL)
    # }
    
    df.A.1   <- df
    df.A.1$A <- 1
    # Prediction for continuous outcome is the same when type is "response" or "link"
    
    # Warning message if there are interaction term and the main effect term has rank deficiency
    # tmp.pred <- tryCatch(predict(fit, newdata = df.A.1, type = "response"), warning = function(w) w)
    # if(is(tmp.pred,"warning")) {
    #   print(df)
    #   print(form.glm)
    #   print(summary(fit))
    #   print(predict(fit, newdata = df.A.1, type = "response"))
    # }
    # pred.A.1 = ressWarnings(predict(fit, newdata = df.A.1, type = "response"))
    pred.A.1 = predict(fit, newdata = df.A.1, type = "response")
    
        
    df.A.0   <- df
    df.A.0$A <- 0
    # pred.A.0 = suppressWarnings(predict(fit, newdata = df.A.0, type = "response"))
    pred.A.0 = predict(fit, newdata = df.A.0, type = "response")
    
  }  
  
  if (method == "GLM"){
    return(list(pred.A.0 = pred.A.0,
                pred.A.1 = pred.A.1, 
                var.rb   = var.rb,
                w        = data.frame(w)))
  } else {
    return(list(pred.A.0 = pred.A.0,
                pred.A.1 = pred.A.1))
  }
  
}

est.prop.sc <- function(df.noy, method = "RF", form.true = NULL){
  
  # Function to predict the propensity scores
  # df:     a data frame with only the ordered columns A, Y, X, where A, Y are the first 2 columns 
  #         and other covariates start from the 3rd column.
  # method: one of "RF", "GLM", or "GAM". The default is "RF".
  # form.true: formula for the true model; if "GAM" used for the method, need to include s() functions.
  
  # Need to move into is.na(form.true) since is.na() cause problem when the parmaeter is formula
  # if (!is.null(form.true)){
  #   form.true <- as.formula(form.true)
  # }
  
  if (method == "RF"){
    if (is.null(form.true)){
      fit <- rfsrc(A ~ ., data = df.noy)
    } else {
      form.true <- as.formula(form.true)
      fit <- rfsrc(form.true, data = df.noy)
    }
    prop.sc <- predict(fit)$predicted
    
  } else{
    
    if (method == "GLM"){
      if (is.null(form.true)){
        
        # fit <- suppressWarnings(glm(A ~ ., data = df.noy, family = binomial(link = "logit")))
        fit <- glm(A ~ ., data = df.noy, family = binomial(link = "logit"))
        # Warning message if there are fitted probability 0 / 1
        # tmp.fit <- tryCatch(glm(A ~ ., data = df.noy, family = binomial(link = "logit")), warning = function(w) w)
        # if(is(tmp.fit, "warning")) {
        #   print(df.noy)
        #   print(summary(fit))
        #   print(predict(fit, type = "response"))
        # }
        
        w   <- model.matrix(A ~ .,
                            data = df.noy)[, !is.na(fit$coefficients)] 
      } else {
        
        form.true <- as.formula(form.true)
        # fit <- suppressWarnings(glm(form.true, data = df.noy, family = binomial(link = "logit")))
        fit <- glm(form.true, data = df.noy, family = binomial(link = "logit"))
        w   <- model.matrix(form.true, 
                            data = df.noy)[, !is.na(fit$coefficients)]
      }
      
    } else {
      
      if (is.null(form.true)){
        form.gam <- "A ~ "
        
        if (class(df.noy[, 2]) == "numeric"){
          form.gam <- paste(form.gam, paste("s(", colnames(df.noy)[2], ")", sep = ""), sep = "")
        } else {
          form.gam <- paste(form.gam, colnames(df.noy)[2], sep = "")
        }
        
        for (i in 3:dim(df.noy)[2]){
          if (class(df.noy[, i]) == "numeric") {
            form.gam <- paste(form.gam, paste("s(", colnames(df.noy)[i], ")", sep = ""), sep = "+")
          } else {
            form.gam <- paste(form.gam, colnames(df.noy)[i], sep = "+")
          }
        }
        form.gam <- as.formula(form.gam)
        fit <- gam(form.gam, family = binomial(link = "logit"), data = df.noy)
      } else {
        form.true <- as.formula(form.true)
        fit <- gam(form.true, family = binomial(link = "logit"), data = df.noy)
      }
      
    }
    
    prop.sc = predict(fit, type = "response")
  }  
  
  prop.sc[prop.sc < 0.1] <- 0.1
  prop.sc[prop.sc > 0.9] <- 0.9
  
  if (method == "GLM"){
    return(list(prop.sc = as.numeric(prop.sc),
                w       = data.frame(w)))
  } else {
    return(list(prop.sc = as.numeric(prop.sc)))
  }
  
}

