itemp <- function(y, offset, parms, wt) {
  
  sfun <- function(yval, dev, wt, ylevel, digits ) {
    paste(" mean=", format(signif(yval, digits)),
          ", MSE=" , format(signif(dev/wt, digits)),
          sep = '')
  }
  
  environment(sfun) <- .GlobalEnv
  list(y = c(y), parms = parms, numresp = 1, numy = 1, summary = sfun)
}
