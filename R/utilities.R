#https://stackoverflow.com/a/4947528/1412059

withWarnings <- function(expr) {
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(expr, warning = wHandler)
  list(value = val, warnings = myWarnings)
} 


#' @importFrom stats logLik
aicc <- function(fit) {
  ll <- logLik(fit)
  n <- attr(ll, "nobs")
  if (n < 3) return(NA_real_)
  m <- attr(ll, "df")
  if (n == m) return(NA_real_)
  AIC <- -2 * as.numeric(ll) + 2 * m
  AICc <- AIC + (2 * m ^ 2 + 2 * m) / (n - m - 1)
  AICc
}

