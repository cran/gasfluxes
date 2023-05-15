#' @title Robust linear concentration - time model
#' 
#' @description
#'   Fit a linear model to concentration - time data using robust methods.
#' @param t time values (usually in hours)
#' @param C concentration values
#' @param A area covered by the chamber
#' @param V effective volume of the chamber
#' @param serie id of the flux measurement
#' @param verbose logical, TRUE prints message after each flux calculation
#' @param plot logical, mainly intended for use in \code{\link{gasfluxes}}
#' @param \dots further parameters, currently none
#'  
#' @return
#'  A list of
#'    \item{f0}{flux estimate}
#'    \item{f0.se}{standard error of flux estimate}
#'    \item{f0.p}{p-value of flux estimate}
#'    \item{C0}{estimated concentration at t = 0 (intercept)}  
#'    \item{weights}{robustness weights}
#'    \item{diagnostics}{error or warning messages}
#'    
#' @details
#' This is basically a wrapper of \code{\link{rlm}} using the Huber M estimator. This function never weights the first or last time point with zero with very few data points. However, there might exist "better" robust regression methods for flux estimation.
#'      
#' @examples
#' #a single fit
#' t <- c(0, 1/3, 2/3, 1)
#' C <- c(320, 330, 315, 351)
#' print(fit <- rlin.fit(t, C, 1, 0.3, "a"))
#' plot(C ~ t)
#' curve({fit$f0/0.3 * x + fit$C0}, from = 0, to = 1, add = TRUE)

#' 
#' @importFrom MASS rlm
#' @importFrom sfsmisc f.robftest
#' @importFrom graphics curve
#' @importFrom stats predict 
#' 
#' @export

rlin.fit <- function (t, C, A = 1, V, serie = "", verbose = TRUE, plot = FALSE, ...) {
  tryCatch({
    stopifnot(length(t) > 3)
    fit <- withWarnings(rlm(C ~ t, maxit = 200))
    w <- if (is.null(fit$warnings)) "" else fit$warnings[[1]]$message
    fit <- fit$value
    fitsum <- withWarnings(summary(fit))
    if (w == "") w <- if (is.null(fitsum$warnings)) "" else fitsum$warnings[[1]]$message
    fitsum <- fitsum$value
    fitsumCoef <- fitsum$coef
    try({
      if (plot) {
        curve(predict(fit, newdata = data.frame(t = x)), 
              from = min(t), to = max(t), add = TRUE, col = "green", lty = 2)
      }}, silent = TRUE)
    
    res <- list(
      f0 = fitsumCoef["t", "Value"] * V/A, 
      f0.se = fitsumCoef["t", "Std. Error"] * V/A, 
      f0.p = f.robftest(fit)[["p.value"]], 
      C0 = fitsumCoef["(Intercept)", "Value"],
      weights = fit$w,
      diagnostics = w)
    if (verbose) message(serie, if (w == "") ": rlm fit successful" else ": rlm fit warning")
    res
  },
  error = function(cond) {
    if (verbose) message(serie, ": rlm fit failed")
    list(
      f0 = NA_real_, 
      f0.se = NA_real_, 
      f0.p = NA_real_, 
      C0=NA_real_,
      weights=NA_real_,
      diagnostics=cond$message)
  }
  )  
} 

utils::globalVariables("x") #non-standard evaluation in curve