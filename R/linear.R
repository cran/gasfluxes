#' @title Linear concentration - time model
#'
#' @description
#'   Fit a linear model to concentration - time data.
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
#'    \item{AIC}{Akaike information criterion} 
#'    \item{AICc}{Akaike information criterion with small sample correction}
#'    \item{RSE}{residual standard error (sigma from summary.nls)}
#'    \item{r}{Pearson's correlation coefficient}
#'    \item{diagnostics}{error or warning messages}
#'    
#' @details
#' This is basically a wrapper of R's OLS fitting facilities. For now \code{\link{lm}} (and methods for objects of class "lm") is used, 
#' but this may change to more efficient alternatives in later versions. 
#'      
#' @examples
#' #a single fit
#' t <- c(0, 1/3, 2/3, 1)
#' C <- c(320, 341, 352, 359)
#' print(fit <- lin.fit(t, C, 1, 0.3, "a"))
#' plot(C ~ t)
#' curve({fit$f0/0.3 * x + fit$C0}, from = 0, to = 1, add = TRUE)
#' 
#' @importFrom stats lm predict AIC cor
#' @importFrom graphics lines
#' @export


lin.fit <- function (t, C, A = 1, V, serie = "", verbose = TRUE, plot = FALSE, ...) {
  tryCatch({
    stopifnot(length(t) > 2)
    fit <- withWarnings(lm(C ~ t))
    w <- if (is.null(fit$warnings)) "" else fit$warnings[[1]]$message
    fit <- fit$value
    r <- cor(C, t)
    fitsum <- withWarnings(summary(fit))
    if (w == "") w <- if (is.null(fitsum$warnings)) "" else fitsum$warnings[[1]]$message
    fitsum <- fitsum$value
    fitsumCoef <- fitsum$coef
    try({
      if (plot) {
        lines(t, predict(fit), col = "black")
      }
    }, silent = TRUE)
        
    res <- list(
      f0 = fitsumCoef["t", "Estimate"] * V/A, 
      f0.se = fitsumCoef["t", "Std. Error"] * V/A, 
      f0.p = fitsumCoef["t", "Pr(>|t|)"], 
      C0 = fitsumCoef["(Intercept)", "Estimate"],
      AIC = AIC(fit),
      AICc = aicc(fit),
      RSE = fitsum$sigma,
      r = r,
      diagnostics = w)
    if (verbose) message(serie, if (w == "") ": lm fit successful" else ": lm fit warning")
    res
  },
  error = function(cond) {
    if (verbose) message(serie, ": linear fit failed")
    list(
      f0 = NA_real_, 
      f0.se = NA_real_, 
      f0.p = NA_real_, 
      C0 = NA_real_,
      AIC = NA_real_,
      AICc = NA_real_,
      RSE = NA_real_,
      r = NA_real_,
      diagnostics=cond$message)
  }
  )  
} 