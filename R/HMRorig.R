#' @title HMR fit using the Pedersen fitting algorithm 
#'   
#' @description
#'   Fit the HMR model using the algorithm from the HMR package.
#' @param t time values (usually in hours)
#' @param C concentration values
#' @param A area covered by the chamber
#' @param V effective volume of the chamber
#' @param serie id of the flux measurement
#' @param verbose logical, TRUE prints message after each flux calculation
#' @param plot logical, mainly intended for use in \code{\link{gasfluxes}}
#' @param ngrid see the HMR documentation
#' @param \dots further parameters, currently none
#'  
#' @return
#'  A list of
#'    \item{f0}{flux estimate}
#'    \item{f0.se}{standard error of flux estimate}
#'    \item{f0.p}{p-value of flux estimate}
#'    \item{kappa, phi}{other parameters of the HMR model}
#'    \item{AIC}{Akaike information criterion}
#'    \item{AICc}{Akaike information criterion with small sample correction}
#'    \item{diagnostics}{error or warning messages}
#'    
#' @details
#' The HMR model (Pedersen et al., 2010) is \eqn{C(t)=\phi+f_0 \frac{e^{-\kappa t}}{-\kappa V/A}}{C(t) = \phi +f0 exp(-\kappa t)/(-\kappa V/A)}.
#' The algorthm from the HMR package version 0.3.1 is used for fitting. Note that this is very inefficient and standard errors and 
#' p-values are over-estimated. \code{\link{HMR.fit}} is recommended instead and this function is only provided to be able to reproduce results obtained with
#' older versions of the HMR package.
#'
#' @author 
#' Asger R. Pedersen for code copied from the not-exported HMR:::.HMR.fit1 function,
#' Roland Fuss 
#'
#' @references
#' Pedersen, A.R., Petersen, S.O., Schelde, K., 2010. A comprehensive approach to soil-atmosphere trace-gas flux estimation with static chambers. European Journal of Soil Science 61(6), 888-902.
#'  
#' @examples
#' #a single fit
#' t <- c(0, 1/3, 2/3, 1)
#' C <- c(320, 341, 352, 359)
#' print(fit <- HMR.orig(t, C, 1, 0.3, "a"))
#' plot(C ~ t)
#' curve({fit$phi + fit$f0 * exp(-fit$kappa * x)/(-fit$kappa*0.3)}, 
#'       from = 0, to = 1, add = TRUE)
#' 
#' @importFrom AICcmodavg AICcCustom
#' @export

HMR.orig <- function (t, C, A = 1, V, serie = "", 
                     verbose = TRUE, ngrid = 1000, plot = FALSE, ...) {
  
  tryCatch({
    stopifnot(length(t) > 3)
    fit <- withWarnings(.HMR.fit1(t, C ,A, V, serie,
                     ngrid, LR.always = FALSE,
                     FollowHMR = TRUE,JPG = FALSE, PS = FALSE,
                     PHMR = FALSE,npred = 500, 
                     xtxt = "", ytxt = "", pcttxt ="",
                     MSE.zero = 10*max(.Machine$double.eps,.Machine$double.neg.eps),
                     bracketing.tol = 1e-7,
                     bracketing.maxiter = 1000))
    w <- if (is.null(fit$warnings)) "" else fit$warnings[[1]]$message
    fit <- fit$value
    
    HMR.fun <-  function (phi,kappa,f0,t, V, A) {
      phi+f0*exp(-kappa*t)/(-kappa*V/A)
    }
    SSE <- sum((C-HMR.fun(fit[["phi"]],
                          fit[["kappa"]],
                          fit[["f0"]],
                          t, V, A))^2)
    n <- length(t)
    k <- 3
    logLik <- -n * (log(2 * pi) + 1 - log(n) +  log(SSE))/2
    
    AIC_HMR <- AICcCustom(logLik, K = k + 1, second.ord = FALSE, nobs = n)
    AICc_HMR <- AICcCustom(logLik, K = k + 1, second.ord = TRUE, nobs = n)
    sigma <- sqrt(SSE/(n-k))
    try({
      if (plot) {
        curve(HMR.fun(fit[["phi"]], fit[["kappa"]], fit[["f0"]], x, V, A), 
              from = min(t), to = max(t), add = TRUE, col = "yellow", lty = 2)
      }}, silent = TRUE)
    
    
    res <- list(
      f0 = fit[["f0"]], 
      f0.se = fit[["f0.se"]], 
      f0.p =  fit[["f0.p"]], 
      kappa = fit[["kappa"]],
      phi = fit[["phi"]],
      AIC = AIC_HMR,
      AICc = AICc_HMR,
      RSE = sigma,
      diagnostics = w)
    if (verbose) message(serie, if (w == "") ": (orig.) HMR fit successful" else ": (orig.) HMR fit warning")
    res
  },
  error = function(cond) {
    if (verbose) message(serie, ": (orig.) HMR fit failed")
    list(
      f0 = NA_real_, 
      f0.se = NA_real_, 
      f0.p = NA_real_, 
      kappa=NA_real_,
      phi=NA_real_,
      AIC = NA_real_,
      AICc = NA_real_,
      RSE=NA_real_,
      diagnostics=cond$message)
  }
  )  
} 

utils::globalVariables("x") #non-standard evaluation in curve