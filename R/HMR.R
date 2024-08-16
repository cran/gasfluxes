#' @title HMR fit 
#'  
#' @description
#'   Fit the HMR model using the Golub-Pereyra algorithm for partially linear least-squares models.
#'   
#' @param t time values (usually in hours)
#' @param C concentration values
#' @param A area covered by the chamber
#' @param V effective volume of the chamber
#' @param serie id of the flux measurement
#' @param k starting value for nls function
#' @param verbose logical, TRUE prints message after each flux calculation
#' @param plot logical, mainly intended for use in \code{\link{gasfluxes}}
#' @param maxiter see \code{\link{nls.control}}
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
#'    \item{RSE}{residual standard error (sigma from summary.nls)}
#'    \item{diagnostics}{error or warning messages}
#'    
#'    
#' @details
#' The HMR model (Pedersen et al., 2010) is \eqn{C(t)=\phi+f_0 \frac{e^{-\kappa t}}{-\kappa \frac{V}{A}}}{C(t) = \phi + f0 exp(-\kappa t)/(-\kappa V/A)}.
#' To ensure the lower bound \eqn{\kappa > 0}, the substitution \eqn{\kappa = e^k}{\kappa = exp(k)} is used. The resulting reparameterized model is then 
#' fit using \code{\link{nls}} with \code{algorithm = "plinear"}. This is computationally more efficient than the manual implementation in the HMR package and results 
#' in almost identical flux values. Flux standard errors and p-values differ strongly from those reported by the HMR package <= version 0.3.1, 
#' but are equal to those reported by later versions.
#' 
#' The default starting value \eqn{k = log(\kappa)} assumes that time is in hours. If you use a different time unit, you should adjust it accordingly.
#' 
#' There have been demands to return the initial concentration as predicted by the model as this is useful for checking plausibility. However, 
#' this can be easily calculated from the parameters and the equation of the model by setting t = 0, i.e., \eqn{C_0=\phi-\frac{f_0}{\kappa \frac{V}{A}}}{C0 = \phi - f0/(\kappa V/A)}.
#' 
#' Note that \code{nls} is used internally and thus this function should not be used with artificial "zero-residual" data.
#' 
#'
#' @references
#' Pedersen, A.R., Petersen, S.O., Schelde, K., 2010. A comprehensive approach to soil-atmosphere trace-gas flux estimation with static chambers. European Journal of Soil Science 61(6), 888-902.
#'      
#' @examples
#' #a single fit
#' t <- c(0, 1/3, 2/3, 1)
#' C <- c(320, 341, 352, 359)
#' print(fit <- HMR.fit(t, C, 1, 0.3, "a"))
#' plot(C ~ t)
#' curve({fit$phi + fit$f0 * exp(-fit$kappa * x)/(-fit$kappa*0.3)}, 
#'       from = 0, to = 1, add = TRUE)
#' 
#' \dontrun{
#' #a dataset of 1329 chamber N2O flux measurements
#' data(fluxMeas)
#' fluxMeas[, n := length(time), by=serie]
#' print(fluxMeas)
#' fluxes <- fluxMeas[n > 3, HMR.fit(time, C, A, V, serie), by=serie]
#' print(fluxes)
#' plot(f0.se ~ f0, data = fluxes)
#' #one very large f0.se value (and several infinite ones not shown in the plot)
#' fluxes[is.finite(f0.se),][which.max(f0.se),]
#' plot(C~time, data=fluxMeas[serie=="ID940",])
#' print(tmp <- fluxes[is.finite(f0.se),][which.max(f0.se),])
#' curve({tmp[, phi] + tmp[, f0] * exp(-tmp[, kappa] * x)/
#'       (-tmp[, kappa]*fluxMeas[serie=="ID940", V[1]]/
#'       fluxMeas[serie=="ID940",A[1]])}, 
#'       from = 0, to = 1, add = TRUE)
#' plot(f0.se ~ f0, data = fluxes[f0.se < 1e4,], pch = 16)
#' boxplot(fluxes[f0.se < 1e4, sqrt(f0.se)])
#' }
#' 
#' @importFrom stats nls nls.control predict AIC
#' @importFrom graphics curve 
#' @export

HMR.fit <- function (t, C, A = 1, V, serie = "", k = log(1.5), verbose = TRUE, plot = FALSE, maxiter = 100, ...) {

  tryCatch({
    stopifnot(length(t) > 3)
    fit <- withWarnings(nls(C ~ cbind(1, exp(-exp(k)*t)/(-exp(k)*V/A)), 
               start= list(k=k), algorithm = "plinear",
               control=nls.control(maxiter=maxiter, minFactor=1e-10, scaleOffset = 1)))
    w <- if (is.null(fit$warnings)) "" else fit$warnings[[1]]$message
    fit <- fit$value
    fitsum <- withWarnings(summary(fit))
    if (w == "") w <- if (is.null(fitsum$warnings)) "" else fitsum$warnings[[1]]$message
    fitsum <- fitsum$value
    fitsumCoef <- fitsum$coef
    try({
      if (plot) {
        curve(predict(fit, newdata = data.frame(t = x)), 
              from = min(t), to = max(t), add = TRUE, col = "red")
      }}, silent = TRUE)
    
    res <- list(
      f0 = fitsumCoef[".lin2", "Estimate"], 
      f0.se = fitsumCoef[".lin2", "Std. Error"], 
      f0.p = fitsumCoef[".lin2", "Pr(>|t|)"], 
      kappa = exp(fitsumCoef["k", "Estimate"]),
      phi = fitsumCoef[".lin1", "Estimate"],
      AIC = AIC(fit),
      AICc = aicc(fit),
      RSE = fitsum$sigma,
      diagnostics = w)
    if (res$f0.p > (1 - .Machine$double.neg.eps)) stop("Dubious fit (extreme standard error)")
    if (verbose) message(serie, if (w == "") ": HMR fit successful" else ": HMR fit warning")
    res
  },
  error = function(cond) {
    if (verbose) message(serie, ": HMR fit failed")
    list(
      f0 = NA_real_, 
      f0.se = NA_real_, 
      f0.p = NA_real_, 
      kappa = NA_real_,
      phi = NA_real_,
      AIC = NA_real_,
      AICc = NA_real_,
      RSE = NA_real_,
      diagnostics=cond$message)
  }
  )  
} 

utils::globalVariables("x") #non-standard evaluation in curve