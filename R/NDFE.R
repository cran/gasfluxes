#' @title erfc
#' 
#' @description
#'   This is the complementary error function.
#' @param x a numeric vector
#'  
#' @return
#'  A numeric vector, i.e., the erfc values.
#'
#' @importFrom stats pnorm
#' @export


#the complementary error function
erfc <- function(x) 2*pnorm(-sqrt(2)*x)


#' @title NDFE fit 
#' 
#' @description
#'   Fit the the non-steady-state diffusive flux extimator model using the Golub-Pereyra algorithm for partially linear least-squares models.
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
#'    \item{C0, tau}{other parameters of the NDFE model}
#'    \item{AIC}{Akaike information criterion}
#'    \item{AICc}{Akaike information criterion with small sample correction}
#'    \item{RSE}{residual standard error (sigma from summary.nls)}
#'    \item{diagnostics}{error or warning messages}
#'    
#' @details
#' The NDFE model (Livingston et al., 2006) is \eqn{C(t)=C_0+f_0\tau \frac{A}{V}\left [\frac{2}{\sqrt{\textup{pi}}}\sqrt{t/\tau}+e^{t/\tau}\textup{erfc}(\sqrt{t/\tau})-1\right ]}{C(t)=C0+f0\tau A/V[2/\sqrt(pi)\sqrt{t/\tau}+exp(t/\tau)erfc(\sqrt(t/\tau)-1]}.
#' To ensure the lower bound \eqn{\tau > 0}, the substituion \eqn{\tau = e^k}{\tau = exp(k)} is used. The resulting reparameterized model is then 
#' fit using \code{\link{nls}} with \code{algorithm = "plinear"}. 
#' 
#' Note that according to the reference the model is not valid for negative fluxes. Warning: This function does not check if fluxes are positive. It's left to the user to handle negative fluxes.
#'
#' The default starting value \eqn{k = log(\tau)} assumes that time is in hours. If you use a different time unit, you should adjust it accordingly.
#' 
#' Note that \code{nls} is used internally and thus this function should not be used with artificial "zero-residual" data.
#'
#' @references
#' Livingston, G.P., Hutchinson, G.L., Spartalian, K., 2006. Trace gas emission in chambers: A non-steady-state diffusion model. Soil Sci. Soc. Am. J. 70(5), 1459-1469.
#'      
#' @examples
#' #a single fit
#' t <- c(0, 1/3, 2/3, 1)
#' C <- c(320, 340, 355, 362)
#' print(fit <- NDFE.fit(t, C, 1, 0.3, "a"))
#' plot(C ~ t)
#' curve({fit$C0+fit$f0*fit$tau*1/0.3*(2/sqrt(pi)*sqrt(x/fit$tau)+
#'       exp(x/fit$tau)*erfc(sqrt(x/fit$tau))-1)}, 
#'       from = 0, to = 1, add = TRUE)
#' #note that the flux estimate is very uncertain because 
#' #there are no data points in the region of high curvature
#' 
#' @export



NDFE.fit <- function (t, C, A = 1, V, serie = "", k = log(0.01), verbose = TRUE, plot = FALSE, maxiter = 100, ...) {
  
  tryCatch({
    stopifnot(length(t) > 3)
    fit <- withWarnings(nls(C ~ cbind(1, exp(k)*(A/V)*(2/sqrt(pi)*sqrt(t/exp(k))+exp(t/exp(k))*erfc(sqrt(t/exp(k)))-1)), 
               start= list(k=k), algorithm = "plinear",
               control=nls.control(maxiter=maxiter, minFactor=1e-10)))
    w <- if (is.null(fit$warnings)) "" else fit$warnings[[1]]$message
    fit <- fit$value
    fitsum <- withWarnings(summary(fit))
    if (w == "") w <- if (is.null(fitsum$warnings)) "" else fitsum$warnings[[1]]$message
    fitsum <- fitsum$value
    fitsumCoef <- fitsum$coef
    
    try({
      if (plot) {
        curve(predict(fit, newdata = data.frame(t = x)), 
              from = min(t), to = max(t), add = TRUE, col = "blue")
    }}, silent = TRUE)
        
    res <- list(
      f0 = fitsumCoef[".lin2", "Estimate"], 
      f0.se = fitsumCoef[".lin2", "Std. Error"], 
      f0.p = fitsumCoef[".lin2", "Pr(>|t|)"], 
      tau = exp(fitsumCoef["k", "Estimate"]),
      C0 = fitsumCoef[".lin1", "Estimate"],
      AIC = AIC(fit),
      AICc = aicc(fit),
      RSE = fitsum$sigma,
      diagnostics = w)
    if (verbose) message(serie, if (w == "") ": NDFE fit successful" else ": NDFE fit warning")
    res
  },
  error = function(cond) {
    if (verbose) message(serie, ": NDFE fit failed")
    list(
      f0 = NA_real_, 
      f0.se = NA_real_, 
      f0.p = NA_real_, 
      tau = NA_real_,
      C0 = NA_real_,
      AIC = NA_real_,
      AICc = NA_real_,
      RSE = NA_real_,
      diagnostics=cond$message)
  }
  )  
} 

utils::globalVariables("x") #non-standard evaluation in curve
