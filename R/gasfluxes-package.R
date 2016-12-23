#' Calculate greenhouse gas flux calculation from chamber measurements.
#'
#' Gasfluxes provides functions for fitting non-linear concentration - time models as well as convenience functions
#' for checking data and combining different calculation methods.
#'
#' The wrapper function for convenient flux calculation is
#' \code{\link{gasfluxes}}. Several concentration - time models are implemented
#' \itemize{
#' \item \code{\link{HMR.orig}}: The original implementation of HMR.
#' \item \code{\link{HMR.fit}}: A new implementation of HMR using partially linear least-squares. This is recommended over \code{HMR.orig}.
#' \item \code{\link{NDFE.fit}}: An implementation of the NDFE model using partially linear least-squares.
#' \item \code{\link{lin.fit}}: A simple linear model.
#' \item \code{\link{rlin.fit}}: A simple linear model fit using robust regression.
#' }
#'
#' @docType package
#' @name gasfluxes-package
#' @aliases gasfluxes-package
NULL