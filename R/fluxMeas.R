#' Data from chamber N2O flux measurements.
#'
#' A dataset containing data from 1329 chamber N2O flux measurements.
#'
#' @format A data.table with 5300 rows and 5 variables:
#' \itemize{
#'   \item serie: ID of flux measurement
#'   \item V: Volume (normalized by area, i.e., the height in m)
#'   \item A: Area (always 1)
#'   \item time: closing time in h
#'   \item C: N2O concentration in mg N / m^3
#' }
#' @source own data (anonymized by not including site and treatment information)
#' @name fluxMeas
NULL