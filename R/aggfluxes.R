#' @title Accumulation of fluxes
#'
#' @description
#'   Aggregate a time series of fluxes to a cummulative flux value.
#'   
#' @param fluxes flux values
#' @param datetimes datetime values (POSIXct or POSIXlt)
#' @param timeunit the unit of time (denominator of the flux unit), supported are the explicit units supported by \code{difftime} 
#'  
#' @return
#'  A one-row data.frame with columns
#'    \item{flux}{the cumulative flux}
#'    \item{from}{the start of the cumulation period}
#'    \item{to}{the end of the cumulation period}
#'  The return value being a data.frame is useful, when the function is used for "split-apply-combine" type operations to calculate groupwise cumulated values, 
#'  e.g., using package data.table.
#'    
#' @details
#' The function uses linear interpolation. The unit of the cumulative flux is [fluxes] * timeunit. 
#' NA values are removed and values sorted according to time order. If less then two non-NA value pairs are provided, 
#' NA is returned for the cumlative flux.
#'      
#' @examples
#' #Some random example data
#' datetimes <- Sys.time() + (1:20)/2*24*3600
#' set.seed(42)
#' fluxes <- rlnorm(20, 5)
#' agg.fluxes(fluxes, datetimes)
#' 
#' @importFrom stats filter
#' @export

agg.fluxes <- function(fluxes, datetimes, timeunit = "hours"){
  stopifnot(timeunit %in% c("secs", "mins", "hours",
                            "days", "weeks"))
  stopifnot(is.numeric(fluxes))
  stopifnot("POSIXt" %in% class(datetimes))
  keep <- !is.na(fluxes) & !is.na(datetimes)
  fluxes <- fluxes[keep]
  datetimes <- datetimes[keep]
  fluxes <- fluxes[order(datetimes)]
  datetimes <- datetimes[order(datetimes)]
  
  if (length(datetimes)>1){
    #mean of subsequent values
    .mean <- filter(fluxes,c(1,1))/2
    .mean <- .mean[-length(.mean)]
    #time difference between subsequent measurements
    dtime <- diff(datetimes)
    units(dtime) <- timeunit
    dtime <- as.numeric(dtime)
    #aggregated flux
    res <- sum(.mean*dtime)
    
    #return values
    data.frame(flux=res,from=min(datetimes),to=max(datetimes))
  } else{
    data.frame(flux=NA_real_,from=min(datetimes),to=max(datetimes)) 
  }
}