#' @title Flux calculation
#' 
#' @description
#'  A wrapper function for convenient flux calculation.
#' @param dat a data.frame or data.table with data from flux measurements.
#' @param .id character vector specifying the columns to be used as ID, multiple ID columns are possible.
#' @param .V character specifying the column containing chamber volume values.
#' @param .A character specifying the column containing chamber area values.
#' @param .times character specifying the column containing chamber closing time values.
#' @param .C character specifying the column containing concentration values.
#' @param methods character; which methods to use for flux estimation. See details for available methods.
#' @param k_HMR starting value for \code{\link{HMR.fit}}.
#' @param k_NDFE starting value for \code{\link{NDFE.fit}}.
#' @param verbose logical; print progress messages?
#' @param plot create a PDF with plots in the working directory if \code{TRUE} (the default). The IDs are used as plot names. The plots are only intended to facilitate quick checking, not for publication quality graphs.
#' @param select deprecated; please use function \code{\link{selectfluxes}}.
#' @param maxiter see \code{\link{nls.control}}
#' @param \dots further parameters
#'  
#' @return
#'  A data.table with the results of the flux calculation. See the documentation of the fitting functions for details. 
#'  If a selection algorithm has been specified, the last columns are the selected flux estimate, the corresponding standard error and p-value and the method with which the selected flux was estimated.
#'    
#' @details
#' Available methods are
#' \tabular{ll}{
#' "linear":\tab \code{\link{lin.fit}}\cr
#' "robust linear":\tab \code{\link{rlin.fit}}\cr
#' "HMR":\tab \code{\link{HMR.fit}}\cr
#' "NDFE":\tab \code{\link{NDFE.fit}}\cr
#' }
#' Specifying other methods results in an error. 
#' 
#' The default starting values for "HMR" and "NDFE", \eqn{k = log(\kappa)} and \eqn{k = log(\tau)}, resp., assume that time is in hours. If you use a different time unit, you should adjust them accordingly.
#' Note that \code{nls} is used internally by these functions and thus they should not be used with artificial "zero-residual" data.
#' 
#' The input data.frame or data.table should be in the following format:
#' 
#' \preformatted{
#'    serie        V A      time         C
#'1:    ID1 0.522625 1 0.0000000 0.3317823
#'2:    ID1 0.522625 1 0.3333333 0.3304053
#'3:    ID1 0.522625 1 0.6666667 0.3394311
#'4:    ID1 0.522625 1 1.0000000 0.4469102
#'5:    ID2 0.523625 1 0.0000000 0.4572708
#'}
#'
#'However, more than one ID column are possible. E.g., the first ID column could be the plot and a second ID column could be the date. Keep in mind that the combination of IDs must be a unique identifier for each flux measurement.
#'
#'Units of the output depend on input units. 
#'It's recommended to use [V] = m^3, [A] = m^2, [time] = h, [C] = [mass or mol]/m^3, which results in [f0] = [mass or mol]/m^2/h.
#'Since all algorithms use V/A, A can be input as 1 and V as the chamber height.
#' 
#' @seealso 
#' \code{\link{selectfluxes}} for flux selection
#'           
#' @examples
#' \dontrun{
#' #compare result of original HMR with plinear HMR
#' data(fluxMeas)
#' res <- gasfluxes(fluxMeas[1:400,], 
#'                  .id = "serie", .V = "V", .A = "A",
#'                  .times = "time", .C = "C",
#'                  methods = c("HMR"), verbose = TRUE)
#'                  
#' #number of successful fits
#' res[, sum(!is.na(HMR.kappa))]
#' 
#' 
#' 
#' res <- gasfluxes(fluxMeas, 
#'                  .id = "serie", .V = "V", .A = "A",
#'                  .times = "time", .C = "C",
#'                  methods = "HMR", verbose = TRUE)
#' # Error: time not sorted in flux ID ID556. 
#' # Investigate the problem:
#' fluxMeas[serie %in% c("ID555", "ID556", "ID557")]
#' #    serie        V A      time         C
#' # 1: ID555 0.551625 1 0.0000000 0.3884388
#' # 2: ID555 0.551625 1 0.3333333 0.4125270
#' # 3: ID555 0.551625 1 0.6666667 0.3714207
#' # 4: ID555 0.551625 1 1.0000000 0.3735092
#' # 5: ID556 0.524250 1 0.0000000 0.3638239
#' # 6: ID556 0.524250 1 0.3333333 0.3520481
#' # 7: ID556 0.524250 1 0.6666667 0.3551644
#' # 8: ID557 0.528375 1 0.0500000 0.3954601
#' # 9: ID556 0.524250 1 0.0000000 0.3839834
#' #10: ID557 0.528375 1 0.3333333 0.3967269
#' #11: ID557 0.528375 1 0.6666667 0.3764967
#' #12: ID557 0.528375 1 1.0000000 0.3973055
#' 
#' # some mixup of IDs and times
#' # usually an input or Excel error during data preparation
#' # investigate and fix
#' } 
#' 
#' @import data.table
#' @importFrom grDevices pdf dev.off dev.cur
#' @importFrom graphics layout plot.new legend par
#' @importFrom stats setNames
#' @export

gasfluxes <- function (dat, .id = "ID", .V = "V", .A = "A", .times = "time", .C = "C", 
                       methods = c("linear", "robust linear", "HMR", "NDFE"),
                       k_HMR = log(1.5), k_NDFE = log(0.01), verbose = TRUE,
                       plot = TRUE, select, maxiter = 100, ...) {
  
  #close pdf device
  on.exit({
    if (plot && exists("curdev") && curdev == dev.cur()) invisible(dev.off(which = curdev))
  })
  
  stopifnot(is.data.frame(dat))
  stopifnot(all(lengths(list(.V, .A, .times, .C)) == 1L))
  if (!missing(select)) stop("select has been deprecated, please refer to the documentation.")
  
  setDT(dat)
  
  
  checkInput(dat, .id, .V, .A, .times, .C)
  
  #create and open PDF file for plots
  if (plot) {
    pdf(file = sprintf("gasfluxes_plots_%s.pdf", format(Sys.time(), "%Y%m%d_%H%M%S")),
        paper = "default", width = 0, height = 0,
        )
    
    curdev <- dev.cur()
    
    layout(matrix(1:6, ncol = 2, byrow = TRUE))
      cols <- c("linear" = "black", "robust linear" = "green",
                "HMR" = "red", "original HMR" = "yellow",
                "NDFE" = "blue")
      ltys <- c("linear" = 1, "robust linear" = 2,
                "HMR" = 1, "original HMR" = 2,
                "NDFE" = 1)
      
      plot.new()
      legend("center", legend = methods, col = cols[methods], lty = ltys[methods],
             cex = 1.5)

  }
  
  
  funs <- list("linear" = lin.fit,
               "robust linear" = function(...) {
                 res <- rlin.fit(...)
                 res[["weights"]] <- paste(signif(res[["weights"]], 2), collapse="|")
                 res
               },
               "HMR" = function(...) HMR.fit(..., k = k_HMR, maxiter = maxiter),
               "NDFE" = function(...) NDFE.fit(..., k = k_HMR, maxiter = maxiter))
  stopifnot(all(methods %in% names(funs)))
  funs <- funs[methods]
  
  callFun <- function(fun, .t, .C, .A, .V, .ID, verbose, plot) {
    
    force(fun(t = .t, C = .C, A = .A, V = .V, serie = .ID, verbose = verbose, plot = plot))
  }
  
  
  res <- dat[, {
    .ID <- do.call(paste, c(mget(.id), sep = "_"))[1]
    #plot for each ID
    if (..plot) {
      plot(get(.times), get(.C), pch = 16,
           xlab = "time", ylab = "concentration", main = .ID)
    }
    #fits for each ID
    tmp <- unlist(
      setNames(
        lapply(funs, callFun, .t=get(.times), .C=get(.C), .A=get(.A)[1], .V=get(.V)[1], .ID=.ID, verbose=verbose, plot = plot), 
        make.names(names(funs))),
      recursive=F) 

    tmp
  },   
  by=.id]
  
  res[]
} 

#' @title Select a flux estimate 
#' 
#' @description
#'  Selects the appropriate flux estimate from linear, robust linear and non-linear calculated fluxes.
#' @param dat a data.table as returned by \code{\link{gasfluxes}}. The function modifies it by reference.
#' @param select character; specify a ruleset for selection of the final flux value, see details.
#' @param f.detect detection limit for HMR method. This can be determined by a simple simulation (see examples) or for four data points the approximation in Parkin et al. (2012) can be used. 
#' @param t.meas a vector or single value giving the measurement time factor that relates to kappa.max. It is suggested to use the time difference between the first and last sample taken from the closed chamber. The unit should be consistent with the units of \code{f.detect} and kappa (e.g., h if kappa is in 1/h).
#' @param tol the relative tolerance \code{abs((linear.f0 - HMR.f0)/HMR.f0)} below which the linear flux estimate and the HMR flux estimate are considered equal in the "kappa.max" algorithm. This is to protect against HMR fits that equal the linear fit and have extremely high standard errors. Defaults to \code{tol = 5e-5}. 
#' @param \dots further parameters
#' 
#'  
#' @return
#'  A data.table with the with following columns added to the function input: selected flux estimate, the corresponding standard error and p-value and the method with which the selected flux was estimated. 
#'  For the  "kappa.max" method the "kappa.max" values are included. These columns are also added to the input data.table by reference. 
#'    
#' @details
#' Available selection algorithms currently are
#' \describe{
#' \item{"kappa.max"}{The selection algorithm restricts the use of HMR by imposing a maximal value for kappa "kappa.max", 
#' depending on the quotient of the linear flux estimate and the minimal detectable flux (f.detect), as well as 
#' the chamber closure time (t.meas). kappa.max = f.lin/f.detect/t.meas. This is currently the recommended algorithm. 
#' Note that the algorithm was developed for predominantly positive fluxes (such as N2O fluxes). If data with considerable gas uptake is analyzed, the algorithm needs to be modified, which currently means the user needs to implement it themselves.}
#' }
#' 
#' Other selection algorithms could be implemented, but selection can always be done as a postprocessing step. E.g., if many data points are available for each flux measurement it is probably most sensible to use AICc.
#' @references
#' Parkin, T.B., Venterea, R.T., Hargreaves, S.K., 2012. Calculating the Detection Limits of Chamber-based Soil Greenhouse Gas Flux Measurements. Journal of Environmental Quality 41, 705-715.
#' 
#' Hueppi, R., Felber, R., Krauss, M., Six, J., Leifeld, J., Fuss, R., 2018. Restricting the nonlinearity parameter in soil greenhouse gas flux calculation for more reliable flux estimates. PLOS ONE 13(7): e0200876. https://doi.org/10.1371/journal.pone.0200876
#'  
#' @examples
#' \dontrun{
#' res <- gasfluxes(fluxMeas[1:499], 
#'                  .id = "serie", .V = "V", .A = "A",
#'                  .times = "time", .C = "C",
#'                  methods = c("linear", "robust linear", "HMR"), verbose = FALSE, plot = FALSE)
#' 
#' ### estimate f.detect by simulation ###
#' #ambient concentration:
#' C0 <- 320/1000 * 28 * 273.15 / 22.4 / (273.15 + 15) #mg N / m^3
#' #uncertainty of GC measurement:
#' sdGC <- 5/1000 * 28 * 273.15 / 22.4 / (273.15 + 15) #mg N / m^3 
#' #create simulated concentrations corresponding to 1 hour flux measurements with zero fluxes:
#' set.seed(42)
#' sim <- data.frame(t = seq(0, 1, length.out = 4), C = rnorm(4e3, mean = C0, sd = sdGC), 
#'                   id = rep(1:1e3, each = 4), A = 1, V = 0.52)
#' #fit HMR model:                  
#' simflux <- gasfluxes(sim, .id = "id", .times = "t", methods = c("HMR", "linear"), plot = FALSE) 
#' simflux[, f0 := HMR.f0]
#' simflux[is.na(f0), f0 := linear.f0]
#' #dection limit as 97.5 % quantile (95 % confidence):
#' f.detect <- simflux[, quantile(f0, 0.975)] #0.03 mg N / m^2 / h
#' 
#' # example using the kappa.max (ref. Hueppi et al., 2018) with a single t.meas value
#' t.meas <- max(fluxMeas$time[1:499]) #1
#' selectfluxes(res, "kappa.max", f.detect = f.detect, t.meas = t.meas)
#' res[method == "HMR", .N] # 11        
#' 
#' # example using the kappa.max with a vector for t.meas
#' t.meas <- fluxMeas[1:499][, max(time), by = serie][["V1"]]
#' selectfluxes(res, "kappa.max", f.detect = f.detect, t.meas = t.meas)
#' res[method == "HMR", .N] # 10                  
#' }
#' 
#' @export

selectfluxes <- function(dat, select, f.detect = NULL, t.meas = NULL, tol = 5e-5, ...) {
  stopifnot(is.data.table(dat))
  if (!(select %in% c("kappa.max"))) stop('Please specify an implemented algorithm, see help("select.fluxes").')
  
  if (select == "kappa.max") {
    if (!all(c("robust.linear.f0",
               "robust.linear.f0.se",
               "robust.linear.f0.p",
               "linear.f0",
               "linear.f0.se",
               "linear.f0.p",
               "linear.r",
               "HMR.f0",
               "HMR.f0.se",
               "HMR.f0.p",
               "HMR.kappa") %in% names(dat))) stop("Please use return value of gasfluxes function containing results of linear, robust linear and HMR regression.")
    
    if (is.null(f.detect))  stop("Please specify the minimal detectable flux (f.detect) of your system in the selectflux function (i.e. by using the f.detect simulation function described in the help file).")
    if (is.null(t.meas))  stop("Please specify the chamber measurement time of the system that was used to calculte the minimal detectable flux (f.detect) i.e. max(t).")
    if (any(t.meas <= 0))   stop("Measurement time factor t.meas needs to be always > 0.")
    if (length(t.meas) != 1 && length(t.meas) != nrow(dat)) stop("Measurement time factor t.meas needs to have the length 1 (same value for all measurements) or the length of the number of fluxes to be calculated (individual value for each measurement).")
    
    kappa.max <- function(f.lin, t.meas, f.detect) f.lin/f.detect/t.meas
    
    dat[, kappa.max := kappa.max(linear.f0, t.meas, f.detect)]
    
    dat[, c("flux", "flux.se", "flux.p", "method") := list(robust.linear.f0,
                                                           robust.linear.f0.se,
                                                           robust.linear.f0.p,
                                                           "robust linear")]
    dat[!is.finite(flux), c("flux", "flux.se", "flux.p", "method") := list(linear.f0,
                                                                           linear.f0.se,
                                                                           linear.f0.p,
                                                                           "linear")]
    dat[method == "robust linear" & 
          is.finite(HMR.f0) &
          (HMR.kappa < kappa.max) & 
          !(abs(HMR.f0.se/HMR.f0) > 1e10 & abs((linear.f0 - HMR.f0)/HMR.f0) < tol), #protect against linear HMR fit 
        c("flux", "flux.se", "flux.p", "method") := list(HMR.f0,
                                                         HMR.f0.se,
                                                         HMR.f0.p,
                                                         "HMR")]
    dat[!is.finite(flux), method := "error"]
  }
  dat[]
}

utils::globalVariables(c('robust.linear.f0',
                         'robust.linear.f0.se',
                         'robust.linear.f0.p',
                         'flux', 'linear.f0',
                         'linear.f0.se', 'linear.f0.p',
                         'method', 'original.HMR.f0',
                         'original.HMR.AIC', 'linear.AIC',
                         'original.HMR.f0.p', 'original.HMR.f0.se',
                         'original.HMR.kappa',
                         'HMR.f0', 'HMR.AIC',
                         'HMR.f0.p','HMR.f0.se', "HMR.kappa",
                         'i', 'robust.linear.weights', 'V1',
                         'linear.r',
                         '..plot')) #non-standard evaluation in data.table
