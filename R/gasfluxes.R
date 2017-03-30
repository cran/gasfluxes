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
#' @param select character; specify a ruleset for selection of the final flux value, use NULL for no selection.
#' @param k_HMR starting value for \code{\link{HMR.fit}}.
#' @param k_NDFE starting value for \code{\link{NDFE.fit}}.
#' @param verbose logical; print progress messages?
#' @param plot create plots if TRUE (the default). The IDs are used as file names and should thus not include characters that are not allowed in file names, such as \code{/} or \code{=}. A directory "pics"  is created in the working directory if it doesn't exist. The plots are only intended to facilitate quick checking, not for publication quality graphs.
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
#' "original HMR":\tab \code{\link{HMR.orig}}\cr
#' "NDFE":\tab \code{\link{NDFE.fit}}\cr
#' }
#' Specifying other methods results in an error. 
#' 
#' The default starting values for "HMR" and "NDFE", \eqn{k = log(\kappa)} and \eqn{k = log(\tau)}, resp., assume that time is in hours. If you use a different time unit, you should adjust them accordingly.
#' Note that \code{nls} is used internally by these functions and thus they should not be used with artificial "zero-residual" data.
#' 
#' Available selection algorithms currently are
#' \describe{
#' \item{"RF2011"}{The algorithm used, e.g., in Leiber-Sauheitl 2014 (doi:10.5194/bg-11-749-2014).
#' This overwrites the methods parameter. The factor guarding against degenerate HMR fits
#' can be set via the ellipsis as \code{gfactor}. Default is \code{gfactor = 4}.}
#' \item{"RF2011new"}{The same rules as "RF2011", but using the improved fitting function for HMR,
#' which results in larger SE and p-values. Thus, it is less likely to select the HMR result.
#' This selection algorithm might need further testing and validation.}
#' }
#' 
#' Other selection algorithms could be implemented, but selection can always be done as a postprocessing step. E.g., if many data points are available for each flux measurement it is probably most sensible to use AICc.
#' 
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
#'      
#' @examples
#' \dontrun{
#' #compare result of original HMR with plinear HMR
#' data(fluxMeas)
#' res <- gasfluxes(fluxMeas[1:400,], 
#'                  .id = "serie", .V = "V", .A = "A",
#'                  .times = "time", .C = "C",
#'                  methods = c("HMR", "original HMR"), verbose = TRUE)
#'                  
#' #number of successful fits
#' res[, sum(!is.na(HMR.kappa))]
#' res[, sum(!is.na(original.HMR.kappa))]
#' 
#' #Do the results differ?
#' plot(res[["HMR.f0"]], res[["original.HMR.f0"]])
#' abline(0, 1)
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
#' @importFrom grDevices png dev.off
#' @importFrom graphics layout plot.new legend par
#' @importFrom stats setNames
#' @export

gasfluxes <- function (dat, .id = "ID", .V = "V", .A = "A", .times = "time", .C = "C", 
                       methods = c("linear", "robust linear", "HMR", "NDFE"), 
                       select = NULL,
                       k_HMR = log(1.5), k_NDFE = log(0.01), verbose = TRUE,
                       plot = TRUE, ...) {
  
  stopifnot(is.data.frame(dat))
  stopifnot(all(lengths(list(.V, .A, .times, .C)) == 1L))
  isDT <- is.data.table(dat)
  
  setDT(dat)
  
  
  checkInput(dat, .id, .V, .A, .times, .C)
  
  if (!is.null(select) && select == "RF2011") methods <- c("linear", "robust linear", "original HMR")
  if (!is.null(select) && select == "RF2011new") methods <- c("linear", "robust linear", "HMR")
  
  #create dir for plots
  if (plot) {
    #create folder /pics in working directory if it doesn't exist
    mainDir <- getwd()
    subDir <- "pics"
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)  
  }
  
  
  funs <- list("linear" = lin.fit,
               "robust linear" = function(...) {
                 res <- rlin.fit(...)
                 res[["weights"]] <- paste(signif(res[["weights"]], 2), collapse=",")
                 res
               },
               "HMR" = function(...) HMR.fit(..., k = k_HMR),
               "original HMR" = HMR.orig,
               "NDFE" = function(...) NDFE.fit(..., k = k_HMR))
  stopifnot(all(methods %in% names(funs)))
  funs <- funs[methods]
  
  callFun <- function(fun, .t, .C, .A, .V, .ID, verbose, plot) {
    
    force(fun(t = .t, C = .C, A = .A, V = .V, serie = .ID, verbose = verbose, plot = plot))
  }
  
  
  res <- dat[, {
    .ID <- do.call(paste, c(mget(.id), sep = "_"))[1]
    #plot for each ID
    if (plot) {
      png(filename = paste0(paste(file.path(mainDir, subDir), .ID, sep = "/"), ".png"),
          width = 600, height = 480, units = "px", pointsize = 12)
      
      
      layout(t(c(1, 2)), widths = c(3, 1))
      plot(get(.times), get(.C), pch = 16,
           xlab = "time", ylab = "concentration", main = .ID)
    }
    #fits for each ID
    tmp <- unlist(
      setNames(
        lapply(funs, callFun, .t=get(.times), .C=get(.C), .A=get(.A)[1], .V=get(.V)[1], .ID=.ID, verbose=verbose, plot = plot), 
        make.names(names(funs))),
      recursive=F) 
    if (plot) {
      cols <- c("linear" = "black", "robust linear" = "green",
                "HMR" = "red", "original HMR" = "yellow",
                "NDFE" = "blue")
      ltys <- c("linear" = 1, "robust linear" = 2,
                "HMR" = 1, "original HMR" = 2,
                "NDFE" = 1)
      par(mar = rep(0.1, 4))
      plot.new()
      
      legend("center", legend = methods, col = cols[methods], lty = ltys[methods])
      dev.off()
    }
    tmp
  },   
  by=.id]
  
  #RF2011
  if (!is.null(select) && select == "RF2011") {
    if (!("gfactor" %in% names(list(...)))) gfactor <- 4 else gfactor <- list(...)[["gfactor"]]
    res[, c("flux", "flux.se", "flux.p", "method") := list(robust.linear.f0,
                                                           robust.linear.f0.se,
                                                           robust.linear.f0.p,
                                                           "robust linear")]
    res[!is.finite(flux), c("flux", "flux.se", "flux.p", "method") := list(linear.f0,
                                                                           linear.f0.se,
                                                                           linear.f0.p,
                                                                           "linear")]
    res[method == "robust linear" & 
          is.finite(original.HMR.f0) &
          is.finite(original.HMR.kappa) &
          (original.HMR.AIC < linear.AIC) &
          (original.HMR.f0.p < robust.linear.f0.p) &
          (abs(original.HMR.f0) <= abs(robust.linear.f0) * gfactor), 
        c("flux", "flux.se", "flux.p", "method") := list(original.HMR.f0,
                                                         original.HMR.f0.se,
                                                         original.HMR.f0.p,
                                                         "original HMR")]
    res[!is.finite(flux), method := "error"]
  }
  
  #RF2011new
  if (!is.null(select) && select == "RF2011new") {
    if (!("gfactor" %in% names(list(...)))) gfactor <- 4 else gfactor <- list(...)[["gfactor"]]
    res[, c("flux", "flux.se", "flux.p", "method") := list(robust.linear.f0,
                                                           robust.linear.f0.se,
                                                           robust.linear.f0.p,
                                                           "robust linear")]
    res[!is.finite(flux), c("flux", "flux.se", "flux.p", "method") := list(linear.f0,
                                                                           linear.f0.se,
                                                                           linear.f0.p,
                                                                           "linear")]
    res[method == "robust linear" & 
          is.finite(HMR.f0) &
          (HMR.AIC < linear.AIC) &
          (HMR.f0.p < robust.linear.f0.p) &
          (abs(HMR.f0) <= abs(robust.linear.f0) * gfactor), 
        c("flux", "flux.se", "flux.p", "method") := list(HMR.f0,
                                                         HMR.f0.se,
                                                         HMR.f0.p,
                                                         "HMR")]
    res[!is.finite(flux), method := "error"]
  }
  
  
  if (!isDT) {
    setDF(dat)
  }
  res[]
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
                         'HMR.f0.p','HMR.f0.se')) #non-standard evaluation in data.table