#check the input befor flux calculation
#dat must be data.table and contain columns specified as id, V, A, times, C


checkInput <- function (dat, .id, .V, .A, .times, .C) {
  
  stopifnot(is.data.table(dat))
  
  #check if specified columns exist
  stopifnot(all(c(.id, .V, .A, .times, .C) %in% names(dat)))
  
  #check if input is numeric
  testtype <- vapply(dat[, mget(c(.V, .A, .times, .C))], is.numeric, FUN.VALUE = TRUE)  
  if (!all(testtype)) stop(paste(c(.V, .A, .times, .C)[!testtype], collapse = ", "), 
                           ' not numeric.\nDecimal seperator must be ".". Check for text comments or typos. Missing values must be removed.',
                           call. = FALSE) 
  
  #check for NAs
  testNA <- vapply(dat[, mget(c(.V, .A, .times, .C))], function(x) any(!is.finite(x)), FUN.VALUE = TRUE)  
  if (any(testNA)) stop("Non-finite values (such as NA) in ", 
                         paste(c(.V, .A, .times, .C)[testNA], collapse = ", "), 
                         ". Remove them from input, e.g., using na.omit.",
                           call. = FALSE) 
 
  
  #non-unique values in V or A, time sorted, time values unique
  dat[, {
    if (length(unique(get(.V))) > 1) stop(.V, " not unique in flux ID ", paste(mget(.id), collapse = "_"), ".", call. = FALSE)
    if (length(unique(get(.A))) > 1) stop(.A, " not unique in flux ID ", paste(mget(.id), collapse = "_"), ".", call. = FALSE)
    if (any(diff(order(get(.times))) != 1L)) stop(.times, " not sorted in flux ID ", paste(mget(.id), collapse = "_"), ".", call. = FALSE)
    if (anyDuplicated(get(.times))) stop("Duplicated ", .times, " values in flux ID ", paste(mget(.id), collapse = "_"), ".", call. = FALSE)
    }, by = .id]
  
  #only positive input allowed
  testpos <- vapply(dat[, mget(c(.V, .A, .C))], function(x) any(x <= 0), FUN.VALUE = TRUE)  
  if (any(testpos)) stop("Not all values in ", 
                         paste(c(.V, .A, .C)[testpos], collapse = ", "), 
                         " are > 0.",
                         call. = FALSE) 
  if (any(dat[[.times]] < 0)) stop("Negative time values.", call. = FALSE)
  
  
  invisible(NULL)
} 