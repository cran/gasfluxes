% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roblinear.R
\name{rlin.fit}
\alias{rlin.fit}
\title{Robust linear concentration - time model}
\usage{
rlin.fit(t, C, A = 1, V, serie = "", verbose = TRUE, plot = FALSE, ...)
}
\arguments{
\item{t}{time values (usually in hours)}

\item{C}{concentration values}

\item{A}{area covered by the chamber}

\item{V}{effective volume of the chamber}

\item{serie}{id of the flux measurement}

\item{verbose}{logical, TRUE prints message after each flux calculation}

\item{plot}{logical, mainly intended for use in \code{\link{gasfluxes}}}

\item{\dots}{further parameters, currently none}
}
\value{
A list of
   \item{f0}{flux estimate}
   \item{f0.se}{standard error of flux estimate}
   \item{f0.p}{p-value of flux estimate}
   \item{C0}{estimated concentration at t = 0 (intercept)}  
   \item{weights}{robustness weights}
   \item{diagnostics}{error or warning messages}
}
\description{
Fit a linear model to concentration - time data using robust methods.
}
\details{
This is basically a wrapper of \code{\link[MASS]{rlm}} using the Huber M estimator. This function never weights the first or last time point with zero with very few data points. However, there might exist "better" robust regression methods for flux estimation.
}
\examples{
#a single fit
t <- c(0, 1/3, 2/3, 1)
C <- c(320, 330, 315, 351)
print(fit <- rlin.fit(t, C, 1, 0.3, "a"))
plot(C ~ t)
curve({fit$f0/0.3 * x + fit$C0}, from = 0, to = 1, add = TRUE)

}
