### Changes in v0.7
Bugfix: Whether an HMR fit was successful with difficult input data was platform dependent. Fits are now reported as failed with error message "Dubious fit (extreme standard error)" if the p value of the flux is > (1 - .Machine$double.neg.eps), which ensures that fits with extreme standard error of the flux estimate (in particular SE = Inf) are not reported anymore. This impacts the result of the recommended simulation to determine *f.detect* and thus the flux selection in *selectfluxes* on some platforms, especially on Windows systems. *fdetect* will be lower leading to successful HMR fits being more likely to be selected but there will also be less successful but dubious HMR fits. 

### Changes in v0.6-4
Fixed for CRAN note (regarding a documentation link).

### Changes in v0.6-2
Fixed an embarrassing maths error in the documentation.

### Changes in v0.6

#### Minor changes
Following a suggestion by Thomas Gremmen, diagnostic plots are now in a single PDF file. The (potentially numerous) PNG files created by previous versions caused issues on some network drives.

### Changes in v0.5
The old deprecated original HMR algorithm has been removed. Old deprecated flux selection algorithms have been removed. Removed dependency on package AICcmodavg.


### Changes in v0.4-4

#### Minor changes
1. Warning handlers did not work correctly.

2. There was a scoping issue with the plot parameter in *gasfluxes* if one of the ID columns was called "plot". 


### Changes in v0.4-3

#### Minor changes
1. *selectfluxes* has now gained a tolerance parameter to enable fine-tuning the check introduced in the previous version.

2. CRAN has asked for a package rebuild with the newest knitr version due to encoding issues in the vignette.


### Changes in v0.4-2

#### Minor changes
1. In rare instances, HMR could produce a fit that is practically linear, which is probably not a numerically stable fit. Such fits could be selected by *selectfluxes* with the "kappa.max" algorithm, which is a problem because the standard error of the flux would become huge due to extremely high uncertainty of the curvature parameter. An attempt is made now at identifying such fits and selecting the robust linear fit instead.

### Changes in v0.4-1

#### Minor changes
1. Weights from robust linear regression are now separated by "|" instead of "," to better accommodate international users.


### Changes in v0.4

#### NEW FEATURES
1. Added Pearson's correlation coefficient to output of linear fit.

2. Breaking change: Flux selection has now its own function *selectfluxes*. This makes implementing new algorithms easier. Please consult help("selectfluxes") for usage examples. 

3. Implemented Roman Hueppi's flux selection algorithm.

4. New Vignette included.

### Changes in v0.2-1

#### BUG FIXES

1. File names and titles of plots were not correctly using the flux ID.
