### Changes in v0.4-2

#### Minor changes
1. In rare instances, HMR could produce a fit that is practically linear, which is probably not a numerically stable fit. Such fits could be selected by *selectfluxes* with the "kappa.max" algorithm, which is a problem because the standard error of the flux would become huge due to extremely high uncertainty of the curvature parameter. An attempt is made now at identifying such fits and selecting the robust linear fit instead.

### Changes in v0.4-1

#### Minor changes
1. Weights from robust linear regression are now separated by "|" insteaf of "," to better accomodate international users, [#12](https://bitbucket.org/ecoRoland/gasfluxes/issues/12/robustlinearweights-output-messes-up-comma)


### Changes in v0.4

#### NEW FEATURES
1. Added Pearson's correlation coefficient to output of linear fit.

2. Breaking change: Flux selection has now its own function *selectfluxes*. This makes implementing new algorithms easier. Please consult help("selectfluxes") for usage examples. 

3. Implemented Roman Hueppi's flux selection algorithm.

4. New Vignette included.

### Changes in v0.2-1

#### BUG FIXES

1. File names and titles of plots were not correctly using the flux ID, [#9](https://bitbucket.org/ecoRoland/gasfluxes/issues/9/plot-names-when-id-is-a-factor-variable).
