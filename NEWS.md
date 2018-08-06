### Changes in v0.4

#### NEW FEATURES
1. Added Pearson's correlation coefficient to output of linear fit.

2. Breaking change: Flux selection has now its own function *selectfluxes*. This makes implementing new algorithms easier. Please consult help("selectfluxes") for usage examples. 

3. Implemented Roman Hueppi's flux selection algorithm.

4. New Vignette included.

### Changes in v0.2-1

#### BUG FIXES

1. File names and titles of plots were not correctly using the flux ID, [#9](https://bitbucket.org/ecoRoland/gasfluxes/issues/9/plot-names-when-id-is-a-factor-variable).
