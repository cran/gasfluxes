library(gasfluxes)
context("fitting checks")

test_that("linear regression gives expected result", {
  t <- c(0, 1/3, 2/3, 1)
  C <- c(320, 341, 352, 359)
  fit <- lin.fit(t, C, 1, 0.3, "a", verbose = FALSE)
  expect_equal(fit, structure(list(f0 = 11.52, f0.se = 2.02849698052524, f0.p = 0.0296345042190038, 
                                   C0 = 323.8, AIC = 27.5179162394956, AICc = Inf, RSE = 5.03984126734168, 
                                   r = 0.970365495780996, diagnostics = ""), 
                              .Names = c("f0", "f0.se", "f0.p", "C0", "AIC", "AICc", "RSE", "r","diagnostics")), tolerance = 1e-5)
  
  C <- 320 + 0:3 * 10
  expect_equal(lin.fit(t, C, 1, 0.3, "a", verbose = FALSE)$diagnostics, "essentially perfect fit: summary may be unreliable")
})

test_that("robust linear regression gives expected result", {
  skip_on_cran() #no idea why this test doesn't pass on winbuilder (Win Server 2008)
  
  t <- c(0, 1/3, 2/3, 1)
  C <- c(320, 341, 352, 359)
  fit <- rlin.fit(t, C, 1, 0.3, "a", verbose = FALSE)
  expect_equal(fit, structure(list(f0 = 11.52, f0.se = 2.02849698052524, f0.p = 0.0296345042190038, 
                                   C0 = 323.8, weights = c(1, 1, 1, 1), diagnostics = ""), 
                              .Names = c("f0", "f0.se", "f0.p", "C0", "weights", "diagnostics")), tolerance = 1e-5)
  C[2] <- 400
  fit <- rlin.fit(t, C, 1, 0.3, "a", verbose = FALSE)
  expect_equal(fit, structure(list(f0 = 10.3651999807438, f0.se = 7.74617801784988, 
                                   f0.p = 0.283145377159807, C0 = 328.932444530028, 
                                   weights = c(1, 0.224712479154538, 1, 1), diagnostics = ""), 
                              .Names = c("f0", "f0.se", "f0.p", "C0", "weights", "diagnostics")), tolerance = 1e-5)
  
  C <- 320 + 0:3 * 10
  fit <- rlin.fit(t, C, 1, 0.3, "a", verbose = FALSE)
  expect_equal(fit, structure(list(f0 = 9, f0.se = 5.05560019475045e-16, f0.p = 3.15544362088405e-33, 
                                   C0 = 320, weights = c(1, 1, 1, 1), diagnostics = ""), 
                              .Names = c("f0", "f0.se", "f0.p", "C0", "weights", "diagnostics")))
})

test_that("HMR regression gives expected result", {
  skip_on_cran() #slightly different SE/p value on winbuilder
                 #doesn't return NA on winbuilder
  
  t <- c(0, 1/3, 2/3, 1, 1.2, 1.3)
  C <- c(320, 341, 352, 359, 360, 360)
  fit <- HMR.fit(t, C, 1, 0.3, "a", verbose = FALSE)
  
  expect_equal(fit, list(f0 = 26.1229339419655, f0.se = 1.60625294868361, f0.p = 0.000505786166225468,
                         kappa = 1.97015450227329, phi = 364.122024254071, AIC = 17.4284243644314, 
                         AICc = 57.4284243644314, RSE = 0.750764759562218, diagnostics = ""), tolerance = 1e-7)
 
  C[2] <- 500
  fit <- HMR.fit(t, C, 1, 0.3, "a", verbose = FALSE)
  expect_equal(fit$diagnostics, "singular gradient")
  
  t <- c(0, 1/3, 2/3, 1)
  C <- 320 + 0:3 * 10
  fit <- HMR.fit(t, C, 1, 0.3, "a", verbose = FALSE)
  expect_equal(fit$diagnostics, "Missing value or an infinity produced when evaluating the model")
  
})

test_that("NDFE regression gives expected result", {
  skip_on_cran() #slightly different flux/SE/p/tau value on winbuilder
  t <- c(0, 1/3, 2/3, 1)
  C <- c(320, 341, 352, 359)
  fit <- NDFE.fit(t, C, 1, 0.3, "a", verbose = FALSE)
  expect_equal(fit, structure(list(f0 = 108.747818996036, f0.se = 62.9208765449726, 
                                   f0.p = 0.333927030887926, tau = 0.0111016600093449, C0 = 319.992148880558, 
                                   AIC = 10.7342239838649, AICc = NA_real_, RSE = 0.681122330960559, 
                                   diagnostics = ""), .Names = c("f0", "f0.se", "f0.p", "tau", 
                                                                 "C0", "AIC", "AICc", "RSE", "diagnostics")), tolerance = 1e-5)
  
  C[2] <- 400
  fit <- NDFE.fit(t, C, 1, 0.3, "a", verbose = FALSE)
  expect_equal(fit, structure(list(f0 = NA_real_, f0.se = NA_real_, f0.p = NA_real_, 
                                   tau = NA_real_, C0 = NA_real_, AIC = NA_real_, AICc = NA_real_, 
                                   RSE = NA_real_, diagnostics = "Missing value or an infinity produced when evaluating the model"), 
                              .Names = c("f0", "f0.se", "f0.p", "tau", "C0", "AIC", "AICc", "RSE", "diagnostics"
                                   )))
  
  C <- 320 + 0:3 * 10
  fit <- NDFE.fit(t, C, 1, 0.3, "a", verbose = FALSE)
  expect_equal(fit, structure(list(f0 = NA_real_, f0.se = NA_real_, f0.p = NA_real_, 
                                   tau = NA_real_, C0 = NA_real_, AIC = NA_real_, AICc = NA_real_, 
                                   RSE = NA_real_, diagnostics = "step factor 5.82077e-11 reduced below 'minFactor' of 1e-10"), 
                              .Names = c("f0", "f0.se", "f0.p", "tau", "C0", "AIC", "AICc", "RSE", "diagnostics"
                                   )))
})