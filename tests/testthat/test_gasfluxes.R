library(gasfluxes)
context("gasfluxes checks")

test_that("gasfluxes returns expected values", {
  skip_on_cran() #less NA values on winbuilder
  data(fluxMeas)
  res <- gasfluxes(fluxMeas[1:20,], 
                   .id = "serie", .V = "V", .A = "A",
                   .times = "time", .C = "C", 
                    methods=c("HMR", "original HMR"),verbose=FALSE, plot = FALSE)
  expect_equal(data.table:::setDF(res), structure(list(serie = c("ID1", "ID2", "ID3", "ID4", "ID5"), 
                                                            HMR.f0 = c(NA, NA, -0.123566162944742, NA, 0.0418403302775212), 
                                                            HMR.f0.se = c(NA, NA, 0.39896371493858, NA, 135759548593600), 
                                                            HMR.f0.p = c(NA, NA, 0.808792404502918, NA, 1), 
                                                            HMR.kappa = c(NA, NA, 2.39686189586586, NA, 6.26517905920712e-06), 
                                                            HMR.phi = c(NA, NA, 0.383169270087676, NA, 12445.2133007081), 
                                                            HMR.AIC = c(NA, NA, -6.62832670015326, NA, -19.0235445653063), 
                                                            HMR.AICc = c(NA, NA, -46.6283267001533, NA, -59.0235445653063), 
                                                            HMR.RSE = c(NA, NA, 0.0777441330321253, NA, 0.0165109014333297), 
                                                            HMR.diagnostics = c("singular gradient", "Missing value or an infinity produced when evaluating the model",  
                                                                                "", "Missing value or an infinity produced when evaluating the model", ""), 
                                                            original.HMR.f0 = c(0.0555669866740785, -0.0611626095687462, -0.123568879808343, 0, 0.0418402102972429), 
                                                            original.HMR.f0.se = c(0.028697055158409, 0.0263636317002555, 0.100299647826534, NA, 0.00840546056899412), 
                                                            original.HMR.f0.p = c(0.19245050979196, 0.146138869452879, 0.343139597683113, NA, 0.0380688306059212), 
                                                            original.HMR.kappa = c(NA, NA, 2.39693676099531, NA, NA), 
                                                            original.HMR.phi = c(NA, NA, 0.383170338419227,NA, NA), original.HMR.AIC = c(NA, NA, -6.62832670016226, NA, NA), 
                                                            original.HMR.AICc = c(NA, NA, -46.6283267001623, NA, NA), 
                                                            original.HMR.RSE = c(NA, NA, 0.0777441330320377, NA, NA), 
                                                            original.HMR.diagnostics = c("", "", "", "", "")), 
                                                       .Names = c("serie", "HMR.f0", "HMR.f0.se", "HMR.f0.p", "HMR.kappa", "HMR.phi", "HMR.AIC", 
                                                                  "HMR.AICc", "HMR.RSE", "HMR.diagnostics", "original.HMR.f0", 
                                                                  "original.HMR.f0.se", "original.HMR.f0.p", "original.HMR.kappa", 
                                                                  "original.HMR.phi", "original.HMR.AIC", "original.HMR.AICc", 
                                                                  "original.HMR.RSE", "original.HMR.diagnostics"), row.names = c(NA, -5L), class = "data.frame"))
  res <- gasfluxes(fluxMeas[1:20,], 
                   .id = "serie", .V = "V", .A = "A",
                   .times = "time", .C = "C", 
                   select = "RF2011",verbose=FALSE, plot = FALSE)
  
  expect_equal(data.table:::setDF(res), structure(list(serie = c("ID1", "ID2", "ID3", "ID4", "ID5"), 
                                                       linear.f0 = c(0.0555669866740785, -0.0611626095687461, -0.0428897148924301, -0.0111793261413552, 0.041840210297243), 
                                                       linear.f0.se = c(0.028697055158409, 0.0263636317002555, 0.0425127107690505, 0.0599487242659592, 0.00840546056899412), 
                                                       linear.f0.p = c(0.19245050979196, 0.146138869452879, 0.419251528032237, 0.86926934170296, 0.0380688306059213), 
                                                       linear.C0 = c(0.308970790118557, 0.475855538008889, 0.465418448020012, 0.385341651037846, 0.333253184107449), 
                                                       linear.AIC = c(-10.9887845537839, -11.6825495623863, -8.01509608282695, -5.16760526792864, -21.0235708278954), 
                                                       linear.AICc = c(Inf, Inf, Inf, Inf, Inf), 
                                                       linear.RSE = c(0.040927093104423, 0.0375274115555534, 0.0593530390736887, 0.0847274536777346, 0.0116749320402043), 
                                                       linear.diagnostics = c("", "", "", "", ""), 
                                                       robust.linear.f0 = c(0.0555669866740785, -0.0611626095687461, -0.0303109387388622, 0.00681856647255265, 0.041840210297243), 
                                                       robust.linear.f0.se = c(0.028697055158409, 0.0263636317002555, 0.0153424660030732, 0.0204971767735417, 0.00840546056899412), 
                                                       robust.linear.f0.p = c(0.19245050979196, 0.146138869452879, 0.163090704934176, 0.750421672669818, 0.0380688306059213), 
                                                       robust.linear.C0 = c(0.308970790118557, 0.475855538008889, 0.473272206610976, 0.339838564756016, 0.333253184107449), 
                                                       robust.linear.weights = c("1,1,1,1", "1,1,1,1", "1,1,0.16,1", "1,0.15,1,1", "1,1,1,1"), 
                                                       robust.linear.diagnostics = c("", "", "", "", ""), 
                                                       original.HMR.f0 = c(0.0555669866740785, -0.0611626095687462, -0.123568879808343, 0, 0.0418402102972429), 
                                                       original.HMR.f0.se = c(0.028697055158409, 0.0263636317002555, 0.100299647826534, NA, 0.00840546056899412), 
                                                       original.HMR.f0.p = c(0.19245050979196, 0.146138869452879, 0.343139597683113, NA, 0.0380688306059212), 
                                                       original.HMR.kappa = c(NA, NA, 2.39693676099531, NA, NA), 
                                                       original.HMR.phi = c(NA, NA, 0.383170338419227, NA, NA), 
                                                       original.HMR.AIC = c(NA, NA, -6.62832670016226, NA, NA), 
                                                       original.HMR.AICc = c(NA, NA, -46.6283267001623, NA, NA), 
                                                       original.HMR.RSE = c(NA, NA, 0.0777441330320377, NA, NA), 
                                                       original.HMR.diagnostics = c("", "", "", "", ""), 
                                                       flux = c(0.0555669866740785, -0.0611626095687461, -0.0303109387388622, 0.00681856647255265, 0.041840210297243), 
                                                       flux.se = c(0.028697055158409, 0.0263636317002555, 0.0153424660030732, 0.0204971767735417, 0.00840546056899412), 
                                                       flux.p = c(0.19245050979196, 0.146138869452879, 0.163090704934176, 0.750421672669818, 0.0380688306059213), 
                                                       method = c("robust linear", "robust linear", "robust linear", "robust linear", "robust linear")), 
                                                  .Names = c("serie", "linear.f0", "linear.f0.se", "linear.f0.p", "linear.C0", "linear.AIC", "linear.AICc", 
                                                             "linear.RSE", "linear.diagnostics", "robust.linear.f0", "robust.linear.f0.se", "robust.linear.f0.p", 
                                                             "robust.linear.C0", "robust.linear.weights", "robust.linear.diagnostics","original.HMR.f0", 
                                                             "original.HMR.f0.se", "original.HMR.f0.p", "original.HMR.kappa", "original.HMR.phi", "original.HMR.AIC", 
                                                             "original.HMR.AICc", "original.HMR.RSE", "original.HMR.diagnostics", "flux", "flux.se", "flux.p", "method"), 
                                                  row.names = c(NA, -5L), class = "data.frame"))
})



