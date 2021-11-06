require(data.table)
require(dplyr)
library(spNNGP)
library(spBayes)
library(sp)
library(MBA)
library(fields)


####################################################################################
#function to run daily wrf data and predict onto station data.
#Do not use this with prediction grid!
####################################################################################
run_Bayesian_model_wrf = function(climate_data_list, test_data_list = NULL, n.samples = 5000,
                                  var_of_interest, logFile, logFile2) {
  predCoords = climate_data_list$predCoords
  predCoords_orig = climate_data_list$predCoords_orig
  predCovars = climate_data_list$predCovars
  mydata_station = climate_data_list$mydata_station
  wrf_coords_original = climate_data_list$wrf_coords_original
  wrf_coords = climate_data_list$wrf_coords
  mydata_wrf = climate_data_list$mydata_wrf
  spLM_coords_wrf = climate_data_list$spLM_coords_wrf
  d.max = climate_data_list$d.max
  
  YEAR = mydata_wrf$YEAR[1]
  DAY = mydata_wrf$DAY[1]
  MONTH = mydata_wrf$MONTH[1]
  logFile = logFile
  logFile2 = logFile2
  n.samples=n.samples
  
  #formula for model
  myformula = formula(paste(var_of_interest,"~ ELEV"))
  
  #fit daily model using exponential covariance function
  spatial_model = spNNGP::spNNGP(myformula,
                                 data = mydata_wrf, coords = spLM_coords_wrf,
                                 starting = list("phi" = 3/100, "sigma.sq" = 5, "tau.sq"= .5),
                                 tuning = list("phi" = 1, "sigma.sq" = 1, "tau.sq" = 1), n.neighbors = 15,
                                 priors = list("phi.Unif" = c(3/d.max, 3/10), "sigma.sq.IG" = c(2, 20), "tau.sq.IG" = c(2, 1)),
                                 cov.model = "exponential", n.samples = n.samples, method="sequential", n.omp.threads=4,
                                 search.type = "brute")
  
  burn.in <- floor(0.75*n.samples):n.samples
  #get posterior median and SD for spatial effects
  w.samples = spatial_model$p.w.samples
  w.hat.mu <- apply(spatial_model$p.w.samples[,burn.in], 1, median)
  w.hat.sd = apply(spatial_model$p.w.samples[,burn.in], 1, sd)
  
  #posterrio median and SD for betas (intercept and elevation coefficient)
  beta.samples = spatial_model$p.beta.samples
  beta.hats = apply(spatial_model$p.beta.samples[ burn.in, ], 2, median)
  
  #this is the spPredict version for spNNGP (arguments are slightly different for spBayes)
  sp_preds_test = spNNGP::spPredict(spatial_model, start = floor(0.75*n.samples), 
                                    thin = 2, X.0 = predCovars, coords.0 = predCoords, n.omp.threads = 4)
  
  #get posterior means
  test.pred.mu = apply(sp_preds_test$p.y.0, 1, mean) 
  
  #get posterior sd's
  test.pred.sd = apply(sp_preds_test$p.y.0, 1, sd)
  
  #calculate RMSE on test data (e.g. station data)
  MSE = mean((mydata_station[[var_of_interest]] - test.pred.mu)^2)
  rmse = round(sqrt(MSE), 2)
  
  #save rmse values to logfile
  cat("RMSE: Y: ", YEAR, " M: ", MONTH, " DAY ", DAY, "=" , rmse, "\n\n", file = logFile, append = TRUE)
  
  #save coeffs for intercept and elev to logfile2
  capture.output(print(beta.hats, print.gap=3), file=logFile2, append = TRUE)
  
  results_dat = data.frame(cbind(mydata_station, test.pred.mu))
  
  
  return(results_dat)
}

#############################################################################################
#Run Bayesian spatial model on EITHER WRF/Station data as training data, then predict onto 
#prediction grid
##############################################################################################

run_Bayesian_model_1km_test = function(climate_data_list, test_data_list, n.samples = 5000,
                                  var_of_interest, logFile, logFile2) {
  #get prediction coordinates and covariates
  predCoords = test_data_list$predCoords
  predCoords_orig = test_data_list$predCoords_orig
  predCovars = test_data_list$predCovars

  #mydata is training data--either WRF or station
  mydata = climate_data_list$mydata
  mydata_coords_original = climate_data_list$coords_original
  spLM_coords_wrf = climate_data_list$spLM_coords
  d.max = climate_data_list$d.max
  
  YEAR = mydata$YEAR[1]
  DAY = mydata$DAY[1]
  MONTH = mydata$MONTH[1]
  logFile = logFile
  logFile2 = logFile2
  n.samples = n.samples
  
  #formula for model . ELEV is the only covariate
  myformula = formula(paste(var_of_interest,"~ ELEV"))
  
  #fit daily model using exponential covariance function
  spatial_model = spNNGP::spNNGP(myformula,
                                 data = mydata, coords = spLM_coords,
                                 starting = list("phi" = 3/100, "sigma.sq" = 5, "tau.sq"= .5),
                                 tuning = list("phi" = 1, "sigma.sq" = 1, "tau.sq" = 1), n.neighbors = 15,
                                 priors = list("phi.Unif" = c(3/d.max, 3/10), "sigma.sq.IG" = c(2, 20), "tau.sq.IG" = c(2, 1)),
                                 cov.model = "exponential", n.samples = n.samples, method="sequential", n.omp.threads=4,
                                 search.type = "brute")
  
  burn.in <- floor(0.75*n.samples):n.samples
  #get posterior median and SD for spatial effects
  w.samples = spatial_model$p.w.samples
  w.hat.mu <- apply(spatial_model$p.w.samples[,burn.in], 1, median)
  w.hat.sd = apply(spatial_model$p.w.samples[,burn.in], 1, sd)
  
  #posterrio median and SD for betas (intercept and elevation coefficient)
  beta.samples = spatial_model$p.beta.samples
  beta.hats = apply(spatial_model$p.beta.samples[ burn.in, ], 2, median)
  
  #this is the spPredict version for spNNGP (arguments are slightly different for spBayes)
  #predict to 1km grid
  sp_preds_test = spNNGP::spPredict(spatial_model, start = floor(0.75*n.samples), 
                                    thin = 2, X.0 = predCovars, coords.0 = predCoords, n.omp.threads = 4)
  
  #get posterior means
  test.pred.mu = apply(sp_preds_test$p.y.0, 1, mean) 
  
  #get posterior sd's
  test.pred.sd = apply(sp_preds_test$p.y.0, 1, sd)
  
  #calculate RMSE on test data (e.g. station data)
  MSE = mean((mydata_station[[var_of_interest]] - test.pred.mu)^2)
  rmse = round(sqrt(MSE), 2)
  
  #save rmse values to logfile
  cat("RMSE: Y: ", YEAR, " M: ", MONTH, " DAY ", DAY, "=" , rmse, "\n\n", file = logFile, append = TRUE)
  
  #save coeffs for intercept and elev to logfile2
  capture.output(print(beta.hats, print.gap=3), file=logFile2, append = TRUE)
  
  results_dat = data.frame(cbind(mydata_station, test.pred.mu))
  
  
  return(results_dat)
}