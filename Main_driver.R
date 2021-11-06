#!/usr/bin/Rscript
args<-commandArgs(trailingOnly=TRUE)

#Functions to create WRF/station yearly data
source("/raid/mholthui/fromBabbage/NewRuns/data_construction.R")
print("sourced the data")
#####################################################################################
#to train on WRF and predict to stations use:
#source("scenario001-UseWRFandStation.R")
#to train on WRF or stations use"
#you must specify train_data = "stations" or "wrf". Default is "stations"
#source("scenario002-Teston1kmGrid.R")
#source('scenario001-UseWRFandStation.R')
#scenario 4 is training on bias corrected data, predict to 1km grid
#scenario = "/data/mholthui/R_parallel/NewRuns/scenario004-Train_on_WRF_predictions.R"
#scenario = "/data/mholthui/R_parallel/NewRuns/scenario004-Train_on_WRF_preds_crossval.R"
#######################################################################################

start_time <- Sys.time()
main_driver = function(VAR, days = 1:31, months = 1:12, years=args[1], scenario, train_data = NULL){

	#specify which scenario to use
	source(scenario)
  #print("train data is ", train_data)
  print("inside main driver")
	#preprocess grid data ONLY ONCE (it will never change)
	#results for prediction coords, covariates are stored in 
  if (scenario == "/raid/mholthui/fromBabbage/NewRuns/scenario005-Train_on_WRF_preds_testset.R" ||
  scenario == "/raid/mholthui/fromBabbage/NewRuns/scenario005-Train_on_WRF_preds_testset_NOELEV.R"){
    print("doing cross validation")
	  test_set = import_test_grid(VAR, k=1)
  }
  else{
    test_set = import_test_grid()
  }
	test_data_list = preprocess_climate_data_test(test_set = test_set)
	
 #import data. Type of data will vary depending on scenario
	df = import_daily_station_data(VAR)

	
	#directory to dump results into
	results_dir = create_results_dir()
	
	#make progress log
	prog_log = make_progress_log(results_dir)

	#directory to dump results into
	create_results_dir()


  	for(y in 1:length(years)) {
      MY_YEAR = as.integer(years[y])
  
  	  cat("Year: ",  MY_YEAR ,"\n", file = prog_log, append = TRUE)
  
  	  Bayesian_Driver(days=days, months = months, n.samples=5000,
     VAR = VAR, MY_YEAR = MY_YEAR, train_data = train_data, test_data_list = test_data_list, 
     df=df, prog_log = prog_log)
  
  	}
}

#scenario = "/data/mholthui/R_parallel/NewRuns/scenario004-Train_on_WRF_predictions.R" #train on bias corrected output from scenario 1, predict to 1km grid
#scenario = "/data/mholthui/R_parallel/NewRuns/scenario001-UseWRFandStation.R" #train on wrf, predict to stations
#scenario = "/data/mholthui/R_parallel/NewRuns/scenario006-Train_on_WRF_preds_crossval.R"
#scenario = "/data/mholthui/R_parallel/NewRuns/scenario005-Train_on_WRF_preds_testset.R"
#scenario = "/data/mholthui/R_parallel/NewRuns/scenario001B-NBkrig_to_1km_grid.R"
scenario = "/raid/mholthui/fromBabbage/NewRuns/scenario005-Train_on_WRF_preds_testset_NOELEV.R"
VAR  = "tmax"
main_driver(VAR = "tmax", scenario = scenario)
end_time <- Sys.time()

totaltime = end_time - start_time
totaltime
resfile = "mytimesBayesCV.txt"
compTime = as.numeric(totaltime, units = "mins")
cat(compTime,"\n",  file = resfile, append = TRUE)