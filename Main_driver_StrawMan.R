#!/usr/bin/Rscript
args<-commandArgs(trailingOnly=TRUE)

#Functions to create WRF/station yearly data
source("/raid/mholthui/fromBabbage/NewRuns/data_construction.R")

#####################################################################################
#driver for various straw man (IDW) downscaling options
#######################################################################################

start_time <- Sys.time() 
main_driver = function(VAR, days = 1:31, months = 1:12, years=args[1], scenario, lapse_params = NULL){
  
  #specify which scenario to use
  source(scenario)
  #print("train data is ", train_data)
  print("inside main driver")
  print(VAR)
  
  if(scenario == "/raid/mholthui/fromBabbage/NewRuns/scenario017-IDW_test_setV2.R" ||
  scenario == "/raid/mholthui/fromBabbage/NewRuns/scenario019-IDW_WRF_to_test_locations.R" || 
  scenario == "/raid/mholthui/fromBabbage/NewRuns/scenario017-IDW_test_setV2_NOELEV.R"){
  print("doing crossval")
    test_stations = import_test_set(VAR, k=5)
  }
  
  #import 1km test grid. 
  test_set = import_test_grid()
 
  #Convert to raster and transform to LCC coord system
  test_data_list = preprocess_climate_data_test(test_set = test_set)
  
  #import data. Type of data will vary depending on scenario
  df = import_daily_station_data(VAR)
  
  
  #directory to dump results into
  results_dir = create_results_dir()
  
  #make progress log
  prog_log = make_progress_log(results_dir)

  for(y in 1:length(years)) {
    MY_YEAR = as.integer(years[y])
    
    cat("Year: ",  MY_YEAR ,"\n", file = prog_log, append = TRUE)
    
    Model_Driver(days=days, months = months, 
                    VAR = VAR, MY_YEAR = MY_YEAR, test_data_list = test_data_list, 
                    df=df, prog_log = prog_log, lapse_params = lapse_params, test_stations = test_stations)
    
  }
  
}

#scenario = "/raid/mholthui/fromBabbage/NewRuns/scenario014-IDW_to_station.R"
scenario = "/raid/mholthui/fromBabbage/NewRuns/scenario018-IDW_WRF_to_1km_grid.R"

#use this for future preds  to stations only! (can also use for tmax / tmin...but look at function. Need to specify GHCN station locations w/ separate dataset in csv file in /Data/GHCN
#scenario = "/raid/mholthui/fromBabbage/NewRuns/scenario-014_IDW_to_station_FuturePRCP.R"
VAR  = "precip"

#use for precip: use bias corrected precip at station locations --> interpolate to 1km grid
#scenario = "/raid/mholthui/fromBabbage/NewRuns/scenario016-IDW_to_pred_grid.R"
#VAR = "precip"

#these are lapse rates for IDW, Leave OUT if you are using any other model than IDW
#there were my original lapse rates. Going to go w/ Winter's lapse rates for RCP85 downscaling
#lapse_params = list("tmaxLR" = 0.0054899011, "tminLR" = 0.0050204795, "precipLR" = 0.0002421679)

#I just used the ones from Winter et al. for RCP85 downscaling
lapse_params = list("tmaxLR" = 0.00592, "tminLR" = 0.00485, "precipLR" = 0.00025)
main_driver(VAR = "precip", scenario = scenario, lapse_params = lapse_params)

