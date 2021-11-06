#scenario 14 -idw just predict to station lcoations
#has been checked 
require(data.table)
require(dplyr)
library(sp)
library(raster)
library(fields)
library(rlist)
library(gstat) 

######################################################################################################
#days: vector with values 1:31 
#months: int vector w/ values 1-12
#VAR: tmax, tmin, or precip. Specify one and specll correctly.
#years: number of years to consider--1980-2014 for historical data
#example: Bayesian_Driver(days = 1:31, months = 1:12, n.samples = 5000,
#######################################################################################################
#train on WRF , and follow same method as scenario 001 except use IDW interpolation to predict to station 
#locations. Also take into account elevation by first bringing meteorological variable to reference elevation 
#(200m) and then backtransforming. After this, you should do bias correction. And then run scenario 015 (Cross val)
#or scenario 016 (predict to prediction grid)
Model_Driver = function(days=days, months=months, 
                           VAR, MY_YEAR=MY_YEAR, 
                           test_data_list=NULL, df=df, prog_log = prog_log, lapse_params=NULL, test_stations=NULL) {
 
  lcc_test_grid_raster = test_data_list
  
  print(paste("running year ", MY_YEAR))

  #create WRF data for a year for VAR (use historical = FALSE for future projections)
  wrfdat80 = get_wrf_data(MY_YEAR, VAR, historical = FALSE)
  
  #create station data for a year for VAR
  d3long = get_station_data(MY_YEAR, VAR, df)
  
  #list to hold results for predictions to station locations
  preds_list = list()
  

  if (VAR == "precip"){
    var_name = "precipdata"
    Lapse_rate = lapse_params$precipLR
  }
  
  if (VAR == "tmax"){
    var_name = "tmaxdata"
    Lapse_rate = lapse_params$tmaxLR
  }
  
  if (VAR == "tmin"){
    var_name = "tmindata"
    Lapse_rate = lapse_params$tminLR
  }
  
  
  climateDF_stations = d3long
  climateDF_wrf = wrfdat80
  
  #loop over months
  for (k in 1:length(months)){
    cur_data_station = filter(climateDF_stations, MONTH == months[k])
    cur_data_wrf = filter(climateDF_wrf, MONTH == months[k])
    
    #loop over days 
    for (i in 1:length(days)){
      
      mydata_station = filter(cur_data_station, DAY == days[i])
      
      var_of_interest = names(mydata_station)[4]
      
      mydata_wrf = filter(cur_data_wrf, DAY == days[i])
      
      cat("month is: ", mydata_station$MONTH[1], "\n", 
          "day is: ", mydata_station$DAY[1],
          file = prog_log, append = TRUE)
      
      
      if (dim(mydata_station)[1] == 0){
        print("there are no observations! Skipping this day.")
        next
      }
      
      
      #skip december 31 for leap years WRF historical/ Jan 29 for future rcp 95 driven data
      if (dim(mydata_wrf)[1] ==0){
        print("No obs for WRF for this day. skipping!")
        next
      }
      
      
      #remove crazy ridiculous values if present
      if (var_of_interest %in% c("TMAX", "TMIN")){
        mydata_station = mydata_station[mydata_station[[var_of_interest]] < 100, ]
      }
      
      if (var_of_interest == "PRCP"){
        mydata_station = mydata_station[mydata_station[[var_of_interest]] >= 0, ]
      }
      
      
      #preprocess wrf/station data:
      #WRF and station data will be returned as data.frames w/ LCC coords
      #transformation to reference elevation occurs in this function
      #a list of dataframes w/ lat/lon and LCC coords is returned (4 data.frames total)
      climate_data_list = preprocess_climate_data_train(mydata_station, mydata_wrf, VAR="tmin", Lapse_rate)
      
      #IDW interpolation to station locations; Tref is temp at reference elevation
      #Tref is created in preprocess_climate_data_train()
      #lcc_test_grid_raster is the 1km prediction grid in lcc coords
      results_dat = IDW_interpolation(climate_data_list, var_of_interest = "Tref",
                                      lcc_test_grid_raster, VAR, Lapse_rate)
      
      
      res = list(results_dat)
      
      # res will be a large list composed of ind. dataframes for each day. Dataframe 
      #will include an extra column 'test.pred.mu' which are the interpolated values from IDW interpolation
      #this is identical to scenario001 except we are using IDW interpolation instead of bayesian kriging to 
      #make predictions at station locations
      preds_list = c(preds_list, res)
    }
  }
  
  print("all done")
  results_dir = create_results_dir()
  saveRDS(preds_list, file = paste(results_dir, "/Year_", MY_YEAR, var_name, "IDW.Rds", sep = ""))
}


#test.pred.mu
####################################################################################
#function to inverse distance weighting to interpolate to test grid
#also, extract the predictions at station locations to validate method
####################################################################################
#testgrid lat/lon is the prediction grid on lat/lon
#test_data_list is the prediciton grid on lcc coords
#both include elevation -- ELEV
IDW_interpolation = function(climate_data_list, var_of_interest = "Tref",
                             lcc_test_grid_raster, VAR, Lapse_rate) {
  
  #WRFdf is wrf data, test is Stationdf data (LCC coordinate system)
  #these are both data.frames w/ LCC coords, NOT spatial objects
  mydata_wrf.lccDF = climate_data_list$mydata_wrf.lccDF
  mydata_station.lccDF = climate_data_list$mydata_station.lccDF
  
  mydata_station_latlonDF = climate_data_list$mydata_station_latlonDF
  mydata_wrfDF_latlon = climate_data_list$mydata_wrfDF_latlon
  
  #want to use mydata_wrf.lccDF to make interpolated projections to mydata_station.lccDF
  #Tref is the temp at reference elevation
  myformula = formula(paste(var_of_interest,"~ 1"))
  interp_to_stations = idw(formula = myformula,  locations = ~LON + LAT, data = mydata_wrf.lccDF, 
                           newdata = mydata_station.lccDF, idp=2, nmax=9) #idw
  
  #Now, BACKTRANSFORM to get actual elevation-adjusted temperature
  interp_to_stations$ELEV = mydata_station.lccDF$ELEV
  interp_to_stations$var1.var = NULL
  
  #apply back transformation to get the actual elevation-adjusted value for VAR
  test.pred.mu = back_transform(VAR, Lapse_rate, interp_to_stations)
  
  mydata_station_latlonDF$test.pred.mu = test.pred.mu
  LCC_coords = mydata_station.lccDF[ ,c("LON", "LAT")]
  names(LCC_coords) = c("LON_LCC", "LAT_LCC")
  results_dat = cbind(mydata_station_latlonDF, LCC_coords)
  return(results_dat)
}

###################################################################################
#import_test_grid -- not necessary for this scenario.
###################################################################################
#don't need this function here, so it's empty.
#We will use it when we train on bias corrected data
import_test_grid = function(){

}



###################################################################################
#import daily station data
#######################################################################################
#import overall D2 daily data
#will load as 'df'. 
import_daily_station_data <- function(VAR = NULL){
  e <- new.env()
  load("/raid/mholthui/fromBabbage/data/GHCND/daily_data2.RData", envir = e)
  df = get("df", envir =  e)
  return(df)
}

#make progress log
make_progress_log = function(results_dir){
  prog_log = paste(results_dir, "/prog.txt", sep = "")
  return(prog_log)
}

#create directory to put results in -- change this so it works for pascal
create_results_dir = function(){

  if(!dir.exists("/raid/mholthui/fromBabbage/IDW_scenario014tmin_RCP45_TOGHCND")){
    dir.create("/raid/mholthui/fromBabbage/IDW_scenario014tmin_RCP45_TOGHCND")
  }
  return("/raid/mholthui/fromBabbage/IDW_scenario014tmin_RCP45_TOGHCND")
}



################################################################
##preprocess_climate_data_train
## transforms WRF data to lcc coordinates so we can interpolate
#returns original WRF/station dataframes, as well as dataframes w/ 
#LCC coordinates (note-these are not spatial objects)
################################################################
preprocess_climate_data_train = function(mydata_station, mydata_wrf, VAR = "tmin", Lapse_rate){
  mydata_station_latlonDF = mydata_station
  mydata_wrfDF_latlon = mydata_wrf
  
  #apply baseline tranformation to get temperature at reference elevation
  mydata_wrf = reference_transform(mydata_wrf, VAR, Lapse_rate)


  coordinates(mydata_wrf) <- c("LON", "LAT")
  proj4string(mydata_wrf) <- CRS("+init=epsg:4326") # set initial projection to lat/lon, WGS 84
  
  #this is the correct lambert conformal proj4string
  CRS.new <- CRS("+proj=lcc
                 +lat_0=40
                 +lat_1=44.09
                 +lat_2=45.15
                 +lon_0=-74.0
                 +x_0=0
                 +y_0=0
                 +a=6370
                 +b=6370
                 +units=m
                 +no_defs")
  mydata_wrf.lcc = spTransform(mydata_wrf, CRS.new)
  mydata_wrf.lccDF = as.data.frame(mydata_wrf.lcc)
  
  #also transform station data #this is for validation purposes only
  coordinates(mydata_station)  = ~ LON +LAT
  proj4string(mydata_station) =  CRS("+init=epsg:4326")
  
  #transfrom to lcc projection
  mydata_station = spTransform(mydata_station, CRS.new)
  mydata_station.lccDF = as.data.frame(mydata_station)
  
  #return spatial points dataframes of wrf and station data (lcc transformed), plus month, day and year
  #plus original WRF/station datasets--denoted w/ DF (lat/lon dataframes)
  return(list(mydata_station.lccDF = mydata_station.lccDF, mydata_wrf.lccDF = mydata_wrf.lccDF,
              mydata_station_latlonDF = mydata_station_latlonDF,
              mydata_wrfDF_latlon = mydata_wrfDF_latlon))
}

#we do not need this function--no need for prediction grid yet.
preprocess_climate_data_test = function(test_set = NULL){
}

#function to backtransform from reference elevation value to actual elevation
back_transform = function(VAR, Lapse_rate, interp_to_stations){
  if (VAR == "tmax"){
    interp_to_stations$TMAX_interp_corr = interp_to_stations$var1.pred - Lapse_rate*(interp_to_stations$ELEV - 200)
    test.pred.mu = interp_to_stations$TMAX_interp_corr
  }
  
  if (VAR == "tmin"){
    interp_to_stations$TMIN_interp_corr = interp_to_stations$var1.pred - Lapse_rate*(interp_to_stations$ELEV - 200)
    test.pred.mu = interp_to_stations$TMIN_interp_corr
  }
  
  if (VAR == "precip"){
    interp_to_stations$PRECIP_interp_corr = interp_to_stations$var1.pred*((1+Lapse_rate*(interp_to_stations$ELEV - 200)) /
                                      (1-Lapse_rate*(interp_to_stations$ELEV - 200)))
    test.pred.mu = interp_to_stations$PRECIP_interp_corr
  }
  return(test.pred.mu)
}

#apply baseline tranformation to get temperature at reference elevation
reference_transform = function(mydata_wrf, VAR, Lapse_rate){
  if (VAR == "tmax"){
    
    mydata_wrf$Tref = mydata_wrf$TMAX - Lapse_rate*(200 - mydata_wrf$ELEV)
  }
  
  if (VAR == "tmin"){
    mydata_wrf$Tref = mydata_wrf$TMIN - Lapse_rate*(200 - mydata_wrf$ELEV)
  }
  
  if (VAR == "precip"){
    mydata_wrf$Tref = mydata_wrf$PRCP *((1 + Lapse_rate*(200 - mydata_wrf$ELEV))/ 
                                          (1 - Lapse_rate*(200 - mydata_wrf$ELEV)))
  }
  
  return(mydata_wrf)
}

