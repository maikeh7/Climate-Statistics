#scenario 16 -idw predict to grid based on bias corrected predictions to station data
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
#years: number of years to consider. Min is 1980, max is 1990.
#example: Bayesian_Driver_wrf(days = 1:31, months = 1:12, n.samples = 5000,
#######################################################################################################
#here, interpolate from bias corrected interpolated values at station locations to 1km grid
Model_Driver = function(days=days, months=months, 
                           VAR, MY_YEAR=MY_YEAR, 
                           test_data_list=NULL, df=df, prog_log = prog_log, lapse_params, test_stations=NULL) {
  
  print(paste("running year ", MY_YEAR))
  
  #this is the prediction raster grid
  lcc_test_grid_raster = test_data_list
  
  #import lat/lon test grid b/c we want that as well
  #this is just a data.frame, NOT a raster
  testgrid_LCC_DF = import_test_grid_LCC_DF()
  
  #import bias corrected data
  dat = filter(df, YEAR == MY_YEAR)
  
  #list to hold results for prediction to 1km grid
  preds_list = list()
  
  #list to hold results for predictions to station locations
  preds_list_stations = list()
  
  
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
  
  
  for (i in 1:length(months)){
    monthly_data = filter(dat, MONTH == i)
    
    #print(dim(monthly_data))
    for (j in 1:length(days)){
      print(paste("day = ", j, " Month = ", i, "\n"))
      mydata_station = filter(monthly_data, DAY == j)
      
      if (dim(mydata_station)[1] == 0){
        next
      }
        
        #transform to baseline elevation in this function  
        #transformed variable will be called 'Tref'
        #we already had lat/lon LCC lat/lon coords so no need to add those again
        climate_data_list = preprocess_climate_data_train(mydata_station, lcc_test_grid_raster, VAR, Lapse_rate)
        
        #IDW interpolation to station locations; Tref is temp at reference elevation
        #Tref is created in preprocess_climate_data_train()
        results_dat = IDW_interpolation(climate_data_list, var_of_interest = "Tref", testgrid_LCC_DF,
                                              VAR,  Lapse_rate)
        
        
        
        resGrid = list(results_dat$interp_data_latlon) #predictions to 1km grid (lat/lon + LCC coords) -dataframe
        resStations = list(results_dat$station_preds) #predictions to stations - dataframe 
        
        #of predictions. res will be a large list composed of ind. dataframes for each day. Dataframe 
        #will include an extra column 'test.pred.mu' which are the predictions from bayesian kriging
        preds_list = c(preds_list, resGrid)
        preds_list_stations = c(preds_list_stations, resStations)
    }
  }
    
  results_dir = create_results_dir()[1]
  results_dir_test = create_results_dir()[2]
  saveRDS(preds_list, file = paste(results_dir, "Year_", MY_YEAR, var_name, "_IDWgrid.Rds", sep = ""))
  saveRDS(preds_list_stations, file = paste(results_dir_test, "Year_", MY_YEAR, var_name, "_IDWTest.Rds", sep = ""))
}



#test.pred.mu
####################################################################################
#function to inverse distance weighting to interpolate to test grid
#also, extract the predictions at station locations to validate method
####################################################################################
IDW_interpolation = function(climate_data_list, var_of_interest = "Tref", testgrid_LCC_DF,
                             VAR, Lapse_rate) {
  YEAR = climate_data_list$YEAR
  MONTH = climate_data_list$MONTH
  DAY = climate_data_list$DAY
  
  #station data w/ reference transformed variable Tref
  #this is a spatial object w/ LCC coords
  mydata_station = climate_data_list$mydata_station
  
  #same as above but non-spatial -- just a data.frame
  mydata_station.lccDF = climate_data_list$mydata_station.lccDF
  
  #1km prediction grid--raster object
  lcc_test_grid_raster = climate_data_list$lcc_test_grid_raster

  #want to use mydata_station, which is a spatial object containing reference variable values to be interpolated 
  #to grid. Tref is still what we want to interpolate -- will back transform once predictions are made
  myformula = formula(paste(var_of_interest,"~ 1"))
  
  gs = gstat(formula = myformula, locations = mydata_station, nmax=9, set=list(idp = 2))
  
  #interpolate makes a raster using the specifications from the gstat function
  nn = interpolate(lcc_test_grid_raster, gs)
  
  #get interpolated data from nn raster
  nnCoords = data.frame(coordinates(nn))
  interp_data = as.data.frame(nn)
  interp_data = cbind(interp_data, nnCoords)
  names(interp_data) = c("var1.pred", "LON", "LAT")
  
  #this is very important. You must arrange in order otherwise elevation from testgridlatlon will 
  #not match interp_data b/c the coordinates will not be sorted correctly!
  #testgrid_LCC_DF has already been sorted by LON, LAT
  # head(testgrid_LCC_DF)
  head(interp_data)
  interp_data = arrange(interp_data, LON, LAT)
  ELEV = testgrid_LCC_DF$ELEV
  interp_data = cbind(interp_data, ELEV)
  
  #back transform to actual elevation value
  interp_data = back_transform(interp_data, VAR, Lapse_rate)
  
  
  YEAR = rep(YEAR, nrow(interp_data))
  DAY = rep(DAY, nrow(interp_data))
  MONTH = rep(MONTH, nrow(interp_data))
  interp_data = cbind(YEAR, MONTH, DAY, interp_data)
  
  #now backtransform from LCC to lat/lon -- lcc coords are still included in this data frame (regular datframe)
  #the variable you will want to bias correct later is VAR_interp_corr
  interp_data_latlon = transform_latlon(interp_data)
  
  #returns a data.frame of orignial station data w/ a new column of predictions at station locations
  #called TMAX_interp_corr (or TMIN... or PRECIP...)
  #mydata_station is a spatial object w/ lcc coord system
  station_preds = make_station_predictions(nn, mydata_station, mydata_station.lccDF,
                                           lcc_test_grid_raster, VAR, Lapse_rate)
  
  #interp_data_latlon is the prediction grid w/ interpolated values in lat/lon
  #station_preds is a dataframe w/ original data and interpolated values so we can validate approach
  results_dat = list(interp_data_latlon = interp_data_latlon, station_preds = station_preds)
  
  return(results_dat)
}

###################################################################################
#import_test_grid -- #import 1km prediction grid (derived from D3 DEM).***LCC TRANSFORMED!!**
###################################################################################
#import 1km prediction grid (derived from D3 DEM).#LCC TRANSFORMED VALUES
import_test_grid = function(){
  testgrid = read.csv("/raid/mholthui/fromBabbage/data/LCC_grid1km/testGrid1km_RACC_LCC.csv")
  testgrid$X = NULL
  return(testgrid)
}

import_test_grid_LCC_DF = function(){
  testgrid = read.csv("/raid/mholthui/fromBabbage/data/LCC_grid1km/testGrid1km_RACC_LCC.csv")
  testgrid$X = NULL
  testgrid = arrange(testgrid, LON, LAT)
  return(testgrid)
}

###################################################################################
#import daily station data /raid/mholthui/fromBabbage/data/WRF_data/
#######################################################################################
#import bias corrected data
import_daily_station_data <- function(VAR){
  if (VAR == "tmax"){  
    df = readRDS("/raid/mholthui/fromBabbage/data/WRF_data/Bias_corrected_RCP85tmax_1976_2099_IDW.Rds")
  }
  if (VAR == "precip"){
    df = readRDS("/raid/mholthui/fromBabbage/data/WRF_data/Bias_corrected_RCP45precip_1976_2099_IDW.Rds")
  }
  
  if (VAR == "tmin"){
    df = readRDS("/raid/mholthui/fromBabbage/data/WRF_data/Bias_corrected_RCP85tmin_1976_2099_IDW.Rds")
  }
  return(df)
}
#make progress log
make_progress_log = function(results_dir){
  prog_log = paste(results_dir[1], "prog.txt", sep = "")
  return(prog_log)
}

#create directories to put results in 
#grid predictions go into ...IDWgrid. Predictions at station locations go IDWtest
create_results_dir = function(){
  #directory to dump results into
  if(!dir.exists("/raid/mholthui/fromBabbage/scenario016_IDWgridPRCP_RCP45")){
    dir.create("/raid/mholthui/fromBabbage/scenario016_IDWgridPRCP_RCP45")
  }
  
  if(!dir.exists("/raid/mholthui/fromBabbage/scenario016_IDWTestPRCP_RCP45")){
    dir.create("/raid/mholthui/fromBabbage/scenario016_IDWTestPRCP_RCP45")
  }
  return(c("/raid/mholthui/fromBabbage/scenario016_IDWgridPRCP_RCP45/",
           "/raid/mholthui/fromBabbage/scenario016_IDWTestPRCP_RCP45/"))
}



################################################################
##preprocess_climate_data_train
## just apply baseline elevation transformation to the training dataset
################################################################ 
#Note: the variable will be called 'T_corr' as a result from doing the bias correction
preprocess_climate_data_train = function(mydata_station, lcc_test_grid_raster, VAR = "precip", Lapse_rate){
  YEAR = mydata_station$YEAR[1]
  MONTH = mydata_station$MONTH[1]
  DAY = mydata_station$DAY[1]
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
  
  if (VAR == "tmax"){
    
    mydata_station$Tref = mydata_station$T_corr - Lapse_rate*(200 - mydata_station$ELEV)
  }
  
  if (VAR == "tmin"){
    mydata_station$Tref = mydata_station$T_corr - Lapse_rate*(200 - mydata_station$ELEV)
  }
  
  if (VAR == "precip"){
    mydata_station$Tref = mydata_station$T_corr *((1 + Lapse_rate*(200 - mydata_station$ELEV))/ 
                                          (1 - Lapse_rate*(200 - mydata_station$ELEV)))
  }
  
  mydata_station.lccDF = mydata_station
  
  coordinates(mydata_station) <- c("LON_LCC", "LAT_LCC")
  proj4string(mydata_station) <- CRS.new

  
  return(list(mydata_station = mydata_station, 
              lcc_test_grid_raster = lcc_test_grid_raster, 
              mydata_station.lccDF = mydata_station.lccDF,
              YEAR = YEAR, MONTH = MONTH, DAY = DAY))
}


#make prediction raster w/ LCC transformation applied
preprocess_climate_data_test = function(test_set){
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
  coordinates(test_set) = ~ LON + LAT
  proj4string(test_set) = CRS.new
  
  #set downscaling grid to be truly gridded
  gridded(test_set) = TRUE
  #transform into a raster for IDW
  r = raster(test_set)
  return(r)
}


transform_latlon = function(interp_data){
  
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
  
  CRS.LatLon = CRS("+init=epsg:4326")
  
  dfsp = interp_data
  LCC_coords = dfsp[, c("LON", "LAT")]
  names(LCC_coords) = c("LON_LCC", "LAT_LCC")
  #dfsp = dfsp[, c("LON", "LAT", "TMAX_interp_corr")]
  coordinates(dfsp) = ~LON + LAT
  proj4string(dfsp) = CRS.new
  
  interp_data_latlon = spTransform(dfsp, CRS.LatLon)
  coordsLatLon = data.frame(coordinates(interp_data_latlon))
  interp_info = interp_data_latlon@data
  interp_data_latlon = cbind(interp_info, coordsLatLon, LCC_coords) 
  
  return(interp_data_latlon)
  
}


reference_transform = function(mydata_station, VAR, Lapse_rate){

    if (VAR == "tmax"){
      
      mydata_station$Tref = mydata_station$TMAX - Lapse_rate*(200 - mydata_station$ELEV)
    }
    
    if (VAR == "tmin"){
      mydata_station$Tref = mydata_station$TMIN - Lapse_rate*(200 - mydata_station$ELEV)
    }
    
    if (VAR == "precip"){
      mydata_station$Tref = mydata_station$PRCP *((1 + Lapse_rate*(200 - mydata_station$ELEV))/ 
                                            (1 - Lapse_rate*(200 - mydata_station$ELEV)))
    }
    
    return(mydata_station)
  }



#function to back transform to actual value based on elevation
#apply back transformation to get the actual elevation-adjusted value for VAR
back_transform = function(interp_data, VAR, Lapse_rate){
  if (VAR == "tmax"){
    interp_data$TMAX_interp_corr = interp_data$var1.pred - Lapse_rate*(interp_data$ELEV - 200)
    interp_data$var1.pred = NULL
    interp_data$test.pred.mu_grid = interp_data$TMAX_interp_corr
  }
  
  if (VAR == "tmin"){
    interp_data$TMIN_interp_corr = interp_data$var1.pred - Lapse_rate*(interp_data$ELEV - 200)
    interp_data$var1.pred = NULL
    interp_data$test.pred.mu_grid = interp_data$TMIN_interp_corr
  }
  
  if (VAR == "precip"){
    interp_data$PRECIP_interp_corr = interp_data$var1.pred *((1+Lapse_rate*(interp_data$ELEV - 200)) /
                                                               (1-Lapse_rate*(interp_data$ELEV - 200)))
    interp_data$var1.pred = NULL
    interp_data$test.pred.mu_grid = interp_data$PRECIP_interp_corr
  }
  
  
  return(interp_data)
}

#testspdf is the spatiopointsdataframe version of mydata_station
#get predictions at station locations so we can validate the approach.
#mydata_station is a spatial pts dataframe w/ lcc transform
#mydata_station.lccDF is the nonspatial version of the same thing
make_station_predictions = function(nn, mydata_station, mydata_station.lccDF,
                                    lcc_test_grid_raster, VAR, Lapse_rate){
  
  #lcc transformation
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
  
  #grab interpolated values -- these are still at baseline elevation
  #station values located in testspdf
  if (VAR == "tmax"){
    Ref_TMAX_interp = extract(nn, mydata_station)
    
    #add them to the original data
    mydata_station.lccDF$Ref_TMAX_interp = Ref_TMAX_interp
    
    #get rid of NA values
    mydata_station.lccDF = mydata_station.lccDF[complete.cases(mydata_station.lccDF$Ref_TMAX_interp), ]
  }
  
  if (VAR == "tmin"){
    Ref_TMIN_interp = extract(nn, mydata_station)
    
    #add them to the original data
    mydata_station.lccDF$Ref_TMIN_interp = Ref_TMIN_interp
    
    #get rid of NA values
    mydata_station.lccDF = mydata_station.lccDF[complete.cases(mydata_station.lccDF$Ref_TMIN_interp), ]
  }
  
  if (VAR == "precip"){
    Ref_PRECIP_interp = extract(nn, mydata_station)
    
    #add them to the original data
    mydata_station.lccDF$Ref_PRECIP_interp = Ref_PRECIP_interp
    
    #get rid of NA values
    mydata_station.lccDF = mydata_station.lccDF[complete.cases(mydata_station.lccDF$Ref_PRECIP_interp), ]
  }
  
  mydata_stationSP = mydata_station.lccDF
  #make into spatialPointsdataframe--again, this time without the NA values
  coordinates(mydata_stationSP)  = ~ LON_LCC +LAT_LCC
  #transform to lat/lon
  proj4string(mydata_stationSP) =  CRS.new
  
  #grab elevation values from the test grid raster and append to non-spatial dataframe
  DEM_ELEV = extract(lcc_test_grid_raster, mydata_stationSP)
  mydata_station.lccDF$DEM_ELEV = DEM_ELEV
  
  #back transform to actual variable value based on elevation from test grid
  if (VAR == "tmax"){
    mydata_station.lccDF$TMAX_interp_corr = mydata_station.lccDF$Ref_TMAX_interp - 
      Lapse_rate*(mydata_station.lccDF$DEM_ELEV - 200)
    mydata_station.lccDF$Ref_TMAX_interp = NULL
  }
  
  if (VAR == "tmin"){
    mydata_station.lccDF$TMIN_interp_corr = mydata_station.lccDF$Ref_TMIN_interp - 
      Lapse_rate*(mydata_station.lccDF$DEM_ELEV - 200)
    mydata_station.lccDF$Ref_TMIN_interp = NULL
  }
  
  
  if (VAR == "precip"){
    mydata_station.lccDF$PRECIP_interp_corr = mydata_station.lccDF$Ref_PRECIP_interp *((1 + Lapse_rate*(mydata_station.lccDF$DEM_ELEV - 200)) / (1 - Lapse_rate*(mydata_station.lccDF$DEM_ELEV)))
    mydata_station.lccDF$Ref_PRECIP_interp = NULL
  }
  
  return(mydata_station.lccDF)
}
