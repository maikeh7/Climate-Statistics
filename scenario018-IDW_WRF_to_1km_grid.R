#scenario 18 - train on WRF , predict to 1km grid
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
#######################################################################################################
#train on WRF , and make predictions to 1km grid without bias correction. Also take into account elevation by 
#first bringing meteorological variable to reference elevation 
#(200m) and then backtransforming. 
Model_Driver = function(days=days, months=months, 
                        VAR, MY_YEAR=MY_YEAR, 
                        test_data_list=test_data_list,
                        df=NULL, prog_log = prog_log, 
                        lapse_params=lapse_params,
                        test_stations=NULL) {
  #this is 1km prediction raster grid
  lcc_test_grid_raster = test_data_list
  
  #this is just a data.frame, NOT a raster
  testgrid_LCC_DF = import_test_grid_LCC_DF()
  
  print(paste("running year ", MY_YEAR))
  
  #create WRF data for a year for VAR
  #"RCP85" for rcp85, "RCP45" for rcp45, use historical = TRUE ONLY for ERA data
  wrfdat80 = get_wrf_data(MY_YEAR, VAR, historical = FALSE, rcp = "RCP85")

  
  #list to hold results for predictions to 1km prediction grid
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
  
  
  climateDF_wrf = wrfdat80
  
  #loop over months
  for (k in 1:length(months)){
    cur_data_wrf = filter(climateDF_wrf, MONTH == months[k])
    
    #loop over days 
    for (i in 1:length(days)){
      var_of_interest = names(cur_data_wrf)[4]
      
      mydata_wrf = filter(cur_data_wrf, DAY == days[i])
      
      cat("month is: ", mydata_wrf$MONTH[1], "\n", 
          "day is: ", mydata_wrf$DAY[1],
          file = prog_log, append = TRUE)
      
      
      if (dim(mydata_wrf)[1] < 1){
        print("there are no observations! Skipping this day.")
        next
      }
      
      #preprocess wrf data:
      #WRF and  data will be returned as data.frames w/ LCC coords
      #transformation to reference elevation occurs in this function
      #a list of dataframes w/ lat/lon and LCC coords is returned (4 data.frames total)
      climate_data_list = preprocess_climate_data_train(mydata_wrf, VAR, Lapse_rate)
      
      #IDW interpolation to station locations; Tref is temp at reference elevation
      #Tref is created in preprocess_climate_data_train()
      #lcc_test_grid_raster is the 1km prediction grid in lcc coords
      results_dat = IDW_interpolation(climate_data_list, var_of_interest = "Tref",
                                      lcc_test_grid_raster, #raster to predict to (LCC)
                                      testgrid_LCC_DF, #same as raster except not a spatial data frame (LCC)
                                      VAR, 
                                      Lapse_rate)
      
      
      res = list(results_dat)
      
      # res will be a large list composed of ind. dataframes for each day. Dataframe 
      #will include the variable 'test.pred.mu' which are the interpolated values from IDW interpolation
      #to the 1km prediction grid
      preds_list = c(preds_list, res)
    }
  }
  results_dir = create_results_dir()
  saveRDS(preds_list, file = paste(results_dir, "/Year_", MY_YEAR, var_name, "IDW_grid.Rds", sep = ""))
}

####################################################################################
#function to inverse distance weighting to interpolate to test grid
####################################################################################
#
IDW_interpolation = function(climate_data_list, var_of_interest = "Tref",
                             lcc_test_grid_raster, testgrid_LCC_DF, VAR, Lapse_rate) {
  
  YEAR = climate_data_list$YEAR
  MONTH = climate_data_list$MONTH
  DAY = climate_data_list$DAY
  #this is WRF data w/ LCC coordinate system - spatial object
  mydata_wrf.lcc = climate_data_list$mydata_wrf.lcc
  
  #This is WRF data.frame w/ LCC coords
  mydata_wrf.lccDF = climate_data_list$mydata_wrf.lccDF
  
  #WRF data w/ lat long coords
  mydata_wrfDF_latlon = climate_data_list$mydata_wrfDF_latlon
  
  #want to use mydata_wrf.lccDF to make interpolated projections to lcc_test_grid_raster
  #Tref is the temp at reference elevation
  myformula = formula(paste(var_of_interest,"~ 1"))
  
  gs = gstat(formula = myformula, locations = mydata_wrf.lcc, nmax=9, set=list(idp = 2))
  
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
  # head(interp_data)
  interp_data = arrange(interp_data, LON, LAT)
  ELEV = testgrid_LCC_DF$ELEV
  interp_data = cbind(interp_data, ELEV)
  
  #back transform to actual elevation value
  interp_data = back_transform(interp_data, VAR, Lapse_rate)
  
  
  YEAR = rep(YEAR, nrow(interp_data))
  DAY = rep(DAY, nrow(interp_data))
  MONTH = rep(MONTH, nrow(interp_data))
  interp_data = cbind(YEAR, MONTH, DAY, interp_data)
  #names(interp_data)[4] = "test.pred.mu"
  results_dat = transform_latlon(interp_data)
  return(results_dat)
}

###################################################################################
#import_test_grid - import 1km prediction grid
###################################################################################
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
#import daily station data - don't need station data here
#######################################################################################
import_daily_station_data <- function(VAR = NULL){}
#   e <- new.env()
#   load("/data/mholthui/StationData/daily_data2.RData", envir = e)
#   df = get("df", envir =  e)
#   return(df)
# }

#make progress log
make_progress_log = function(results_dir){
  prog_log = paste(results_dir, "/prog.txt", sep = "")
  return(prog_log)
}

#create directory to put results in
create_results_dir = function(){
  
  if(!dir.exists("/raid/mholthui/fromBabbage/IDW_scenario018precip_RCP85_TOGRID")){
    dir.create("/raid/mholthui/fromBabbage/IDW_scenario018precip_RCP85_TOGRID")
  }
  return("/raid/mholthui/fromBabbage/IDW_scenario018precip_RCP85_TOGRID")
}



################################################################
##preprocess_climate_data_train
## transforms WRF data to lcc coordinates so we can interpolate
#returns original WRF/station dataframes, as well as dataframes w/ 
#LCC coordinates (note-these are not spatial objects), plus current 
#Month, Day, and Year
################################################################
preprocess_climate_data_train = function(mydata_wrf, VAR = "precip", Lapse_rate){
  YEAR = mydata_wrf$YEAR[1]
  MONTH = mydata_wrf$MONTH[1]
  DAY = mydata_wrf$DAY[1]
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
  
  #return spatial points dataframes of wrf data (lcc transformed)
  #plus original WRF datasets--denoted w/ DF (lat/lon dataframes)
  return(list(mydata_wrf.lccDF = mydata_wrf.lccDF, #data.frame w/ LCC coords
              mydata_wrfDF_latlon = mydata_wrfDF_latlon, #original dataframe w/ Lat lon coords
              mydata_wrf.lcc = mydata_wrf.lcc,
         YEAR = YEAR, MONTH = MONTH, DAY = DAY))#spatial object w/ lCC coords
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


#function to backtransform from reference elevation value to actual elevation value
back_transform = function(interp_data, VAR, Lapse_rate){
  if (VAR == "tmax"){
    interp_data$TMAX_interp_corr = interp_data$var1.pred - Lapse_rate*(interp_data$ELEV - 200)
    interp_data$var1.pred = NULL
  }
  
  if (VAR == "tmin"){
    interp_data$TMIN_interp_corr = interp_data$var1.pred - Lapse_rate*(interp_data$ELEV - 200)
    interp_data$var1.pred = NULL
  }
  
  if (VAR == "precip"){
    interp_data$PRECIP_interp_corr = interp_data$var1.pred *((1+Lapse_rate*(interp_data$ELEV - 200)) /
                                                               (1-Lapse_rate*(interp_data$ELEV - 200)))
    interp_data$var1.pred = NULL
  }
  
  return(interp_data)
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

#this just adds the lat lon coordinates back to the lcc prediction grid in case 
#you want lat lon coords as well
#lat lon coords are 'LON, LAT', LCC coords are 'LON_LCC', 'LAT_LCC'
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
