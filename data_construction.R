#library(spNNGP)
#library(spBayes)
library(sp)
#library(ggpubr)
library(dplyr)
library(data.table)
library(reshape)

##############################################################################################
##############################################################################################
# UPDATED: May 23, 2020
# get_wrf_data
# function to create wrf data by year. Returns 1 year of daily wrf data for either 
#historical or future projections (use historical = FALSE for future projections)
# Accounts for leap years! Day Dec. 31 will be removed for leap years FOR HISTORICAL DATA ONLY
# for future , rcp 85 driven data, there are no leap days 
##############################################################################################
# parameters: 
# MY_YEAR: integer in 1980-2014
# VAR: 'tmax', 'tmin' or 'precip'
# historical: 'TRUE' or 'FALSE' for historical or future projections, respectively
# Value: a data.frame with 1 years of daily data for WRF
###############################################################################################
###############################################################################################
get_wrf_data = function(MY_YEAR, VAR, historical = TRUE, rcp = "RCP45"){
  require(ncdf4)
  if (!(VAR %in% c("tmax","tmin","precip" ))){
    stop("Variable name must be one of tmax, tmin, or precip")
  }
  if (MY_YEAR < 1976 || MY_YEAR > 2099){
    stop("MY_YEAR must be between 1976 and 2099")
  }
  
  if (VAR == "precip"){
    VAR_NAME = "PRCP"
    VAR = "RAIN"
  } 
  
  if (VAR == "tmax"){
    VAR = "TMAX"
    VAR_NAME = "TMAX"
  }
  
  if (VAR == "tmin"){
    VAR = "TMIN"
    VAR_NAME = "TMIN"
  }
  
  if (historical == T){
    type = "Historical simulations"
  }else{
    type = "Future projections"
  }
  
  
  print(paste("processing data for ", type, " for variable ", VAR, " ", MY_YEAR))
  print("Please note: there are no leap days for RCP85/RCP45-driven WRF projections and simulations")
  
  if (!(historical == TRUE)){
    #open the netcdf 
       #open the netcdf 
    if (rcp == "RCP45"){ #/raid/mholthui/fromBabbage/data/WRF_RCP45/
    print("using RCP4.5 data")
    mync = nc_open(paste("/raid/mholthui/fromBabbage/data/WRF_RCP45/cmip5.bruyere-bc.ccsm4.rcp45.4km.wrf.", MY_YEAR, ".nc", sep = ""))
    }
    if (rcp == "RCP85") {
    print("using RCP8.5 data")
    mync = nc_open(paste("/raid/mholthui/fromBabbage/data/WRF_RCP85/cmip5.bruyere-bc.ccsm4.rcp85.4km.wrf.",MY_YEAR, ".nc", sep = ""))
    }
    
    #get lat and lon
    lat <- ncvar_get(mync, varid="LAT")
    LAT=as.vector(lat)
    lon <- ncvar_get(mync, varid="LON")
    LON=as.vector(lon)
    elev = ncvar_get(mync, varid = "HGT")
    ELEV = as.vector(elev)
    
    months_one_year = c(rep(1, 134757), rep(2, 121716), rep(3, 134757), rep(4, 130410), rep(5, 134757), rep(6, 130410),
                        rep(7, 134757), rep(8, 134757), rep(9, 130410), rep(10, 134757), rep(11, 130410), rep(12, 134757))
    
    #label days correctly (NOTE THERE IS NO FEB 29!!!!!!!!!)
    jan_days = rep(1:31, each = 4347)
    feb_days = rep(1:28, each = 4347)
    march_days = rep(1:31, each = 4347)
    april_days = rep(1:30, each = 4347)
    may_days = rep(1:31, each = 4347)
    june_days = rep(1:30, each = 4347)
    july_days = rep(1:31, each = 4347)
    aug_days = rep(1:31, each = 4347)
    sept_days = rep(1:30, each = 4347)
    oct_days = rep(1:31, each = 4347)
    nov_days = rep(1:30, each = 4347)
    dec_days = rep(1:31, each = 4347)
    
    #id for grid points
    ID = rep(1:4347, 365)
    
    
    #YEAR identifier
    YEAR = rep(MY_YEAR, 1586655)
    
    #Day id
    DAY  =  c(jan_days, feb_days, march_days, april_days, may_days, june_days, 
              july_days, aug_days, sept_days, oct_days, nov_days, dec_days)
    
    
    
    ELEMENT =  as.vector(ncvar_get(mync, varid = VAR, start=c(1, 1, 1), count=c(-1, -1, 365)))
    
    daily_data_df = as.data.frame(cbind(LAT, LON, ELEMENT, ELEV))
    
    
    daily_data_df$YEAR = YEAR
    daily_data_df$MONTH = months_one_year
    daily_data_df$DAY = DAY
    daily_data_df$ID = ID
    daily_data_df = daily_data_df[, c(5,6,7,3,8,1,2,4)]
    colnames(daily_data_df)[4] = VAR_NAME
    
    return(daily_data_df)
  }
  
  else{
    
    if (MY_YEAR < 1980 || MY_YEAR > 2014){
      stop("MY_YEAR must be between 1980 and 2014 for historical data!")
    }
    
    #open the netcdf file
    mync = nc_open(paste("/epscorfs/data/BREE_Climate_Projections/erai/erai.4km.wrf/daily/erai.4km.wrf.",
                         MY_YEAR, ".nc", sep = ""))
    #get lat and lon
    lat <- ncvar_get(mync, varid="LAT")
    LAT=as.vector(lat)
    lon <- ncvar_get(mync, varid="LON")
    LON=as.vector(lon)
    elev = ncvar_get(mync, varid = "HGT")
    ELEV = as.vector(elev)
    
    # in the case of a leap year, keep Feb 29 and get rid of December 31
    if (MY_YEAR %in% seq(1980, 2014, by = 4)){
      
      
      
      #Months and days for one year for leap year
      
      months_one_year = c(rep(1, 134757), rep(2, 126063), rep(3, 134757), rep(4, 130410), rep(5, 134757), rep(6, 130410),
                          rep(7, 134757), rep(8, 134757), rep(9, 130410), rep(10, 134757), rep(11, 130410), rep(12, 130410))
      
      #label days correctly (NOTE THERE IS NO DEC 31!!!!!!!!!)
      jan_days = rep(1:31, each = 4347)
      feb_days = rep(1:29, each = 4347)
      march_days = rep(1:31, each = 4347)
      april_days = rep(1:30, each = 4347)
      may_days = rep(1:31, each = 4347)
      june_days = rep(1:30, each = 4347)
      july_days = rep(1:31, each = 4347)
      aug_days = rep(1:31, each = 4347)
      sept_days = rep(1:30, each = 4347)
      oct_days = rep(1:31, each = 4347)
      nov_days = rep(1:30, each = 4347)
      dec_days = rep(1:30, each = 4347)
      
      #id for grid points
      ID = rep(1:4347, 365)
      
      
      #YEAR identifier
      YEAR = rep(MY_YEAR, 1586655)
      
      #Day id
      DAY  =  c(jan_days, feb_days, march_days, april_days, may_days, june_days, 
                july_days, aug_days, sept_days, oct_days, nov_days, dec_days)
      
      
      print("leap year! Please note, Dec 31 will be removed")
      ELEMENT =  as.vector(ncvar_get(mync, varid=VAR, start=c(1,1, 1), count=c(-1,-1,365)))
      #print(num_days[i])
      daily_data_df = as.data.frame(cbind(LAT, LON, ELEMENT, ELEV))
      
      
      daily_data_df$YEAR = YEAR
      daily_data_df$MONTH = months_one_year
      daily_data_df$DAY = DAY
      daily_data_df$ID = ID
      daily_data_df = daily_data_df[, c(5,6,7,3,8,1,2,4)]
      colnames(daily_data_df)[4] = VAR_NAME
      head(daily_data_df)
      
    }else{
      
      
      #Months and days for one year
      print("not a leap year!")
      months_one_year = c(rep(1, 134757), rep(2, 121716), rep(3, 134757), rep(4, 130410), rep(5, 134757), rep(6, 130410),
                          rep(7, 134757), rep(8, 134757), rep(9, 130410), rep(10, 134757), rep(11, 130410), rep(12, 134757))
      
      #label days correctly (NOTE THERE IS NO FEB 29!!!!!!!!!)
      jan_days = rep(1:31, each = 4347)
      feb_days = rep(1:28, each = 4347)
      march_days = rep(1:31, each = 4347)
      april_days = rep(1:30, each = 4347)
      may_days = rep(1:31, each = 4347)
      june_days = rep(1:30, each = 4347)
      july_days = rep(1:31, each = 4347)
      aug_days = rep(1:31, each = 4347)
      sept_days = rep(1:30, each = 4347)
      oct_days = rep(1:31, each = 4347)
      nov_days = rep(1:30, each = 4347)
      dec_days = rep(1:31, each = 4347)
      
      #id for grid points
      ID = rep(1:4347, 365)
      
      
      #YEAR identifier
      YEAR = rep(MY_YEAR, 1586655)
      
      #Day id
      DAY  =  c(jan_days, feb_days, march_days, april_days, may_days, june_days, 
                july_days, aug_days, sept_days, oct_days, nov_days, dec_days)
      
      
      
      ELEMENT =  as.vector(ncvar_get(mync, varid = VAR, start=c(1, 1, 1), count=c(-1, -1, 365)))
      
      daily_data_df = as.data.frame(cbind(LAT, LON, ELEMENT, ELEV))
      
      
      daily_data_df$YEAR = YEAR
      daily_data_df$MONTH = months_one_year
      daily_data_df$DAY = DAY
      daily_data_df$ID = ID
      daily_data_df = daily_data_df[, c(5,6,7,3,8,1,2,4)]
      colnames(daily_data_df)[4] = VAR_NAME
      
      
    }
    
    
    return(daily_data_df)
    
  }
  
  
}

###################################################################################################
#get_station_data()
###################################################################################################
# function to make station data. Returns 1 year of daily station data

# Parameters:
# MY_YEAR: integer in 1980-2014
# VAR: 'tmax', 'tmin', or 'precip'
# df = "/data/mholthui/StationData/daily_data2.RData" --> overall D2 station data -- you need this!
###################################################################################################

get_station_data = function(MY_YEAR, VAR, df){
  

  d2long  = subset(df, YEAR == MY_YEAR )
  d2long$DAY = as.numeric(d2long$DAY)
  
  
  if (VAR == "precip"){
    VAR = "PRCP"
  } 
  
  if (VAR == "tmax"){
    VAR = "TMAX"
  }
  
  if (VAR == "tmin"){
    VAR = "TMIN"
  }
  
  
  if (VAR == "tave"){
    setDT(d2long)
    d2wide = dcast(d2long, YEAR +MONTH + DAY+ID  +LAT+LON +ELEV~ ELEMENT, value.var = "VALUE")
    #head(d2wide)
    
    d2wide$TAVE = (d2wide$TMAX + d2wide$TMIN) / 2
    meas.vars= names(d2wide[, 8:ncol(d2wide)])
    
    d2long = melt(d2wide, id.vars = c("YEAR", "MONTH", "DAY", "ID", "LAT", "LON", "ELEV"), measure.vars = meas.vars)
    names(d2long) = c("YEAR", "MONTH", "DAY", "ID", "LAT", "LON", "ELEV", "ELEMENT", "VALUE")
    
    #if you do this again check below to make sure it matches the column order of orig. d2 data from when you import daily_data.Rdata
    d2long = d2long[ , c(1,2,8,3,9,4,5,6,7)]
  }
  
  d2long = filter(d2long, LAT >= 43.39 & LAT <= 45.861)
  d2long = filter(d2long, LON >= -74.758 & LON <= -71.556 )
  
  #take out na values
  d3long = d2long[complete.cases(d2long), ]
  d3long = filter(d3long, ELEMENT == VAR)
  d3long$ELEMENT = NULL
  colnames(d3long)[4] = VAR
  
  if (VAR == "TMAX"){
    bads1 = which(d3long[[VAR]] > 60)
    bads2 = which(d3long$MONTH == 8 & d3long[[VAR]] < 0)
    bads3 = which(d3long$MONTH == 7 & d3long[[VAR]] < 0)
    d3long = d3long[-c(bads1, bads2, bads3), ]
  }
  
  if (VAR == "TMIN"){
    bads1 = which(d3long[[VAR]] < -50)
    alldat_df = d3long[ -bads1, ]
  }
  
  if (VAR == "PRECIP"){
    bads1 = which(d3long[[VAR]] > 400)
    bads2 = which(d3long[[VAR]] < 0)
    d3long = d3long[-c(bads1, bads2), ]
    }
  
  return(d3long)
}
