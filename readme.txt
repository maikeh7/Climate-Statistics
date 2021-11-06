For use with paper: 
"Holthuijzen, M. F., Beckage, B., Clemins, P. J., Higdon, D., & Winter, J. M. (2021). Constructing High-Resolution, Bias-Corrected Climate Products: A Comparison of Methods. Journal of Applied Meteorology and Climatology, 60(4), 455-475."

This is a list of scenarios (ELEV/NO ELEV & METHOD) used in the paper. 


**scenario 1 =  train on wrf, predict to stations
scenario 1A = train on daily station data, predict to 1km grid, for use in alternative method of bias correction

scenario 2 = training on wrf OR station data and predict to 1km prediction grid

scenario 3 = train on wrf, predict to stations, THEN, predict to 1km grid (I'm not sure what I was thinking here)

**scenario 4 = Train on BIAS CORRECTED DATA, predict to 1km grid

scenario 5 = train on BIAS CORRECTED DATA, predict on 5 test sets. This is a REVISION of scenario 006, which is not really 
true cross-validation. You need to specify k. So, run all years w/ k=1, then w/ k=2, and so forth up to k=5

**scenario 6 = Train on BIAS CORRECTED DATA, predict on subset of station locations via 5 fold cross validation (this is to verify accuracy of downscaling) -- USING BAYESIAN KRIGING/MODEL FITTING

scenario 8 = Train on bias corrected data. But this time data were bias corrected using only years 1980-1986. So the 'future' years are 1987-1994. In this scenario, take the data corrected from 1980-86 and then do the same 

as scenario 5 (5 fold crossval over stations). Then we can see if doing bias correction for only a subset of the years causes a considerable decrease in perkins SS/RMSE. Data was bias corrected with function 

D:\school\Climate_Analysis_Personal\WRF_exploratory_data\bias_correction\bias_correction_functionCV.R. **Make sure to specify correct years for this, as not all years will be in the training set**

MainDriver.R = for all scenarios except 6
MainDriverCV.R = for scenario 6


MainDriverStrawMan = for scearios 009 and 010 and 013
*NEW*
Linear model
scenario 009: "Straw Man" model comparison to V1 downscaling. Using linear regression instead of Bayesian spatial model. May need to use Bayesian MLR....Results will be similar to those that come from scenario 001
scenario-010: Linear model analogue to scenario-004, except use linear model to predict to 1km grid

NonBayesian Kriging
scenario 007: 5 fold cross validation for NON bayesian kriging, analogue scenario 006 
scenario 011: train on WRF, predict to stations, using NON bayesian kriging, analogue to scenario 001
scenario 012: train on BIAS CORRECTED DATA (resulting of course from running scenario 011), predict to 1km grid, analogue to scenario 004

IDW
scenario 013: apply method from Winter et al. 2016 and output hi-res predictions and also predictions at station locations
Method is run using LCC coordinates for everything, because square grids are better for interpolation
scenario 014 Similar to scenario 001 except instead of using Bayesian spatial model , use IDW. 
scenario 014B: same as above but DO NOT USE ELEVATION
scenario 015: similar to scenario 6 (5 fold cross validation) except using IDW (by day, but see revision scenario 017 V2)**
scenario 016: similar to scenario 4, except using idw to predict to 1km grid
scenario 017: similar to scenario 5, except using idw 
scenario 017V2: similar to above, but predict directly to test station locations instead of predicting to the grid and 
extracting values. Probably more fair that way, and this is also what is done in scenario 005. USE THIS VERSION
scenario 018: IDW--WRF to 1km grid using elevation adjusted IDW
scenario018B: IDW--WRD to 1km grid NOT USING ELEVATION
scenario-014_IDW_to_station_FuturePRCP.R: specifically for FUTURE precip projections. For historical precip to station, use scenario-014.
*updated Nov 2020 to work with new WRF data obtained in June 2020

NOTE: LAT/LON after the preprocessing are the transformed aes coordinates, and LAT1/LON1 are the untransformed lat/lon coordinates
ALSO, I added LAT to the Kriging model in scenario 001. This has not been run yet, but you could run it. Otherwise, get rid of the 'LAT' in the formula.