# ====================================================================================================================================================================
#
# the maximum of this vertical profile is selected between 40°N and 80°N (Polar sector)
# it needs to source the script "Meandering.lat_ref.clean.R" 
# it returns the maximum values of the meandering index, the latitude and the isohypse at which this value is found
#
# ====================================================================================================================================================================

# choose whether to plot hgt field
plot_TF = F
# please notice:
# insert missing values from line 22-> 63, 137-141

# =================================================================================
# set library path (if needed)
.libPaths( c( .libPaths(), "/home/your_lib/") )

# =================================================================================
# load libraries 
library(ncdf4)
library("rworldmap")
# =================================================================================
# load functions
source("/path-to-your-folder/MI_LatRef_functions.R") # loads function "MI_lat_ref" and all dendencies 

# =================================================================================
# set in and out paths
output_path_data <- "/output_path_data/Data/"
output_path_plots <- "/output_path_plots/Plots/"
data_name <- "label-your-experiment"

# fill in
start_year <- 1979 # start year
end_year <- 2016 # end year
years <- seq(start_year, end_year, 1)

lat_north <- 89 # northern lat taken into account, 90 must be excluded
lat_south <- 0 # souther lat take into account , the 2 hemispheres must be kept separated
# =================================================================================
# load geopotential height fields (can be daily or running mean fields)
# load all longitudes and latitudes
# needs del29feb!!!
nc <- nc_open(paste("/path-to-your-data/netcdf-hgt500-data.nc", sep=""))


# =================================================================================
# grid must start at -178.5 or similar and end with 180, translate if needed
# remeber to apply translation to the hgt fields
# extract grid info (change)
dim_lon_old <- ncvar_get(nc, "lon") 
dim_lon <- dim_lon_old - 178.5 # translate to a grid which starts from -180 and goes up to 180
dim_lat_old <- ncvar_get(nc, "lat")
dim_lat <- dim_lat_old[dim_lat_old>=lat_south & dim_lat_old<lat_north] # all lats are extracted and only the chosen lat_south and lat_north are taken into account

# =================================================================================
# extract gph field
hgt3D <- ncvar_get(nc, "var129") # [lon,lat,lev,time]
hgt3D <- hgt3D/9.80665 # needed for ERA-I 
hgt3D_NH <- hgt3D[c(122:240, 0:121),dim_lat_old>=lat_south & dim_lat_old<lat_north,] # select only northern hemisphere and traslate to -178.5 to 180

n_days <- 365
n_years <- dim(hgt3D_NH)[3]/n_days

# from 3D to 4D
hgt4D <- array(0, c(length(dim_lon), length(dim_lat), n_years, n_days))

for(y in 0:(n_years -1)){
  hgt4D[,,y+1, ] <- hgt3D_NH[,,(1+n_days*(y)):(n_days+n_days*(y))]
}

# =================================================================================
# check that the dimension are loaded correctly
# example
if((min(dim_lon)!= -178.5)|(min(dim_lat)!=0)|(max(dim_lon)!=180)|(max(dim_lat)!= 88.5)){
  stop("dim_lon and dim_lat extremes contain errors")
}else{
  print("dimension correctly loaded")
}

# =================================================================================
# always plot and check the fields 
if(plot_TF == T){
  #x11()
  
  # open pdf
  fig_path = paste0(output_path_plots,"hgt500_1981_78.pdf")
  pdf(fig_path, width = 18, height = 12)
  
  # plot first day of first year
  field = hgt4D[,,3,78]
  
  min_scale = min(field) 
  max_scale = max(field) 
  
  nlon = length(dim_lon)
  nlat = length(dim_lat)
  # this conversion expects HadGHCN type ordering of the grid
  #dx   <- 360/nlon
  #dy   <- 180/(nlat-1) # points at SP and NP
  dx = (max(dim_lon)-min(dim_lon))/(nlon-1)
  dy = (max(dim_lat)-min(dim_lat))/(nlat-1)
  
  set_Polypath(FALSE)
  
  # make grid for mapping
  offset = c(min(dim_lon),min(dim_lat))
  cellsize = c(dx,dy) 
  cells.dim=c(nlon,nlat)
  # create lon-lat grid with same dimension as input data
  gt<-GridTopology(cellcentre.offset=offset,cellsize=cellsize,cells.dim=cells.dim)

  gridVals<-data.frame(att=as.vector(field))
  
  # dataframe: grid + data
  sGDF<-SpatialGridayataFrame(gt,data=gridVals)
  
  # define color palette
  #colourPalette= c('lightseagreen','palegreen','skyblue1', 'slateblue1','slateblue4', 'burlywood1', 'burlywood4', 'sienna1', 'peachpuff')
  colourPalette = c("blue", "white", "indianred3")
  catMethod=seq(from=min_scale,to=max_scale,by=(max_scale-min_scale)/41) # ranges for label
  catMethod_labels=seq(from=min_scale,to=max_scale,by=round((max_scale-min_scale)/10)) # labels every 2 colors

  mapParams<-mapGridayedayata(sGDF,nameColumnToPlot='att',catMethod=catMethod,colourPalette= colourPalette, #"rainbow",
                            adayLegend=FALSE,#xlim=c(-180,180),ylim=c(0,45),
                            lwd = 0.5,oceanCol="white",landCol="white",borderCol="black")
  
  do.call(adayMapLegend,c(mapParams,legendLabels="all",legendWidth=1,legendMar=5, #legmar 11.5
                         legendShrink=1,labelFontSize=1, tcl=0., digits = 2, horizontal=F)) # legShr
  title("HGT500 1st jan 1979")  #dev.off()
  dev.off()
}

# =================================================================================
source("/home/dicapua/MeanderingIndex/MI_Github/test/MI_LatRef_functions.R") # loads function "MI_lat_ref" and all dendencies 

# =================================================================================
#core of the function
# a hgt4D field in the form [longitude, latidude, year, day ] is needed

# HEIGHT
isolevel <- seq(4900,6200,5) # vertical profile of isohypses

longitudes <- dim_lon
latitudes <- rev(dim_lat) # countourLines needs increasing lats, check whether needed
# =================================================================================
# loop over the number of days 
# =================================================================================

# parameters to be chosen 
ref_lat <- 60 # the default is 60°N
verbose <- 0 # T or F
lat_min <- 50 # default
lat_max <- 75 # defualt

# init vertical profile isolevel curviness and latitude 
MI_vertical_profile <- array(0, c(n_years, n_days, length(isolevel) ))
Orig_lat_vertical_profile <- array(0, c(n_years, n_days, length(isolevel) ))
isolevel_profile <- array(0, c(n_years, n_days, length(isolevel) ))
max_curviness <- array(0, c(n_years, n_days))
max_isohypse <- array(0, c(n_years, n_days))
max_latitude <- array(0, c(n_years, n_days))
# loop over number of years
for(yy in 1){
#for(yy in 1:n_years){
  year = years[yy]
  print(year)
    # loop over number of days
  #for(day in 1:n_days){
  for(day in 1:50){
      
    # choose daily and yearly hgt field 
    hgt_field <- hgt4D[,rev(1:length(dim_lat)), yy, day] # lats need to be increasing, check whether this is needed
    # initialization of vectors
    MI_daily_profile <- array(0, length(isolevel)) # vertical daily profile of curviness
    Lat_daily_profile <- array(0, length(isolevel))  # vertical daily profile of mean latitudes
    isolvl_daily_profile <- array(0, length(isolevel))  # vertical daily profile of isohypse level
    
    # =================================================================================
    # run MI_lat_ref function from MI_LatRef_functions.R
    # =================================================================================
    for(i in 1:length(isolevel)){
      isolvl  = isolevel[i] #choose the height 
      #print(paste("ampl.curv ","**********isolevel************",isolvl,"********",sep=" "))
      output <- MI_lat_ref(longitudes, latitudes, hgt_field, isolvl, ref_lat, verbose, year, day)
      MI_daily_profile[i] <- output[[1]]
      Lat_daily_profile[i] <- output[[2]]
      isolvl_daily_profile[i] <- isolvl
      #print(MI_daily_profile[i])
    }#end for
    
    # store all data
    MI_vertical_profile[yy, day, ] <- MI_daily_profile
    Orig_lat_vertical_profile[yy, day, ] <- Lat_daily_profile
    isolevel_profile[yy, day, ] <- isolvl_daily_profile
    
    # =================================================================================
    # maximum of the profile in the interval 50°N - 75°N (mean original latitude of the isopleth)
    # =================================================================================
    
    
    arraylevel <- which((Lat_daily_profile > lat_min) & (Lat_daily_profile < lat_max))
    lat_profile<- Lat_daily_profile[arraylevel]
    curv_profile <- MI_daily_profile[arraylevel]
    level_profile <- isolevel[arraylevel]
    
    if(length(curv_profile)==0){
      max_curviness[yy, day] <- 0
      max_isohypse[yy, day] <- 0
      max_latitude[yy, day] <- 0
    }else if(length(max(curv_profile))!=1 ){
      print(max(curv_profile))
      max_curviness[yy, day] <- 0
      max_isohypse[yy, day] <- 0
      max_latitude[yy, day] <- 0
    }else if(length(max(curv_profile))==1){
      #print(lat_profile)
      maximum <- which(curv_profile==max(curv_profile)) # find position of the maximum
      #print(max(curv_profile))
      if(length(maximum)>1){
        max_curviness[yy, day] <- 0
        max_isohypse[yy, day] <- 0
        max_latitude[yy, day] <- 0
      }else{
        if(curv_profile[maximum]!=max(curv_profile)){
          print("the maximum is not equal to the value chosen as maximum")
        }
        max_curviness[yy, day] <- curv_profile[maximum] # maximum of the index
        max_isohypse[yy, day] <- level_profile[maximum] # level of the maximum index
        max_latitude[yy, day] <- lat_profile[maximum] # latitude of the maximum index
      }
      
    } # end if we want to plot the most meandering isopeth but we need actually to catch the whole isopleth (which could be cut at the borders)
    
    # create output list max_isohypse
   
  }# loop over days
  
}
if(verbose>=2){
  print(paste("isohypse ",max_isohypse,sep=" "))
  print(paste("max_curviness ",max_curviness, sep="" ))
  print(paste("latitude",max_latitude, sep=" "))
}
# create list to save all data
MI_max_data <- list(years= years, n_days=n_days, max_curviness=max_curviness, max_isohypse= max_isohypse, max_latitude=max_latitude )


# save data (optional)
print('Save data in ')
print(paste(output_path_data,"MI_MI_max_data", data_name,".RData", sep=""))
print(paste(output_path_data,"MI_MI_vertical_profile", data_name,".RData", sep=""))
print(paste(output_path_data,"MI_Lat_daily_profile", data_name,".RData", sep=""))
#save(MI_max_data, file=paste(output_path_data,"MI_MI_max_data", data_name,".RData", sep=""))
#save(MI_vertical_profile, file=paste(output_path_data,"MI_MI_vertical_profile", data_name,".RData", sep=""))
#save(Orig_lat_vertical_profile, file=paste(output_path_data,"Orig_lat_vertical_profile", data_name,".RData", sep=""))
#save(isolevel_profile, file=paste(output_path_data,"isolevel_profile", data_name,".RData", sep=""))















