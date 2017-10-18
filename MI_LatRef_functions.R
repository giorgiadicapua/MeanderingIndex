# this script calculate the meandering index for a single isohypse (contour of geopotential height) 
# the curviness value is calculated at the chosen reference latitude (60° N suggested for the Polar Front region)
# the script needs a daily field of geopotential height as input, latitudes and longitudes vectors, and the isolevel 
# and returns the value of the meandering index and its latitude



MI_lat_ref <- function ( longitudes, latitudes, hgt_field, isolvl, ref_lat, verbose ){
  
  if(verbose == TRUE){
    print("MI_lat_ref")
    print ("latitude of reference") 
    print(ref_lat)
    print("isolvl")
    print(isolvl)
    print("latitudes")
    print(latitudes)
    print("longitudes")
    print(longitudes)
  }

  #init values
  curv <- 0
  curv_lat <- 0
  
  # calculate position of the isopleth
  isopleth = contourLines(longitudes,latitudes,hgt_field,nlevels=1,levels=isolvl)
  
  # longest_iso returns 0 if the isopleth does not exist, 1 if only one isopleth is found or the number of the longest found isopleth 
  i_longest_iso <- longest_iso(isopleth)
  
  if(i_longest_iso == 0){
    # isopleth does not exist
    print("isopleth == 0")
    curv <- 0
    curv_lat <- 0 
    
  }else{
    # isopleth does exist
    longest_isop <- isopleth[[i_longest_iso]]
    # check whether the last and the first point of the isopleth[[i_longest_iso]]y are equal
    print(abs(longest_isop$y[1]-longest_isop$y[length(longest_isop$y)]))
    
    if(abs(longest_isop$y[1]-longest_isop$y[length(longest_isop$y)])>3){
      print("check whether the last and the first point of the isopleth[[i_longest_iso]]y are equal")
      print(paste("isop[1]",longest_isop$y[1],"isop[last]",longest_isop$y[length(longest_isop$y)] ))
      
      # how to solve the problem
      lon_translated <-list(lon1=c(61:240,1:60),lon2=c(121:240,1:120),lon3=c(181:240,1:180))
      for(ii in 1:3){
        print("inside the for")
        isopleth = contourLines(longitudes,latitudes,hgt_field[lon_translated[[ii]],],nlevels=1,levels=isolvl)
        i_longest_iso <- longest_iso(isopleth)
        longest_isop <- isopleth[[i_longest_iso]]
        
        if(abs(longest_isop$y[1]-longest_isop$y[length(longest_isop$y)])<3){
          output <- curviness.ref.lat.function(longitudes, latitudes,longest_isop, ref_lat)
          curv <- output[[1]]
          curv_lat <- output[[2]]
          if(curv!=0){ 
            contour(longitudes, latitudes,hgt_field[lon_translated[[ii]],], nlevel=10,levels=seq(5000,5900,100),main=curv)
            lines(longest_isop$x,longest_isop$y, col=189,lwd=3)
          } 
          break
        }else{
          curv <- 0
          curv_lat <- 0 
        }
      }# end loop over different latitudes
      
    }else{
      
      # calc isopleth curviness at the reference mean latitude and store the original mean latitude 
      # Isopleths that hit the latitude-borders are disgarded
      output <- curviness.ref.lat.function(longitudes, latitudes,longest_isop, ref_lat)
      curv <- output[[1]]
      curv_lat <- output[[2]]
    }# end if which check whether the isopleth is cut at the border of the field
    
  }
  output_list <- list(curv=curv, curv_lat=curv_lat)
  return(output_list)
}

# script for calculate curviness and curviness at a reference lat

curviness.ref.lat.function <- function(lon, lat,isop_old, ref_lat)
{
  #print(isop_old)
  #print(paste(min(isop_old$x),lon[1],max(isop_old$x),lon[length(lon)]))
  if ( 
    # isopleth should not move outside latitudinal range
    ((min(isop_old$y) > lat[2]) & (max(isop_old$y) < lat[length(lat)-1]))
    &
    # isopleth should circle whole globe  
    ((min(isop_old$x) == lon[1]) & (max(isop_old$x) == lon[length(lon)]))
   ){
    # first, the mean latitude at which the original isopleth is found must be calculated because we need to know the difference between the reference
    # latitude at which we want to calculate the isopleth
    print("calculating curviness at the latitude of reference")
    lat_or <- mean(isop_old$y)
    curv_lat <- lat_or
    # the new translated laditude must be defined (nothing happens to the longitudes)
    # pay attention: the reference latitude must be negative in the Southern hemisphere
    if(lat_or > ref_lat){
      lat_diff <- abs(lat_or-ref_lat)
      lat_new <- lat-lat_diff
    }else if(lat_or < ref_lat){
      lat_diff <- abs(lat_or-ref_lat)
      lat_new <- lat+lat_diff
    }
    # calculate isopleth curviness at the reference latitude
    # the isopleth must be recalculated
    isopleth_new = contourLines(lon,lat_new,hgt_field,nlevels=1,levels=isolvl)
    # the isopleth could not to exist: control
    i_longest_iso <- longest_iso(isopleth_new)
    isop_new <- isopleth_new[[i_longest_iso]]
    if(i_longest_iso == 0){
      curv <- 0
    }
    if ( 
      # isopleth should not move outside latitudinal range
      ((min(isop_new$y) > lat_new[2]) & (max(isop_new$y) < lat_new[length(lat_new)-1]))
      &
      # isopleth should circle whole globe  
      ((min(isop_new$x) == lon[1]) & (max(isop_new$x) == lon[length(lon)]))
    ){
      curv = calculate.isopleth.curviness(isop_new)
    }
  }else{
    print("isopleth is outside latitudinal range or doesn't circle the whole globe ")
    curv <- 0
    curv_lat <- 0 
    
  }
  
  output <- list(curv=curv, curv_lat=curv_lat)
  return(output)
  
}# end function


# calculate the isopleth curviness function
#
calculate.isopleth.curviness <- function ( isop )
{
  #print("calculate.isopleth.curviness")
  # determine x / y coordinates in [m] for this isopleth
  x_coor   = isop$x
  y_coor   = isop$y
  # check correct order (i.e. from west to east)
  if(x_coor[1] > x_coor[length(x_coor)]){
    # invert coordinate vectors
    x_coor = rev(x_coor)
    y_coor = rev(y_coor)
  }
  lat_mean = mean(y_coor)
  # shift to meter-based grid
  x_coor   = lon2meter( x_coor, lat_mean)
  y_coor   = lat2meter( y_coor, lat_mean)
  x_coor = duplicate.dateline.x (x_coor, lat_mean)
  y_coor = duplicate.dateline.y (y_coor)
  #plot(x_coor,y_coor)  
  # calculate distance in meters
  distance = isoline.distance(x_coor, y_coor)
  Earth.Radius  = 6367500
  cos_lat       = cos( lat_mean*pi/180. )
  circumference = 2*pi*Earth.Radius * cos_lat
  return( distance / circumference ) # definition of curviness
}

# finde the longest isopleth

longest_iso <- function(isopleth){
  
  if(length(isopleth)==0){
    
    i_longest_iso = 0
    
  }else{
    # get index of longest isopleth
    if(length(isopleth)>1){
      #print(length(isopleth))
      distance = 0
      for(iso in 1:length(isopleth)){
        if(length(isopleth[[iso]]$x) > distance){
          distance = length(isopleth[[iso]]$x)
          i_longest_iso = iso
        } 
      }
    } else {
      #print(length(isopleth))
      i_longest_iso = 1
    }
  }
  return(i_longest_iso)
}

isoline.distance <- function ( x, y )
{
  #print("isoline.distance")
  
  # calculates the total eucledian distance along a line given by coordinate vectors x and y
  # x and y coordinates should come in units of meters
  # check input:
  if( length(x) != length(y) ) {
    stop("Error: isoline.distance: coordinate vectors have unequal size!")
  }
  # loop over coordinates, start at i = 2
  total_distance = 0
  for(i in 2:length(x)){
    # pythagoras
    total_distance = total_distance + sqrt( (x[i]-x[i-1])^2 + (y[i]-y[i-1])^2 )
  }
  # done. return total distance
  return(total_distance)
}

lon2meter <- function ( lon, lat )
{
  
  #print("lon2meter")
  # returns distance in meters for input vector consisting of longitudinal coordinates "lon" at constant latitude "lat" (both in degrees)
  # check input:
  if( length(lat) != 1 ) { 
    stop("Error: lon2meter: single latitude expected!")
  }
  if( min(lon) != lon[1] ) { 
    stop("Error: lon2meter: method expects first coordinate to be the most westerly coordinate!")
  }
  Earth.Radius = 6367500
  cos_lat      = cos( lat*pi/180. )
  x0           = lon[1] * pi * Earth.Radius / 180. * cos_lat # most westerly point
  # most westerly point is set to x = 0 meter
  # thus, most easterly point should have coordinate (close to) circumpherence of Earth at this latitude
  return ( (lon-lon[1]) * pi * Earth.Radius / 180. * cos_lat);
}

lat2meter <- function ( lat, lat_zero )
{
  # print("lat2meter")
  
  # returns distance in meters of input vector "lat" to latitude "lat_zero" (both in degrees)
  # lat_zero is taken as zero-axis for new meter-based coordinate system
  # check input:
  if( length(lat_zero) != 1 ) { 
    stop("Error: lat2meter: single lat_zero latitude expected!")
  }
  Earth.Radius = 6367500
  return ( (lat-lat_zero) * pi * Earth.Radius / 180. );
}

duplicate.dateline.x <- function ( x, lat_mean )
{
  #print("duplicate.dateline.x")
  
  # a vector "x_" is returned with size length(x)+1
  # the last coordinate of "x" is duplicated and stored in x_[1]
  # check input:
  if( length(lat_mean) != 1 ) { 
    stop("Error: duplicate.dateline.x: single lat_mean latitude expected!")
  }
  eq_circ = 2*pi*6367500
  x_ = array(0,dim=c(length(x)+1))
  x_[1:length(x)] = x
  x_[length(x)+1] = cos(lat_mean*pi/180)*eq_circ
  return(x_)
}

duplicate.dateline.y <- function ( x )
{
  #print("duplicate.dateline.y")
  
  # a vector "x_" is returned with size length(x)+1
  # the first coordinate of "x" is duplicated and stored in last coordinate of x_
  x_ = array(0,dim=c(length(x)+1))
  x_[1:length(x)] = x
  x_[length(x)+1] = x[1]
  return(x_)
} 
