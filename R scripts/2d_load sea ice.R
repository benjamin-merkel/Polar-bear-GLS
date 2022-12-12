# https://github.com/CoastWatch-WestCoast/r_code/blob/master/accessing_projected_datasets.md


# Function to check if pkgs are installed, install missing pkgs, and load
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE,repos='http://cran.us.r-project.org')
    if(!require(x,character.only = TRUE)) stop(x, " :Package not found")
  }
}

list.of.packages <- c("ncdf4","openair","ggplot2","reshape2","scales","lubridate",
                      "cmocean","maps","mapdata","rgdal","raster","RColorBrewer", "sp")

# create list of installed packages
pkges = installed.packages()[,"Package"]
for (pk in list.of.packages) {
  pkgTest(pk)
}

library(sf)
library(terra)
library(move)

# Download the lat-lon GRID
# Use the dataset id 'nsidcCDRice_nh_grid' and the default axis extent values

url <- 'https://polarwatch.noaa.gov/erddap/griddap/'

grid_id <- 'nsidcCDRice_nh_grid'

grid_urlcall <- paste0(url,grid_id,'.nc?longitude[(5812500.0):1:(-5337500.0)][(-3837500.0):1:(3737500.0)],latitude[(5812500.0):1:(-5337500.0)][(-3837500.0):1:(3737500.0)]')

grid_nc <- download.file(grid_urlcall,destfile="data/grid.nc",mode='wb')

# Read the grid file
gridFid <- nc_open('data/grid.nc')
ygrid <- ncvar_get(gridFid, varid="ygrid")
xgrid <- ncvar_get(gridFid, varid="xgrid")
longitude <- ncvar_get(gridFid, varid="longitude")
latitude <- ncvar_get(gridFid, varid="latitude")
nc_close(gridFid)

inds = which(latitude > 70, arr.ind=TRUE)
rowrange <- range(inds[,1])
colrange <- range(inds[,2])

#Generate a DATA request URL using the indices from the grid

dataid <- 'nsidcG02202v4nhmday'
varnames <- c('cdr_seaice_conc_monthly','nsidc_bt_seaice_conc_monthly')
datestring <- '[(2012-01-01T00:00:00Z):1:(2021-12-01T00:00:00Z)]'
coordstring <- paste0('[',colrange[1]-1,':1:',colrange[2]-1,'][',rowrange[1]-1,':1:',rowrange[2]-1,']')

for (i in 1:length(varnames)) {
  if (i == 1) {
    urlcall <- paste0(url,dataid,'.nc?',varnames[i],datestring,coordstring)
  } 
  else {
    urlcall <- paste0(urlcall,',',varnames[i],datestring,coordstring)
  }
}

#Download the netCDF file  (this will take a few minutes, 20 years of data)
data_nc <- download.file(urlcall,destfile="data/NSIDC monthly sea ice conc 2012-2021.nc",mode='wb')

# Read the downloaded netCDF file and load the ice data into R variables

dataFid <- nc_open("data/NSIDC monthly sea ice conc 2012-2021.nc")

datatime <- ncvar_get(dataFid, varid="time")
datatime <- as.Date(as.POSIXlt(datatime,origin='1970-01-01',tz= "GMT"))

ygrid <- ncvar_get(dataFid, varid="ygrid")
xgrid <- ncvar_get(dataFid, varid="xgrid")

seaiceCDR <- ncvar_get(dataFid, varid=varnames[1])
seaiceGoddard <- ncvar_get(dataFid, varid=varnames[2])

nc_close(dataFid)

# Request a grid subset using the same coordinate string used for the data download. 
urlcall <- paste0(url,grid_id,'.nc?longitude',coordstring,',latitude',coordstring) 
grid_subset <- download.file(urlcall,destfile="data/grid_subset.nc",mode='wb')

# Read and format the subsetted grid data from the netCDF file  
gridSubsetFid <- nc_open('data/grid_subset.nc')

ygrid <- ncvar_get(gridSubsetFid, varid="ygrid")
xgrid <- ncvar_get(gridSubsetFid, varid="xgrid")
longitudeSubset <- ncvar_get(gridSubsetFid, varid="longitude")
latitudeSubset <- ncvar_get(gridSubsetFid, varid="latitude")

nc_close(gridSubsetFid)


plotdate <- datatime[as.numeric(strftime(datatime, "%m")) %in% c(2,3,4)]
s <- vector(mode = "list")
for(year in 2012:2021){
  idate   <- which((month(datatime) %in% month(plotdate)) & (year(datatime) %in% year))
  arr     <- seaiceCDR[,,idate]
  mat     <- apply(arr,c(1,2),max)
  
  dims    <- dim(xgrid)
  icemap2 <- expand.grid(xgrid=xgrid, ygrid=ygrid)
  icemap2$Seaice <- c(mat)
  icemap2$Seaice[icemap2$Seaice > 2] <- NA 
  s[[year-2011]] <- rasterFromXYZ(icemap2[,c("xgrid","ygrid","Seaice")])
  c <- rasterToContour(s[[year-2011]], levels = c(0.15))
  if(year==2012) spring_c<- c else spring_c <- rbind(spring_c, c)
}
s <- stack(s)
projection(s) <- 3411
spring_c$year <- 2012:2021

plotdate <- datatime[as.numeric(strftime(datatime, "%m")) %in% c(8,9,10)]
a <- vector(mode = "list")
for(year in 2012:2021){
  idate   <- which((month(datatime) %in% month(plotdate)) & (year(datatime) %in% year))
  arr     <- seaiceCDR[,,idate]
  mat     <- apply(arr,c(1,2),max)
  
  dims    <- dim(xgrid)
  icemap2 <- expand.grid(xgrid=xgrid, ygrid=ygrid)
  icemap2$Seaice <- c(mat)
  icemap2$Seaice[icemap2$Seaice > 2] <- NA 
  a[[year-2011]] <- rasterFromXYZ(icemap2[,c("xgrid","ygrid","Seaice")])
  c <- rasterToContour(a[[year-2011]], levels = c(0.15))
  if(year==2012) autumn_c<- c else autumn_c <- rbind(autumn_c, c)
}
a <- stack(a)
projection(a) <- 3411
autumn_c$year <- 2012:2021 


plot(mean(s))
plot(spring_c,add=T)

plot(mean(a))
plot(autumn_c,add=T)

autumn_c <- st_as_sf(autumn_c)
spring_c <- st_as_sf(spring_c)
st_crs(spring_c) <- st_crs(autumn_c) <- 3411

saveRDS(autumn_c, file = "data/Autumn sea ice concentration 15% contour 2012-2021.RDS")
saveRDS(spring_c, file = "data/Spring sea ice concentration 15% contour 2012-2021.RDS")
