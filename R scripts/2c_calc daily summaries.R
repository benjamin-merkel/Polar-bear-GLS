#library(devtools)
#install_github("SLisovski/twGeos")

library(GeoLight)
library(RColorBrewer)
library(TwGeos)
library(stringr)

load("data/GLS files.RData")
load("data/DEG files.RData")
load("data/Full GLS meta data.RData")
load("output/LYR lufthavn temperature data 2001-2021.RData")

meta2$Capture.GLS.off.lat <- as.numeric(meta2$Capture.GLS.off.lat)
meta2$Capture.GLS.off.lon <- as.numeric(meta2$Capture.GLS.off.lon)

path <- "data/"



t3 = NULL
for(m in 1:nrow(meta2)){
  # load light data
  
  cat("\r",m)
  
  meta3 <- meta2[m,]
  meta3$year.off <- strftime(meta3$GLS_off_date,"%Y")
  
  file_select <- deg.ids$file[deg.ids$gls %in% meta3$GLS_ID]
  file_select <- file_select[grepl(meta3$Bear_ID, file_select)]
  file_select <- file_select[grepl('driftadj', file_select)]
  if(length(file_select)==0) file_select <- deg.ids$file[deg.ids$gls %in% meta3$GLS_ID]
  file_select <- file_select[grepl(meta3$year.off, file_select)]
  file_select <- file_select[1]
  
  t1          <- read.table(paste0(path,file_select[1]), skip=19, header = T)
  t1$datetime <- as.POSIXct(strptime(paste(t1[,1], t1[,2]), "%d/%m/%Y %H:%M:%S"), tz = "UTC")
  t1$date     <- as.Date(t1$datetime)
  t1          <- t1[as.Date(t1$datetime) > as.Date(meta3$GLS_on_date) & 
                      as.Date(t1$datetime) < as.Date(meta3$GLS_off_date),]
  t1$temp     <- t1[,3]
  t1$temp.var <- apply(t1[,5:8],1,var)
  t1$burst.diff <- abs(t1[,5]-t1[,6])+abs(t1[,7]-t1[,6])+abs(t1[,7]-t1[,8])
  t1$temp.diff<- t1$Tmax..C. - t1$Tmin..C.
  
  t1$wet      <- t1[,9]
  t1          <- t1[order(t1$date),]
  if(nrow(t1)>0){
    l2 <- data.frame(mean.temp =  tapply(t1$temp, t1$date, mean),
                     min.temp  =  tapply(t1$temp, t1$date, min),
                     max.temp  =  tapply(t1$temp, t1$date, max),
                     mean.temp.var =  tapply(t1$temp.var, t1$date, mean),
                     min.temp.var  =  tapply(t1$temp.var, t1$date, min),
                     max.temp.var  =  tapply(t1$temp.var, t1$date, max),
                     min.temp.diff= tapply(t1$temp.diff, t1$date, min),
                     max.temp.diff= tapply(t1$temp.diff, t1$date, max),
                     mean.temp.diff= tapply(t1$temp.diff, t1$date, mean),
                     min.burst.diff= tapply(t1$burst.diff, t1$date, min),
                     max.burst.diff= tapply(t1$burst.diff, t1$date, max),
                     mean.burst.diff= tapply(t1$burst.diff, t1$date, mean),
                     sum.wet   =  tapply(t1$wet,  t1$date, sum),
                     date      =  unique(t1$date), 
                     bear.gls  =  meta2$bear.gls[m])
    if(is.null(t3)) t3 <- l2 else t3 <- rbind(t3,l2)
  }
}

t3$doy      <- as.numeric(strftime(t3$date,"%j"))
t3$doy2     <- t3$doy - 181
t3$doy2[t3$doy2 < 1] <- t3$doy2[t3$doy2 < 1] +366
t3$year     <- as.numeric(strftime(t3$date,"%Y"))
t3$year2    <- t3$year
t3$year2[t3$doy < 182]    <- t3$year[t3$doy < 182] -1
t3 <- t3[order(t3$year2),]
t3 <- t3[order(t3$doy2),]
t3$year.id <- paste(t3$year2,t3$bear.gls)
t3$animal_id <- str_split_fixed(t3$bear.gls, " ", 2)[,1]

save(t3,file="output/bear temp and act data summarized daily.RData")



l3 = NULL
for(m in 1:nrow(meta2)){
  # load light data
  
  cat("\r",m)
  
  meta3 <- meta2[m,]
  meta3$year.off <- strftime(meta3$GLS_off_date,"%Y")
  
  file_select <- gls.ids$file[gls.ids$gls %in% meta3$GLS_ID]
  file_select <- file_select[grepl(meta3$Bear_ID, file_select)]
  file_select <- file_select[grepl('driftadj', file_select)]
  if(length(file_select)==0) file_select <- gls.ids$file[gls.ids$gls %in% meta3$GLS_ID]
  file_select <- file_select[grepl(meta3$year.off, file_select)]
  file_select <- file_select[1]
  
  gls          <- read.table(paste0(path,file_select[1]), skip=19, header = T)
  gls$datetime <- as.POSIXct(strptime(paste(gls[,1], gls[,2]), "%d/%m/%Y %H:%M:%S"), tz = "UTC")
  gls$date     <- as.Date(gls$datetime)
  gls          <- gls[as.Date(gls$datetime) > as.Date(meta3$GLS_on_date) & 
                        as.Date(gls$datetime) < as.Date(meta3$GLS_off_date),]
  gls$light    <- gls[,3]
  gls          <- gls[order(gls$date),]
  if(nrow(gls)>0){
    l2 <- data.frame(mean.light=  tapply(gls$light, gls$date, mean),
                     min.light =  tapply(gls$light, gls$date, min),
                     max.light =  tapply(gls$light, gls$date, max),
                     date      =  unique(gls$date), 
                     bear.gls  =  meta2$bear.gls[m])
    if(is.null(l3)) l3 <- l2 else l3 <- rbind(l3,l2)
  }
}


l3$doy      <- as.numeric(strftime(l3$date,"%j"))
l3$doy2     <- l3$doy - 181
l3$doy2[l3$doy2 < 1] <- l3$doy2[l3$doy2 < 1] +366
l3$year     <- as.numeric(strftime(l3$date,"%Y"))
l3$year2    <- l3$year
l3$year2[l3$doy < 182]    <- l3$year[l3$doy < 182] -1
l3 <- l3[order(l3$year2),]
l3 <- l3[order(l3$doy2),]
l3$year.id <- paste(l3$year2,l3$bear.gls)
l3$animal_id <- str_split_fixed(l3$bear.gls, " ", 2)[,1]

save(l3,file="output/bear light data summarized daily.RData")
