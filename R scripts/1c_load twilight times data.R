library(GeoLight)
library(TwGeos)
library(SGAT)
library(lubridate)
source("R scripts/FUNCTIONS.r")

load("data/Full GLS meta data.RData")
load("data/GLS files.RData")
load("data/Daily gps location.RData")

# assess how many file were with or without clock drift adjustment.
gls.ids$selected <- 0
for(m in 1:nrow(meta2)){
  
  cat("\r",m)
  
  meta3 <- meta2[m,]
  meta3$year.off <- strftime(meta3$GLS_off_date,"%Y")
  
  file_select <- gls.ids$file[gls.ids$gls %in% meta3$GLS_ID]
  file_select <- file_select[grepl(meta3$Bear_ID, file_select)]
  file_select <- file_select[grepl('driftadj', file_select)]
  if(length(file_select)==0) file_select <- gls.ids$file[gls.ids$gls %in% meta3$GLS_ID]
  file_select <- file_select[grepl(meta3$year.off, file_select)]
  file_select <- file_select[1]
  
  gls.ids$selected[gls.ids$file == file_select] <- 1
}
table(gls.ids$driftadj, gls.ids$selected)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load twilight data ------

twl6 <- NULL
for(m in 1:nrow(meta2)){
  
  cat("\r",m, " of ",nrow(meta2),"  ")
  
  meta3 <- meta2[m,]
  meta3$year.off <- strftime(meta3$GLS_off_date,"%Y")
  
  file_select <- gls.ids$file[gls.ids$gls %in% meta3$GLS_ID]
  file_select <- file_select[grepl(meta3$Bear_ID, file_select)]
  file_select <- file_select[grepl('driftadj', file_select)]
  if(length(file_select)==0) file_select <- gls.ids$file[gls.ids$gls %in% meta3$GLS_ID]
  file_select <- file_select[grepl(meta3$year.off, file_select)]
  file_select <- file_select[1]
  
  gls          <- read.table(paste("data",file_select[1],sep="/"), skip=19, header = T)
  gls$datetime <- as.POSIXct(strptime(paste(gls[,1], gls[,2]), "%d/%m/%Y %H:%M:%S"), tz = "UTC")
  gls$date     <- as.Date(gls$datetime)
  gls          <- gls[as.Date(gls$datetime) > as.Date(meta3$GLS_on_date) & 
                        as.Date(gls$datetime) < as.Date(meta3$GLS_off_date),]
  
  twl <- twilightCalc(datetime = gls$datetime, light = log(gls$light.lux.), ask = F, maxLight = 10)
  twl$Rise    <- T
  twl$Rise[twl$type==2] <- F
  twl3b <- NULL
  
  if(nrow(twl)>6){
    
    # split into chunks divided by longer gaps of two weeks or more
    twl$diff  <- as.numeric(difftime(twl$tSecond, twl$tFirst, units = "hour"))
    twl$chunk <- 1
    for(i in 1:(nrow(twl)-1)) if(twl$diff[i]> 24 * 14) twl$chunk[(i+1):nrow(twl)] <- twl$chunk[(i+1):nrow(twl)] + 1
    twl <- twl[twl$chunk %in% names(table(twl$chunk))[table(twl$chunk)>14],]
    twl <- twl[twl$diff < 24 *14,]
    
    twl.switch <- twl
    twl.switch$tFirst <- twl.switch$tSecond
    twl.switch$Rise[twl.switch$type==2] <- T
    twl.switch$Rise[twl.switch$type==1] <- F
    twl2 <- rbind(twl, twl.switch)
    twl2 <- twl2[!duplicated(twl2[,1]),]
    twl2 <- twl2[,c(1,4,6)]
    colnames(twl2) <- c("Twilight","Rise","chunk")
    
    # remove nights or days less than 3 hours long
    twl2$diff  <- as.numeric(difftime(c(twl2$Twilight[2:nrow(twl2)],NA), twl2$Twilight, units = "hour"))
    twl2 <- twl2[!(twl2$diff<3 & twl2$Rise==F),]
    twl2 <- twl2[!is.na(twl2$Twilight),]
    
    offset <- -2 # adjusts the y-axis to put night (dark shades) in the middle
    # par(mar=c(4,4,1,1))
    for(chunks in unique(twl2$chunk)){
      try(twl3b <- MytwilightEdit(twilights = twl2[twl2$chunk==chunks,1:2],
                                  offset = offset,
                                  window = 4,           # two days before and two days after
                                  outlier.mins = 30,    # difference in mins
                                  stationary.mins = 15,  # are the other surrounding twilights within x mins of one another
                                  plot = F), silent = T)
      if(is.null(twl3b)) twl3b <- cbind(twl2[twl2$chunk==chunks,1:2], Deleted = F, Edited = F, Twilight0 = twl2[twl2$chunk==chunks,1])
      twl3b$chunk <- chunks
      if(chunks == unique(twl2$chunk)[1]) twl3 <- twl3b else twl3 <- rbind(twl3, twl3b)
    }
    
    
    twl3$year  <- as.numeric(strftime(twl3$Twilight,'%Y'))
    twl3$month <- as.numeric(strftime(twl3$Twilight,'%m'))
    
    twl4 <- twl3#[twl3$Deleted==F,]
    
    twl4$Twilight <- with_tz(twl4$Twilight, "UTC")
    twl4$date     <- as.Date(twl4$Twilight)
    
    twl4$bear.id <- meta3$Bear_ID
    twl4$gls.id  <- meta3$GLS_ID
    twl4$doy     <- as.numeric(strftime(twl4$Twilight, "%j"))
    
    twl4$deploy_timediff   <- difftime(twl4$Twilight, meta3$GLS_on_date,units = "days") 
    twl4$retrieve_timediff <- difftime(meta3$GLS_off_date, twl4$Twilight,units = "days") 
    
    
    if(is.null(twl6)) twl6 <- twl4 else twl6 <- rbind(twl6, twl4)
  }
}
twl6$bear.gls <- paste(twl6$bear.id, twl6$gls.id)
twl6          <- twl6[!is.na(twl6$year),]
twl6$bgy      <- paste(twl6$bear.id, twl6$gls.id,twl6$year)
metasub       <- meta2[,c("bear.gls","Capture.GLS.on.lon",'Capture.GLS.on.lat','Capture.GLS.off.lon','Capture.GLS.off.lat')]
twl6b         <- merge(twl6, metasub, by="bear.gls",all.x=T)
twl6b$Capture.GLS.on.lon  <- as.numeric(twl6b$Capture.GLS.on.lon)
twl6b$Capture.GLS.on.lat  <- as.numeric(twl6b$Capture.GLS.on.lat)
twl6b$Capture.GLS.off.lon <- as.numeric(twl6b$Capture.GLS.off.lon)
twl6b$Capture.GLS.off.lat <- as.numeric(twl6b$Capture.GLS.off.lat)
twl6b$bgc    <- paste(twl6b$bear.gls,twl6b$chunk)
twl6b$season <- "spring"
twl6b$season[twl6b$month > 7] <- "fall"
twl6b$bys       <- paste(twl6b$bear.id, twl6b$year, twl6b$season)
twl6b$bear.date <- paste(twl6b$bear.id, twl6b$date)
twl6b           <- merge(twl6b, gps, by="bear.date", all.x=T)

save(twl6b ,file="data/Full twilight times dataset GLS.RData")

