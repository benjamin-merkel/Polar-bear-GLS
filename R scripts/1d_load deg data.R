library(GeoLight)
library(TwGeos)
library(SGAT)
library(lubridate)


load("data/Full GLS meta data.RData")
load("data/GLS files.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load additional data ------

deg6 <- NULL
for(m in 1:nrow(meta2)){
  
  cat("\r",m)
  meta3 <- meta2[m,]
  meta3$year.off <- strftime(meta3$GLS_off_date,"%Y")
  
  file_select <- deg.ids$file[deg.ids$gls %in% meta3$GLS_ID]
  file_select <- file_select[grepl(meta3$GLS_ID, file_select) | grepl(meta3$Bear_ID, file_select)]
  file_select <- file_select[grepl(meta3$year.off, file_select)]
  file_select <- file_select[1]
  
  deg          <- read.table(paste("data",file_select[1],sep="/"), skip=19, header = T)
  deg$datetime <- as.POSIXct(strptime(paste(deg[,1], deg[,2]), "%d/%m/%Y %H:%M:%S"), tz = "UTC")
  deg$date     <- as.Date(deg$datetime)
  deg          <- deg[as.Date(deg$datetime) > as.Date(meta3$GLS_on_date) & 
                        as.Date(deg$datetime) < as.Date(meta3$GLS_off_date),]
  deg$year  <- as.numeric(strftime(deg$datetime,'%Y'))
  deg$month <- as.numeric(strftime(deg$datetime,'%m'))
  deg$doy   <- as.numeric(strftime(deg$datetime, "%j"))
  
  deg$bear.id <- meta3$Bear_ID
  deg$gls.id  <- meta3$GLS_ID
  if(is.null(deg6)) deg6 <- deg else deg6 <- rbind(deg6, deg)
}
save(deg6 ,file="data/Full deg dataset GLS.RData")
