library(GeoLight)
library(TwGeos)
library(SGAT)
library(lubridate)
library(MASS)
source("R scripts/FUNCTIONS.r")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### calibrate solar angle -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


load("data/Full twilight times dataset GLS.RData")

# assess noisy twilight seasons (due to dirty loggers?)

bgc <- table(twl6b$bgc,twl6b$Deleted)
bgc2<- bgc[,1]/(bgc[,1]+bgc[,2])


png("figures/Proportion of seasons with clean or noisy data.png",units = "cm",width=13,height=13,res = 500)
opar <- par(mar=c(4,4,0.8,0.8))
hist(1-bgc2, xlim=c(0,1), ylim=c(0,100), freq = T, las = 1, xaxs="i", yaxs="i", 
     col=grey(0.2), border = grey(0.2), 
     main="" , xlab="Proportion of twilight events labeled as likely erroneous")
abline(v=0.3, col = 2, lty = 2, lwd = 2)
par(opar)
dev.off()


# only consider seasons with relatively clean data to calibrate solar angle
twl6c <- twl6b[twl6b$bgc %in% names(bgc2[bgc2>0.7]),]
twl6c <- twl6c[!is.na(twl6c$Twilight) & twl6c$Deleted==F,]

# only consider twilight events up to 0.5 years away from a capture event to 
# minimize potential error by using inhabitat calibration with capture positions
twl6d <- twl6c[twl6c$deploy_timediff < 182,]
twl6d$cap.lon <- twl6d$Capture.GLS.on.lon
twl6d$cap.lat <- twl6d$Capture.GLS.on.lat
twl6e <- twl6c[twl6c$retrieve_timediff < 182,]
twl6e$cap.lon <- twl6e$Capture.GLS.off.lon
twl6e$cap.lat <- twl6e$Capture.GLS.off.lat
twl6d <- rbind(twl6d, twl6e)

# remove months characterized by polar night or midnight sun as well as equinoxes
days <- c(1:65,95:250,278:366)
twl6f <- twl6d[twl6d$doy %in% days & twl6d$month %in% c(2:4,8:11),]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### calibrate solar angle based on capture locations ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


### calibrate based on all locations relative to capture location
png("figures/solar calibration/2022 solar angle calib based on all inds and capture locations.png",units = "cm",width=15,height=15,res = 500)
cap.calib <- MythresholdCalibration(twilight = twl6f$Twilight, 
                                    rise     = twl6f$Rise, 
                                    lon      = twl6f$cap.lon, 
                                    lat      = twl6f$cap.lat, 
                                    method   = "gamma")
dev.off()

# only fall
twl6f <- twl6d[twl6d$doy %in% days & twl6d$month %in% c(8:11),]
png("figures/solar calibration/2022 solar angle calib based on all inds and capture locations fall.png",units = "cm",width=15,height=15,res = 500)
cap.calibf <- MythresholdCalibration(twilight = twl6f$Twilight, 
                                     rise     = twl6f$Rise, 
                                     lon      = twl6f$cap.lon, 
                                     lat      = twl6f$cap.lat, 
                                     method   = "gamma")
dev.off()

# only spring
twl6f <- twl6d[twl6d$doy %in% days & twl6d$month %in% c(2:4),]
png("figures/solar calibration/2022 solar angle calib based on all inds and capture locations spring.png",units = "cm",width=15,height=15,res = 500)
cap.calibs <- MythresholdCalibration(twilight = twl6f$Twilight, 
                                     rise     = twl6f$Rise, 
                                     lon      = twl6f$cap.lon, 
                                     lat      = twl6f$cap.lat, 
                                     method   = "gamma")
dev.off()


# calibration for each individual separately
twl6f <- twl6d[twl6d$doy %in% days & twl6d$month %in% c(2:4,8:11),]
for(ii in 1:length(unique(twl6f$bear.gls))){
  calib1 <- NULL
  bg <- unique(twl6f$bear.gls)[ii]
  tt7 <- twl6f[twl6f$bear.gls == bg,]
  png(paste0("figures/solar calibration/ind capture location/CAP_",bg,".png"),units = "cm",width=15,height=15,res = 500)
  try(calib1 <- MythresholdCalibration(twilight = tt7$Twilight, 
                                       rise     = tt7$Rise, 
                                       lon      = tt7$cap.lon, 
                                       lat      = tt7$cap.lat, 
                                       method   = "gamma"),silent = T)
  dev.off()
  if(!is.null(calib1)) {
    cal1 <- data.frame(calib1[1],calib1[2],calib1[3],calib1[4],calib1[5],calib1[6],calib1[7],bg)
    names(cal1) <- c(names(calib1),"bear.gls")
    cal1$months <- paste(sort(unique(tt7$month)),collapse = "-")
    
    if(ii==1) cal <- cal1 else cal <- rbind(cal,cal1)
  }
}
cap.cal <- cal[cal$a1>90,]
cap.cal$effective.n <- cap.cal$n - cap.cal$n.NA
cap.cal$type="CAPTURE"
cap.cal
summary(cap.cal)
write.table(cap.cal, file = "data/solar angle calibration with capture locations.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### calibrate solar angle based on capture locations ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

twl6f <- twl6b[twl6b$doy %in% days & twl6b$month %in% c(2:4,8:11) & 
                 !is.na(twl6b$gps.lon) & twl6b$Deleted==F &
                 twl6b$bgc %in% names(bgc2[bgc2>0.7]),]


png("figures/solar calibration/2022 solar angle calib based daily gps locations.png",units = "cm",width=15,height=15,res = 500)
gps.calib <- MythresholdCalibration(twilight = twl6f$Twilight, 
                                    rise     = twl6f$Rise, 
                                    lon      = twl6f$gps.lon, 
                                    lat      = twl6f$gps.lat, 
                                    method   = "gamma")
dev.off()

twl6f <- twl6c[twl6c$doy %in% days & twl6c$month %in% c(8:11) & !is.na(twl6c$gps.lon),]
png("figures/solar calibration/2022 solar angle calib based daily gps locations fall.png",units = "cm",width=15,height=15,res = 500)
gps.calibf <- MythresholdCalibration(twilight = twl6f$Twilight, 
                                     rise     = twl6f$Rise, 
                                     lon      = twl6f$gps.lon, 
                                     lat      = twl6f$gps.lat, 
                                     method   = "gamma")
dev.off()

twl6f <- twl6c[twl6c$doy %in% days & twl6c$month %in% c(2:4) & !is.na(twl6c$gps.lon),]
png("figures/solar calibration/2022 solar angle calib based daily gps locations spring.png",units = "cm",width=15,height=15,res = 500)
gps.calibs <- MythresholdCalibration(twilight = twl6f$Twilight, 
                                     rise     = twl6f$Rise, 
                                     lon      = twl6f$gps.lon, 
                                     lat      = twl6f$gps.lat, 
                                     method   = "gamma")
dev.off()



# calibration for each individual separately
twl6f <- twl6c[twl6c$doy %in% days & twl6c$month %in% c(2:4,8:11) & !is.na(twl6c$gps.lon),]
# 3, 24, 27 throw errors
for(ii in 1:length(unique(twl6f$bear.gls))){ #c(1,2,4:23,25,26)
  calib1 <- NULL
  bg <- unique(twl6f$bear.gls)[ii]
  tt7 <- twl6f[twl6f$bear.gls==bg,]
  png(paste0("figures/solar calibration/ind gps location/GPS_",bg,".png"),units = "cm",width=15,height=15,res = 500)
  try(calib1 <- MythresholdCalibration(twilight = tt7$Twilight, 
                                       rise     = tt7$Rise, 
                                       lon      = tt7$Capture.GLS.on.lon, 
                                       lat      = tt7$Capture.GLS.on.lat, 
                                       method   = "gamma"), silent = T)
  dev.off()
  if(!is.null(calib1)) {
    
    cal1 <- data.frame(calib1[1],calib1[2],calib1[3],calib1[4],calib1[5],calib1[6],calib1[7],bg)
    names(cal1) <- c(names(calib1),"bear.gls")
    cal1$months <- paste(sort(unique(tt7$month)),collapse = "-")
    
    if(ii==1) cal <- cal1 else cal <- rbind(cal,cal1)
  }
}
gps.cal <- cal
gps.cal$effective.n <- gps.cal$n - gps.cal$n.NA
gps.cal$type="GPS"
gps.cal
summary(gps.cal)
write.table(gps.cal, file = "data/solar angle calibration with gps locations.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### combine calibration -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


cal <- merge(cap.cal[,c('a1',"bear.gls","effective.n")], gps.cal[,c('a1',"bear.gls","effective.n")], by= 'bear.gls',all.y = T)
cal <- cal[!is.na(cal$a1.x),]
cal$a1.diff <- cal$a1.x-cal$a1.y
write.table(cal, file = "data/solar angle calibration with catpure and gps locations.txt")

plot((90-cap.cal$a1), (cap.cal$effective.n), xlim=c(-2.5,-4.5),pch=19,
     ylab="# of twilight sets used to calibrate", xlab="solar calibration angle")
points((90-gps.cal$a1), (gps.cal$effective.n),pch=19,col=2)
for(i in 1:nrow(cal)) lines(c(90-cal$a1.x[i],90-cal$a1.y[i]), c(cal$effective.n.x[i],cal$effective.n.y[i]),col=grey(0.5))
points((90-cap.cal$a1), (cap.cal$effective.n),pch=19,col=1)
abline(v=90-median(cap.cal$a1),lty=2)
abline(v=90-median(gps.cal$a1),lty=2,col=2)

cal2 <- rbind(gps.cal, cap.cal)
png("figures/solar calibration/2022 Ind solar angle calibration distribution based on caputre and gps locations.png",res=500,units="cm",width=15,height=20)
boxplot(90-cal2$a1~cal2$type,col=c("white"),border=c(1,2),lty=1,range=0,ylab="solar angle", las=1,xlab="")
for(i in 1:nrow(cal)) lines(c(1,2),c(90-cal$a1.x[i],90-cal$a1.y[i]),col=grey(0.7))
points(rep(1,nrow(cap.cal)),90-cap.cal$a1,pch=19)
points(rep(2,nrow(gps.cal)),90-gps.cal$a1,pch=19,col=2)
axis(at=c(1,2),side=3,labels = table(cal2$type))
abline(h = 90-cap.calib[1], lty=2)
abline(h = 90-gps.calib[1], lty=2, col=2)
dev.off()


