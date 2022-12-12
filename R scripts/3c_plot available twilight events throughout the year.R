library(ggplot2)
library(GeoLight)
library(RColorBrewer)
library(TwGeos)




light.lon <- 15.6333; light.lat <- 78.2166
light.lon <- 15; light.lat <- 76
dates <- as.Date(0:365, "2016-01-01")
ll <- data.frame(date = dates,
                 rise=twilight(tm = as.POSIXct(paste(dates,"12:00:00"),tz="UTC"), 
                               lon = light.lon, lat = light.lat, rise = T, zenith = 93.55),
                 set=twilight(tm = as.POSIXct(paste(dates,"12:00:00"),tz="UTC"), 
                              lon = light.lon, lat = light.lat, rise = F, zenith = 93.55))
ll$length <- abs(as.numeric(difftime(ll$rise,ll$set,units="hours")))
ll$doy <- as.numeric(strftime(ll$date,"%j"))
ll$length[is.na(ll$length) & (ll$doy < 265 & ll$doy > 81)] <- 24
ll$length[is.na(ll$length) & (ll$doy > 265 | ll$doy < 81)] <- 0
ll$length[ll$length < 12 & (ll$doy < 265 & ll$doy > 81)] <- 24 - ll$length[ll$length < 12 & (ll$doy < 265 & ll$doy > 81)] 


load("data/Full GLS location dataset.RData")
dd <- data.frame(table(factor(path2$doy, levels = c(1:365))))
colnames(dd) <- c("doy","freq")


m <- data.frame(date=dates)
m$month <- as.numeric(strftime(m$date,"%m"))
m$doy      <- as.numeric(strftime(m$date,"%j"))
m <- m[!duplicated(m$month),]
m <- rbind(m,data.frame(date=as.Date("2017-01-01"),month=1,doy=366))




season.col   <- brewer.pal(4, "Set2")[c(2,3)]
no.tw.colour <- '#FFD580'
no.tw.colour <- brewer.pal(8, "Set2")[6]
alpha.value  <- 0.9

Sys.setlocale("LC_ALL", "English")
png(paste0("figures/available twilight events 2.png"), res = 800, width=19, height = 16, units="cm")
opar <- par(mar=c(2,4,2,4))
plot(as.numeric(dd$doy), dd$freq, type="l", xaxt="n", col=grey(1), lty=3,ylab='',axes=F,yaxs="i", xaxs="i", ylim=c(0, 350),xlim=c(1,366))
y <- c(dd$freq,0,0)
x <- c(dd$doy,1,365)

polygon(x,y,col=grey(0.5),border=grey(0.5))

rect(0, 0, 30, 350 ,col=alpha(no.tw.colour,1), density = 10 ,border="transparent", lwd=2)
rect(316, 0, 370, 350 ,col=alpha(no.tw.colour,1), density = 10 ,border="transparent", lwd=2)
rect(min(ll$doy[ll$length==24]), 0, max(ll$doy[ll$length==24]), 350 ,col=alpha(no.tw.colour,1), density = 10 ,border="transparent", lwd=2)

fall.doy <- quantile(path2$doy[is.na(path2$lat) & path2$doy > 182], c(0.05,0.95),na.rm=T)
rect(fall.doy[1], 0, fall.doy[2], 350, border="transparent", angle = -45, density=10, col=alpha(season.col[1],alpha.value), lwd=3)

spring.doy <- quantile(path2$doy[is.na(path2$lat) & path2$doy < 182], c(0.05,0.95),na.rm=T)
rect(spring.doy[1], 0, spring.doy[2], 350, border="transparent", angle = -45, density=10, col=alpha(season.col[2],alpha.value), lwd=3)

lines(ll$doy,ll$length/24*350,col=no.tw.colour,lwd=3.5)

abline(h=0)
axis(1,at=m$doy,labels = strftime(m$date,"%b"),cex.axis=1)
axis(2, las=2)
axis(4, at = c(0, 350/4, 350/2, 350/2+350/4, 350), labels = seq(0,24,6), las=2,
     col=colorspace::darken(no.tw.colour, 0.1), col.axis=colorspace::darken(no.tw.colour, 0.1))

stage.ticks <- c(-10, 30, spring.doy, range(ll$doy[ll$length==24]), fall.doy, 316, 380)
stage.labs <- stage.tick.labs <- c(1, 30, spring.doy, range(ll$doy[ll$length==24]), fall.doy, 316, 365)
for(i in 1:(length(stage.tick.labs) - 1)) stage.tick.labs[i]  <- mean(stage.tick.labs[i:(i+1)])
stage.tick.labs <- stage.tick.labs[1:(length(stage.tick.labs)-1)]
for(i in 1:(length(stage.labs) - 1)) stage.labs[i]  <- stage.labs[(i+1)] - stage.labs[i] + 1
stage.labs <- stage.labs[1:(length(stage.labs)-1)]

axis(3, at = stage.ticks, labels = F)
# axis(3, at = stage.tick.labs, labels = stage.labs, tick = F, col.axis = c("orange", 1, 4, 1, "orange", 1, 2, 1, "orange"))
axis(3, at = stage.tick.labs[c(2,4,6,8)], labels = stage.labs[c(2,4,6,8)], tick = F, col.axis = 1)
axis(3, at = stage.tick.labs[3], labels = stage.labs[3], tick = F, col.axis = season.col[2])
axis(3, at = stage.tick.labs[7], labels = stage.labs[7], tick = F, col.axis = season.col[1])
mtext("Frequency of twilight events recorded", 2, line =3)
mtext("Day length [h]", 4, line =2.5, col = colorspace::darken(no.tw.colour, 0.1))

par(opar)
dev.off()

