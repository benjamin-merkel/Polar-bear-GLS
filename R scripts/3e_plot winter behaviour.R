#library(devtools)
#install_github("SLisovski/twGeos")

library(GeoLight)
library(RColorBrewer)
library(TwGeos)
library(stringr)

# load("data/Full deg dataset GLS 2022.RData")
load("data/GLS files.RData")
load("data/DEG files.RData")
load("data/Full GLS meta data.RData")
load("output/LYR lufthavn temperature data 2001-2021.RData")

meta2$Capture.GLS.off.lat <- as.numeric(meta2$Capture.GLS.off.lat)
meta2$Capture.GLS.off.lon <- as.numeric(meta2$Capture.GLS.off.lon)

path <- "data/"
load("output/bear temp and act data summarized daily.RData")
load("output/bear light data summarized daily.RData")



# plot id by year timeline plots
ids <- unique(l3$animal_id)
for(id in ids){
  # id = 23732 # mistet all unge 
  #id = 23980 # fått 2 kull
  
  l4 <- l3[l3$animal_id == id,]
  l4 <- l4[order(l4$date),]
  t4 <- t3[t3$animal_id == id,]
  
  # plot(l4$date,l4$mean.light,type="l")
  years <- unique(l4$year2)
  
  for(year in years){
    
    m <- data.frame(date=as.Date(0:365, paste0(year,"-01-01")))
    m$month <- as.numeric(strftime(m$date,"%m"))
    m$doy      <- as.numeric(strftime(m$date,"%j"))
    m$doy2     <- m$doy - 182
    m$doy2[m$doy2 < 1] <- m$doy2[m$doy2 < 1] +366
    m <- m[!duplicated(m$month),]
    
    tlys <- l4[l4$year2 == year,]
    ll <- data.frame(date = tlys$date,
                     rise=twilight(tm = as.POSIXct(paste(tlys$date,"12:00:00"),tz="UTC"), 
                                   lon = meta2$Capture.GLS.off.lon[meta2$Bear_ID == tlys$animal_id[1]][1], 
                                   lat = meta2$Capture.GLS.off.lat[meta2$Bear_ID == tlys$animal_id[1]][1], 
                                   rise = T, zenith = 93.55),
                     set=twilight(tm = as.POSIXct(paste(tlys$date,"12:00:00"),tz="UTC"), 
                                  lon = meta2$Capture.GLS.off.lon[meta2$Bear_ID == tlys$animal_id[1]][1], 
                                  lat = meta2$Capture.GLS.off.lat[meta2$Bear_ID == tlys$animal_id[1]][1], 
                                  rise = F, zenith = 93.55))
    ll$length <- abs(as.numeric(difftime(ll$rise,ll$set,units="hours")))
    ll$doy <- as.numeric(strftime(ll$date,"%j"))
    ll$doy2     <- ll$doy - 182
    ll$doy2[ll$doy2 < 1] <- ll$doy2[ll$doy2 < 1] +366
    ll$length[is.na(ll$length) & (ll$doy < 265 & ll$doy > 81)] <- 24
    ll$length[is.na(ll$length) & (ll$doy > 265 | ll$doy < 81)] <- 0
    ll$length[ll$length < 12 & (ll$doy < 265 & ll$doy > 81)] <- 24 - ll$length[ll$length < 12 & (ll$doy < 265 & ll$doy > 81)] 
    # ll$length[ll$length > 12 & (ll$doy > 265 | ll$doy < 81)] <- 24 - ll$length[ll$length > 12 & (ll$doy > 265 | ll$doy < 81)]
    
    ll <- ll[order(ll$doy2),]
    
    png(paste('figures/ind timelines/',id,' ',year, '.png',sep=""), width = 15, height = 19, units="cm",res=600)
    
    par(mfrow=c(6,1),oma=c(2,0,0,0),mar=c(0.1,5,0.1,3))
    t4l <- l4[l4$year2 == year,]
    plot(t4l$doy2,t4l$max.light,type="l",xlim=c(40,304),
         xaxt="n",ylim=c(0,1),col=grey(0.8),lty=3,ylab="light intensity",yaxt="n",axes=F)
    y <- c(t4l$max.light,0,0)/max(t4l$max.light)
    x <- c(t4l$doy2,max(t4l$doy2),min(t4l$doy2))
    polygon(x,y,col=grey(0.8),border=grey(0.8))
    
    axis(2,at=c(0,1),labels = c(0,1),las=1)
    
    y <- c(t4l$mean.light,0,0)/max(t4l$max.light)
    x <- c(t4l$doy2,max(t4l$doy2),min(t4l$doy2))
    polygon(x,y,col=1,border=1)
    
    lines(ll$doy2,ll$length/24,col=tw.colour,lwd=2.5)
    
    t4t <- t4[t4$year2 == year,]
    plot(t4t$doy2,t4t$mean.temp,type="l",xlim=c(40,304),axes=F,
         xaxt="n",ylim=range(t3$mean.temp),col=1,ylab="temperature")
    y <- c(t4t$min.temp,t4t$max.temp[order(t4t$doy2,decreasing = T)])
    x <- c(t4t$doy2,t4t$doy2[order(t4t$doy2,decreasing = T)])
    polygon(x,y,col=grey(0.8),border=grey(0.8))
    lines(t4t$doy2,t4t$mean.temp,type="l",lwd=2)
    
    lines(df$doy2[df$year2==year], df$value[df$year2==year], type="l",col=tw.colour,lwd=2.5)
    axis(2,las=1)
    
    plot(t4t$doy2, t4t$mean.burst.diff,type="l",xlim=c(40,304),axes=F,
         xaxt="n",ylim=range(t3$mean.burst.diff),col=1,ylab="temp burst diff")
    y <- c(t4t$min.burst.diff,t4t$max.burst.diff[order(t4t$doy2,decreasing = T)])
    x <- c(t4t$doy2,t4t$doy2[order(t4t$doy2,decreasing = T)])
    polygon(x,y,col=grey(0.8),border=grey(0.8))
    lines(t4t$doy2,t4t$mean.burst.diff,type="l",lwd=2)
    axis(2,las=1)
    
    plot(t4t$doy2, t4t$mean.temp.diff,type="l",xlim=c(40,304),axes=F,
         xaxt="n",ylim=range(t3$mean.temp.diff),col=1,ylab="temp max-min")
    y <- c(t4t$min.temp.diff,t4t$max.temp.diff[order(t4t$doy2,decreasing = T)])
    x <- c(t4t$doy2,t4t$doy2[order(t4t$doy2,decreasing = T)])
    polygon(x,y,col=grey(0.8),border=grey(0.8))
    lines(t4t$doy2,t4t$mean.temp.diff,type="l",lwd=2)
    axis(2,las=1)
    
    
    # timing of birth with temp data? -> refer to whiteman 2012 science
    t4t <- t4[t4$year2 == year,]
    plot(t4t$doy2,t4t$mean.temp.var,type="l",xlim=c(40,304),axes=F,
         xaxt="n",ylim=range(t3$mean.temp.var),col=1,ylab="temp var")
    y <- c(t4t$min.temp.var,t4t$max.temp.var[order(t4t$doy2,decreasing = T)])
    x <- c(t4t$doy2,t4t$doy2[order(t4t$doy2,decreasing = T)])
    polygon(x,y,col=grey(0.8),border=grey(0.8))
    lines(t4t$doy2,t4t$mean.temp.var,type="l",lwd=2)
    axis(2,las=1)
    
    plot(t4t$doy2, t4t$sum.wet ,type="l",xlim=c(40,304),axes=F,
         xaxt="n",ylim=c(0,1),ylab="daily prop in seawater",col="transparent", xlab="")
    y <- c(t4t$sum.wet/max(t4$sum.wet),0,0)
    x <- c(t4t$doy2,max(t4t$doy2),min(t4t$doy2))
    polygon(x,y,col=1,border=1)
    axis(2,at=c(0,1),labels = c(0,1),las=1)
    axis(1,at=m$doy2,labels = strftime(m$date,"%b"),cex=2)
    
    dev.off()
  }
}




# plot paper figure with example patterns
png(paste('figures/figure 2 comp v3.png',sep=""), width = 40, height = 30, units="cm", res=600)

par(mfrow=c(3,3),oma=c(2, 5, 3, 0),mar=rep(0.8, 4))
id.years <- c("23980.2014", '23980.2013', '23882.2013')
xlims <- c(45,305) #xlims <- c(40,304)
cex.ax<- 1.6
main.col = grey(0.25)
tw.colour <- brewer.pal(8, "Set2")[6]

d2 <- c(105, 238, 270)
# d2 <- c(105, 270)
d3 <- c(216, 252)

for(iy in id.years){
  
  id   <- str_split_fixed(iy, "[.]", 2)[,1]
  year <- as.numeric(str_split_fixed(iy, "[.]", 2)[,2])
  
  l4 <- l3[l3$animal_id == id,]
  l4 <- l4[order(l4$date),]
  
  tlys <- l4[l4$year2 == year,]
  ll <- data.frame(date = tlys$date,
                   rise=twilight(tm = as.POSIXct(paste(tlys$date,"12:00:00"),tz="UTC"), 
                                 lon = meta2$Capture.GLS.off.lon[meta2$Bear_ID == tlys$animal_id[1]][1], 
                                 lat = meta2$Capture.GLS.off.lat[meta2$Bear_ID == tlys$animal_id[1]][1], 
                                 rise = T, zenith = 93.55),
                   set=twilight(tm = as.POSIXct(paste(tlys$date,"12:00:00"),tz="UTC"), 
                                lon = meta2$Capture.GLS.off.lon[meta2$Bear_ID == tlys$animal_id[1]][1], 
                                lat = meta2$Capture.GLS.off.lat[meta2$Bear_ID == tlys$animal_id[1]][1], 
                                rise = F, zenith = 93.55))
  ll$length <- abs(as.numeric(difftime(ll$rise,ll$set,units="hours")))
  ll$doy <- as.numeric(strftime(ll$date,"%j"))
  ll$doy2     <- ll$doy - 182
  ll$doy2[ll$doy2 < 1] <- ll$doy2[ll$doy2 < 1] +366
  ll$length[is.na(ll$length) & (ll$doy < 265 & ll$doy > 81)] <- 24
  ll$length[is.na(ll$length) & (ll$doy > 265 | ll$doy < 81)] <- 0
  ll$length[ll$length < 12 & (ll$doy < 265 & ll$doy > 81)] <- 24 - ll$length[ll$length < 12 & (ll$doy < 265 & ll$doy > 81)] 
  # ll$length[ll$length > 12 & (ll$doy > 265 | ll$doy < 81)] <- 24 - ll$length[ll$length > 12 & (ll$doy > 265 | ll$doy < 81)]
  ll <- ll[order(ll$doy2),]
  
  
  m <- data.frame(date=as.Date(0:365, paste0(year,"-01-01")))
  m$month <- as.numeric(strftime(m$date,"%m"))
  m$doy      <- as.numeric(strftime(m$date,"%j"))
  m$doy2     <- m$doy - 182
  m$doy2[m$doy2 < 1] <- m$doy2[m$doy2 < 1] +366
  m <- m[!duplicated(m$month),]
  
  t4l <- l4[l4$year2 == year,]
  plot(t4l$doy2,t4l$max.light,type="l",xlim=xlims,
       xaxt="n",ylim=c(0,1),col=grey(0.8),lty=3,ylab='',yaxt="n",axes=F,yaxs="i")
  y <- c(t4l$max.light,0,0)/max(t4l$max.light)
  x <- c(t4l$doy2,max(t4l$doy2),min(t4l$doy2))
  polygon(x,y,col=grey(0.8),border=grey(0.8))
  
  if(iy == "23980.2014") axis(2,at=c(0,1),labels = c(0,1),las=1,cex.axis=cex.ax) else axis(2,at=c(0,1),labels = rep('',2),las=1)
  # mtext("light intensity", side = 2, outer = T, line = 2)
  y <- c(t4l$mean.light,0,0)/max(t4l$max.light)
  x <- c(t4l$doy2,max(t4l$doy2),min(t4l$doy2))
  polygon(x,y,col=main.col,border=main.col)
  
  lines(ll$doy2,ll$length/24,col=tw.colour,lwd=3.5)
  
  
  axis(1,at=m$doy2,labels = rep("", length(m$doy2)),cex=2)
  if(iy == '23980.2013') abline(v = d2, lty = c(2,3,2), lwd = c(2,1,2), col = "darkred")
  # if(iy == '23980.2013') abline(v = d2, lty = c(2,2), lwd = 2, col = "darkred")
  if(iy == '23882.2013') abline(v = d3, lty = c(2,2),   lwd = 2, col = "darkred")
  if(iy == '23980.2013') axis(at = d2[c(1,3)], lab= rep("",2),lwd = 0, lwd.ticks = 2, lty = c(2), col.ticks = "darkred", side = 1, tck=-0.07)
  if(iy == '23980.2013') axis(at = d2[c(2)], lab= rep("",1),lwd = 0, lwd.ticks = 1, lty = c(3), col.ticks = "darkred", side = 1, tck=-0.07)
  if(iy == '23882.2013') axis(at = d3, lab= rep("",length(d3)),lwd = 0, lwd.ticks = 2, lty = c(2,2), col.ticks = "darkred", side = 1, tck =-0.07)
  
}
for(iy in id.years){
  
  id   <- str_split_fixed(iy, "[.]", 2)[,1]
  year <- as.numeric(str_split_fixed(iy, "[.]", 2)[,2])
  t4   <- t3[t3$animal_id == id,]
  
  m <- data.frame(date=as.Date(0:365, paste0(year,"-01-01")))
  m$month <- as.numeric(strftime(m$date,"%m"))
  m$doy      <- as.numeric(strftime(m$date,"%j"))
  m$doy2     <- m$doy - 182
  m$doy2[m$doy2 < 1] <- m$doy2[m$doy2 < 1] +366
  m <- m[!duplicated(m$month),]
  
  t4t <- t4[t4$year2 == year,]
  plot(t4t$doy2,t4t$mean.temp,type="l",xlim=xlims,axes=F,
       xaxt="n",ylim=range(t3$mean.temp),col=1,ylab='')
  y <- c(t4t$min.temp,t4t$max.temp[order(t4t$doy2,decreasing = T)])
  x <- c(t4t$doy2,t4t$doy2[order(t4t$doy2,decreasing = T)])
  polygon(x,y,col=grey(0.8),border=grey(0.8))
  lines(t4t$doy2,t4t$mean.temp,type="l",lwd=2, col=main.col)
  
  lines(df$doy2[df$year2==year], df$value[df$year2==year], type="l",col=tw.colour,lwd=2.5)
  
  if(iy == '23980.2013') abline(v = d2, lty = c(2,3,2), lwd = c(2,1,2), col = "darkred")
  # if(iy == '23980.2013') abline(v = d2, lty = c(2,2), lwd = 2, col = "darkred")
  if(iy == '23882.2013') abline(v = d3, lty = c(2,2),   lwd = 2, col = "darkred")
  
  if(iy == "23980.2014") axis(2, at = seq(-40,40, 10), labels = seq(-40,40, 10), las=1,cex.axis=cex.ax) else axis(2, at = seq(-40,40, 10), labels = rep('',9),las=1)
  axis(1,at=m$doy2,labels = rep("", length(m$doy2)),cex=2)
  
  if(iy == '23980.2013') axis(at = d2[c(1,3)], lab= rep("",2),lwd = 0, lwd.ticks = 2, lty = c(2), col.ticks = "darkred", side = 1, tck=-0.07)
  if(iy == '23980.2013') axis(at = d2[c(2)], lab= rep("",1),lwd = 0, lwd.ticks = 1, lty = c(3), col.ticks = "darkred", side = 1, tck=-0.07)
  if(iy == '23882.2013') axis(at = d3, lab= rep("",length(d3)),lwd = 0, lwd.ticks = 2, lty = c(2,2), col.ticks = "darkred", side = 1, tck =-0.07)
  
}
for(iy in id.years){
  
  
  id   <- str_split_fixed(iy, "[.]", 2)[,1]
  year <- as.numeric(str_split_fixed(iy, "[.]", 2)[,2])
  t4   <- t3[t3$animal_id == id,]
  t4t  <- t4[t4$year2 == year,]
  
  m <- data.frame(date=as.Date(0:365, paste0(year,"-01-01")))
  m$month <- as.numeric(strftime(m$date,"%m"))
  m$doy      <- as.numeric(strftime(m$date,"%j"))
  m$doy2     <- m$doy - 182
  m$doy2[m$doy2 < 1] <- m$doy2[m$doy2 < 1] +366
  m <- m[!duplicated(m$month),]
  
  plot(t4t$doy2, t4t$sum.wet ,type="l",xlim=xlims,axes=F,yaxs="i",
       xaxt="n",ylim=c(0,1),ylab='',col="transparent", xlab="")
  y <- c(t4t$sum.wet/max(t4$sum.wet),0,0)
  x <- c(t4t$doy2,max(t4t$doy2),min(t4t$doy2))
  polygon(x,y,col=main.col,border=main.col)
  
  if(iy == '23980.2013'){
    text(d2[1],0.96, "den",   pos = 2, cex = 1.6, col = "darkred")
    text(d2[1],0.90, "entry", pos = 2, cex = 1.6, col = "darkred")
    text(d2[2],0.96, "first", pos = 2, cex = 1.6, col = "darkred")
    text(d2[2],0.90, "opening", pos = 2, cex = 1.6, col = "darkred")
    text(d2[3],0.96, "den",   pos = 4, cex = 1.6, col = "darkred")
    text(d2[3],0.90, "exit",  pos = 4, cex = 1.6, col = "darkred")
  }
  
  if(iy == '23882.2013'){
    text(d3[1],0.96, "den",   pos = 2, cex = 1.6, col = "darkred")
    text(d3[1],0.90, "entry", pos = 2, cex = 1.6, col = "darkred")
    text(d3[2],0.96, "den",   pos = 4, cex = 1.6, col = "darkred")
    text(d3[2],0.90, "exit",  pos = 4, cex = 1.6, col = "darkred")
  }
  if(iy == '23980.2013') abline(v = d2, lty = c(2,3,2), lwd = c(2,1,2), col = "darkred")
  # if(iy == '23980.2013') abline(v = d2, lty = c(2,2), lwd = 2, col = "darkred")
  if(iy == '23882.2013') abline(v = d3, lty = c(2,2),   lwd = 2, col = "darkred")
  
  if(iy == "23980.2014") axis(2,at=c(0, 0.5, 1),labels = c(0, 9, 18),las=1, cex.axis=cex.ax) else axis(2,at=c(0,1),labels = rep('',2),las=1)
  axis(1,at=m$doy2,labels = strftime(m$date,"%b"),cex.axis=cex.ax)
}
mtext("Light intensity",               at=.83,side=2,outer=T,cex=1.5, line = 2.5)  
mtext("Temperature (\u00B0C)",         at=.5 ,side=2,outer=T,cex=1.5, line = 2.5)  
mtext("Daily count of 'wet'",at=.17,side=2,outer=T,cex=1.5, line = 2.5)  

mtext("Non-denning",      at=.17,side=3,outer=T,cex=1.5, line = 0)  
mtext("Maternity denning",at=.5 ,side=3,outer=T,cex=1.5, line = 0)  
mtext("Shelter denning",  at=.83,side=3,outer=T,cex=1.5, line = 0)  

dev.off()
