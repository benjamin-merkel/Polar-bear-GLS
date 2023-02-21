library(readxl)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(lubridate)
library(geosphere)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("data/Full GLS meta data.RData")
gpsdata <- readRDS("data/Bear collar data provided by Clement Dec 2020.rds")

# calc how many ST collars were  deployed on how many females for how long
gps     <- gpsdata[!duplicated(paste(gpsdata$ID.NR, gpsdata$date)),]
gps     <- gps[as.character(gps$ID.NR) %in% meta2$Bear_ID,]
gps     <- gps[gps$year %in% 2012:2021,]
gps$bear_year <- paste(gps$ID.NR, gps$year, sep = "_")
gps$ID.NR <- droplevels(gps$ID.NR)
gps$ctn   <- droplevels(gps$ctn)
summary(as.numeric(table(gps$ID.NR)))/365
length(unique(as.character(gps$ID.NR)))
length(unique(as.character(gps$ctn)))

# How many years did ST collars collected data?
summary(unlist(lapply(tapply(gps$date, gps$ctn, range), FUN = function(x) as.numeric(difftime(x[2],x[1],units="days"))))/365)
hist(unlist(lapply(tapply(gps$date, gps$ctn, range), FUN = function(x) as.numeric(difftime(x[2],x[1],units="days"))))/365, breaks = 5)


# How many years did GLS collected data?
load("data/Full twilight times dataset GLS dusk adjusted.RData")
summary(unlist(lapply(tapply(twl6b$date, twl6b$gls.id, range), FUN = function(x) as.numeric(difftime(x[2],x[1],units="days"))))/365)
hist(unlist(lapply(tapply(twl6b$date, twl6b$gls.id, range), FUN = function(x) as.numeric(difftime(x[2],x[1],units="days"))))/365, breaks = 5)


# calc how many GLS were  deployed on how many females for how long
on         <- meta2[order(meta2$GLS_on_date),]
on         <- on[!duplicated(on$Bear_ID),]
off        <- meta2[order(meta2$GLS_off_date, decreasing = T),]
off        <- off[!duplicated(off$Bear_ID),]
onoff      <- merge(on[,c("GLS_on_date","Bear_ID")], off[,c("GLS_off_date","Bear_ID")], by = "Bear_ID")
onoff$diff <- as.numeric(onoff$GLS_off_date-onoff$GLS_on_date)
summary(onoff$diff/365)
length(unique(meta2$Bear_ID))

# calc how many females have been double tagged
double     <- meta2[meta2$overlap == 1,]
double     <- double[!duplicated(double$Bear_ID),]
length(unique(double$Bear_ID))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("data/Full GLS location dataset.RData")
load("data/Seasonal GLS centroids.RData")
path2$bear_year <- paste0(path2$bear_id,"_", year(path2$datetime))
path2$season    <- str_split_fixed(path2$bys, " ", 3)[,3]
m2$bear_id      <- str_split_fixed(m2$bear.gls," ", 2)[,1]
m2$year         <- str_split_fixed(m2$bys," ", 3)[,2]
m2$bear_year    <- paste0(m2$bear_id,"_", m2$year)



# % location coverage throughout the year
path_year <- path2[!is.na(path2$lat) & !duplicated(paste(path2$bear.gls, as.Date(path2$datetime))),]
path_year$year <- as.numeric(strftime(path_year$datetime, "%Y"))
meta_year <- meta2[meta2$bear.gls %in% path_year$bear.gls,]
meta_year <- meta_year[order(meta_year$bear.gls),]
summary(as.numeric(table(paste(path_year$bear.gls,path_year$year))))/360
summary(as.numeric(table(paste(path_year$bear.gls))/as.numeric(meta_year$GLS_off_date-meta_year$GLS_on_date)))
#% using just longitudes
path_year <- path2[!duplicated(paste(path2$bear.gls, as.Date(path2$datetime))),]
path_year$year <- as.numeric(strftime(path_year$datetime, "%Y"))
meta_year <- meta2[meta2$bear.gls %in% path_year$bear.gls,]
meta_year <- meta_year[order(meta_year$bear.gls),]
summary(as.numeric(table(paste(path_year$bear.gls,path_year$year))))/360
summary(as.numeric(table(paste(path_year$bear.gls))/as.numeric(meta_year$GLS_off_date-meta_year$GLS_on_date)))


# seasonal centroid error rate
m2g              <- m2[!is.na(m2$gps.lat),]
m2g$distance.gps <- distCosine(matrix(c(m2g$gls.lon, m2g$gls.lat), ncol=2), matrix(c(m2g$gps.lon, m2g$gps.lat), ncol=2))/1000
summary(m2g)
summary(m2g[m2g$del.ratio>0.6,])
length(unique(m2$bear_id))
length(unique(m2g$bear_id))


# daily error rate
path3              <- path2[!is.na(path2$gps.lat) & !is.na(path2$lat),]
path3$distance.gps <- distCosine(matrix(c(path3$lon, path3$lat), ncol=2), matrix(c(path3$gps.lon, path3$gps.lat), ncol=2))/1000
path3$season       <- "autumn"
path3$season[path3$doy < 150] <- "spring"
length(unique(path2$bear_id))
length(unique(path3$bear_id))
boxplot(path3$distance.gps ~ path3$season)
summary(path3)



# mdb <- read_excel("data/Repr hist 2022.xls")
mdb <- read_excel("data/Repr hist 2020c.xls")
mdb$DenGLT <- as.numeric(mdb$DenGLT)#[mdb$DenGLT == 'NA'] <- NA
mdb$DenCol <- as.numeric(mdb$DenCol)#mdb$DenTDR[mdb$DenTDR == 'NA'] <- NA
mdb$DenCap <- as.numeric(mdb$DenCap)#mdb$DenCap[mdb$DenCap == 'NA'] <- NA
mdb$Den    <- as.numeric(mdb$Den)#mdb$Den[mdb$Den == 'NA']       <- NA
mdb        <- mdb[,c("DenGLT","DenCap","DenCol","Den","YEAR","age","Age group","ID NR...6","STED")]
mdb$DenGLT[is.na(mdb$DenGLT) & is.na(mdb$DenCap) & is.na(mdb$DenCol) & !is.na(mdb$Den)] <- mdb$Den[is.na(mdb$DenGLT) & is.na(mdb$DenCap) & is.na(mdb$DenCol) & !is.na(mdb$Den)]
mdb$DenCapGLT                             <- mdb$DenCap
mdb$DenCapGLT[is.na(mdb$DenCapGLT)]       <- mdb$DenGLT[is.na(mdb$DenCapGLT)]
mdb$DenCapCol                             <- mdb$DenCap
mdb$DenCapCol[is.na(mdb$DenCapCol)]       <- mdb$DenCol[is.na(mdb$DenCapCol)]
mdb$DenCapColGLT                          <- mdb$DenCap
mdb$DenCapColGLT[is.na(mdb$DenCapColGLT)] <- mdb$DenCol[is.na(mdb$DenCapColGLT)]
mdb$DenCapColGLT[is.na(mdb$DenCapColGLT)] <- mdb$DenGLT[is.na(mdb$DenCapColGLT)]


# mdb$DenCol[is.na(mdb$DenCol)]  <- mdb$DenTDR[is.na(mdb$DenCol)]
mdb$YEAR      <- as.numeric(mdb$YEAR)
mdb$bear_year <- paste0(mdb$`ID NR...6`,"_",mdb$YEAR)
mdb <- mdb[!is.na(mdb$`ID NR...6`) & !is.na(mdb$YEAR),]
mdb$age <- as.numeric(mdb$age)
mdb     <- mdb[order(mdb$YEAR),]
mdb     <- mdb[order(mdb$`ID NR...6`),]
ids      <- unique(mdb$`ID NR...6`)
mdb$age_diff <- NA
for(id in ids) mdb$age_diff[mdb$`ID NR...6` == id] <- mdb$age[mdb$`ID NR...6` == id] - c(NA, head(mdb$age[mdb$`ID NR...6` == id], -1))

summary(mdb)



# frequency of lon and lat during each season and bear
lat_all <- table(factor(path2$bys[!is.na(path2$lat)], level = unique(path2$bys)))
lat_sp <- table(factor(path2$bys[!is.na(path2$lat) & path2$season=="spring"], level = unique(path2$bys[path2$season=="spring"])))
summary(as.numeric(lat_all))
lon_all <- table(factor(path2$bys, level = unique(path2$bys)))
summary(as.numeric(lon_all))

# frequency of lon and lat during spring for denning bears
path3 <- path2[path2$bear_year %in% mdb$bear_year[mdb$Den %in% 1] & str_split_fixed(path2$bys, " ", 3)[,3] == "spring",]
lat_den <- table(factor(path3$bys[!is.na(path3$lat)], level = unique(path3$bys)))
lat_sp_den <- table(factor(path3$bys[!is.na(path3$lat) & path3$season=="spring"], level = unique(path3$bys[path3$season=="spring"])))
summary(as.numeric(lat_den))
lon_den <- table(factor(path3$bys, level = unique(path3$bys)))
summary(as.numeric(lon_den))

hist(lon_all)


hist(lat_all,border='transparent',col=alpha(1,0.2))
hist(lat_sp,border='transparent',col=alpha(4,0.5),add=T)
hist(lat_sp_den,border='transparent',col=alpha(2,0.5),add=T)


# lon and lat during spring for denning bears
path3 <- path2[path2$bear_year %in% mdb$bear_year[mdb$Den %in% 1] & str_split_fixed(path2$bys, " ", 3)[,3] == "spring",]
path3 <- st_as_sf(path3[!is.na(path3$lat),], coords = c("lon", "lat"), crs = 4326)
gps3  <- gps[gps$bear_year %in% unique(path3$bear_year) & as.numeric(gps$month) %in% c(3,4),]
gps3  <- st_as_sf(gps3[!is.na(gps3$x),], coords = c("x", "y"), crs = 4326)
path3 <- path3[path3$bear_year %in% unique(gps$bear_year),]


m3 <- m2[m2$bear_year %in% mdb$bear_year[mdb$Den %in% 1] & m2$season == "spring" & !is.na(m2$gps.lon),]
m3$distance.gps <- distCosine(matrix(c(m3$gls.lon, m3$gls.lat), ncol=2), matrix(c(m3$gps.lon, m3$gps.lat), ncol=2))/1000
m3 <- st_as_sf(m3, coords = c("gls.lon", "gls.lat"), crs = 4326)


plot(m3$geometry)
plot(land,add=T)
for(i in 1:nrow(m3)) lines(c(m3$gps.lon[i], st_coordinates(m3[i,])[,1]), c(m3$gps.lat[i], st_coordinates(m3[i,])[,2]))
plot(m3$geometry,add=T, bg=1, pch=21)
points(m3$gps.lon, m3$gps.lat, bg=2, pch=21)

mdbx          <- mdb[!is.na(mdb$STED),]
mdbx          <- mdbx[mdbx$`ID NR...6` %in% meta2$Bear_ID,]
mdbx          <- mdbx[mdbx$YEAR %in% c(2013:2021),]
summary(as.numeric(table(mdbx$`ID NR...6`))[as.numeric(table(mdbx$`ID NR...6`))!=1])

mdb2          <- mdb[mdb$`Age group` %in% "ad" | mdb$age > 5,]
mdb2          <- mdb2[mdb2$`ID NR...6` %in% meta2$Bear_ID,]

summary(mdb2)

mdb2$DenGLT[is.na(mdb2$DenGLT)] <- -1
mdb2          <- mdb2[order(mdb2$DenGLT, decreasing = T),]
mdb2          <- mdb2[!duplicated(mdb2$bear_year),]
mdb2$DenGLT[mdb2$DenGLT == -1] <- NA
mdb3          <- mdb2[mdb2$YEAR %in% c(2013:2021),]
mdb3          <- mdb3[order(mdb3$YEAR),]

summary(mdb3)
summary(mdb3[is.na(mdb3$DenCap) & is.na(mdb3$DenCol) & is.na(mdb3$DenGLT),])

mdb_GLS_CAP    <- data.frame(mdb3[!is.na(mdb3$DenGLT) & !is.na(mdb3$DenCap),])
# remove 23992_2015 as it is odd in Jons file
mdb_GLS_CAP    <- mdb_GLS_CAP[mdb_GLS_CAP$bear_year != "23992_2015",]
mdb_COL_CAP    <- mdb3[!is.na(mdb3$DenCol) & !is.na(mdb3$DenCap),] 
mdb_GLS_COL    <- data.frame(mdb3[!is.na(mdb3$DenGLT) & is.na(mdb3$DenCap) & !is.na(mdb3$DenCol),])
# remove 23992_2018
mdb_GLS_COL    <- mdb_GLS_COL[mdb_GLS_COL$bear_year != "23992_2018",]


all(mdb_GLS_CAP$DenGLT==mdb_GLS_CAP$DenCap)
all(mdb_COL_CAP$DenCol==mdb_COL_CAP$DenCap)
all(mdb_GLS_COL$DenGLT==mdb_GLS_COL$DenCol)


ids <- unique(mdb3$`ID NR...6`)
for(id in ids){
  mm <- as.data.frame(mdb3[mdb3$`ID NR...6` == id, c("YEAR", "DenCap", "DenCol", "DenGLT")])
  mm$yeardiff <- mm$YEAR - c(mm$YEAR[1]-1, head(mm$YEAR, - 1))
  mm$recCap <- mm$recCol <- mm$recTDR <- mm$recGLT <- NA
  mm$recCap[mm$yeardiff==1 & !is.na(mm$DenCap)] <- 1
  mm$recCol[mm$yeardiff==1 & !is.na(mm$DenCol)] <- 1
  mm$recTDR[mm$yeardiff==1 & !is.na(mm$DenTDR)] <- 1
  mm$recGLT[mm$yeardiff==1 & !is.na(mm$DenGLT)] <- 1
  
  m2 <- mm[1, 6:9]
  m2$recCapGLT <- m2$recCapColTDR <- m2$recCapColGLT <- m2$recCapCol <- m2$recCap <- m2$recCol <- m2$recTDR <- m2$recGLT <- 0  
  m3 <- m2
  # mcap <- mm[mm$DenCap %in% c(0,1),]
  for(l in 1:nrow(mm)){
    
    m3$recCap    <- max(c(m2$recCap, m3$recCap))
    m3$recCol    <- max(c(m2$recCol, m3$recCol))
    m3$recTDR    <- max(c(m2$recTDR, m3$recTDR))
    m3$recGLT    <- max(c(m2$recGLT, m3$recGLT))
    m3$recCapCol <- max(c(m2$recCapCol, m3$recCapCol))
    m3$recCapGLT <- max(c(m2$recCapGLT, m3$recCapGLT))
    m3$recCapColGLT <- max(c(m2$recCapColGLT, m3$recCapColGLT))
    m3$recCapColTDR <- max(c(m2$recCapColTDR, m3$recCapColTDR))
    
    if(mm$recCap[l] %in% 1) m2$recCap    <- m2$recCap+1    else m2$recCap    <- 0 
    if(mm$recCap[l] %in% 1 |
       mm$recGLT[l] %in% 1) m2$recCapGLT <- m2$recCapGLT+1 else m2$recCapGLT <- 0 
    if(mm$recCap[l] %in% 1 |
       mm$recCol[l] %in% 1) m2$recCapCol <- m2$recCapCol+1 else m2$recCapCol <- 0 
    if(mm$recCap[l] %in% 1 |
       mm$recCol[l] %in% 1 |
       mm$recGLT[l] %in% 1) m2$recCapColGLT <- m2$recCapColGLT+1 else m2$recCapColGLT <- 0 
    if(mm$recCap[l] %in% 1 | mm$recTDR[l] %in% 1 |
       mm$recCol[l] %in% 1) m2$recCapColTDR <- m2$recCapColTDR+1 else m2$recCapColTDR <- 0 
    if(mm$recCol[l] %in% 1) m2$recCol    <- m2$recCol+1    else m2$recCol    <- 0 
    if(mm$recTDR[l] %in% 1) m2$recTDR    <- m2$recTDR+1    else m2$recTDR    <- 0 
    if(mm$recGLT[l] %in% 1) m2$recGLT    <- m2$recGLT+1    else m2$recGLT    <- 0 
  }
  m3$recCap    <- max(c(m2$recCap, m3$recCap))
  m3$recCol    <- max(c(m2$recCol, m3$recCol))
  m3$recTDR    <- max(c(m2$recTDR, m3$recTDR))
  m3$recGLT    <- max(c(m2$recGLT, m3$recGLT))
  m3$recCapCol <- max(c(m2$recCapCol, m3$recCapCol))
  m3$recCapGLT <- max(c(m2$recCapGLT, m3$recCapGLT))
  m3$recCapColGLT <- max(c(m2$recCapColGLT, m3$recCapColGLT))
  m3$recCapColTDR <- max(c(m2$recCapColTDR, m3$recCapColTDR))
  m3$bear_id   <- id
  m3$n_years   <- nrow(mm)
  if(id == ids[1]) m4 <- m3 else m4 <- rbind(m4, m3)
}
summary(m4)



m5a <- m4[,c("recGLT","bear_id","n_years")]
colnames(m5a) <- c("freq","bear_id","n_years")
m5a$type <- "GLS"
m5b <- m4[,c("recCap","bear_id","n_years")]
colnames(m5b) <- c("freq","bear_id","n_years")
m5b$type <- "Capture"
m5c <- m4[,c("recCol","bear_id","n_years")]
colnames(m5c) <- c("freq","bear_id","n_years")
m5c$type <- "Collar"
m5d <- m4[,c("recCapCol","bear_id","n_years")]
colnames(m5d) <- c("freq","bear_id","n_years")
m5d$type <- "Cap+Col"
m5e <- m4[,c("recCapGLT","bear_id","n_years")]
colnames(m5e) <- c("freq","bear_id","n_years")
m5e$type <- "Cap+GLS"
m5f <- m4[,c("recCapColGLT","bear_id","n_years")]
colnames(m5f) <- c("freq","bear_id","n_years")
m5f$type <- "Cap+GLS+GPS"
m5 <- rbind(m5a, m5b, m5c, m5d, m5e, m5f)


png(paste0("figures/years of unbroken records boxplot.png"), res = 800, width=8, height = 26, units="cm")
opar <- par(mar=c(8,4,2,1))
boxplot(m4$recCapGLT, ylim = c(0, 8), boxwex=0.15, xlim = c(0.85, 1.45), col=grey(1), ylab = "years of unbroken record", las = 2)
boxplot(m4$recCapColGLT,add=T,at=0.9, boxwex=0.15, col=grey(1), axes = F)
boxplot(m4$recCapCol,add=T,at=1.1, boxwex=0.15, col=grey(1), axes = F)
boxplot(m4$recGLT,add=T,at=1.2, boxwex=0.15, col=grey(1), axes = F)
boxplot(m4$recCap,add=T,at=1.3, boxwex=0.15, col=grey(1), axes = F)
boxplot(m4$recCol,add=T,at=1.4, boxwex=0.15, col=grey(1), axes = F)
axis(1, at = seq(0.9,1.4, 0.1), labels = c("Capture +\nGLS + ST","Capture +\nGLS", "Capture +\nST collar", "GLS", "Capture", "ST collar"), las=3)
par(opar)
dev.off()


png(paste0("figures/years of unbroken records density.png"), res = 800, width=20, height = 20, units="cm")
ggplot(m5, aes(x=freq, fill=type), xlim = c(0,8)) +
  geom_density(alpha=0.2, adjust = 1.2) +
  theme_classic()#+  
dev.off()


table(mdb3$DenGLT)
table(mdb3$DenCap)
table(mdb3$DenCol)
table(mdb3$DenCapGLT)
table(mdb3$DenCapCol)
table(mdb3$DenCapColGLT)
table(mdb3$Den)

table(mdb3$DenCap, mdb3$DenGLT)
table(mdb3$DenCap, mdb3$DenCol)

mdb4 <- mdb3[is.na(mdb3$DenCap) & !is.na(mdb3$DenGLT) & !is.na(mdb3$DenCol),]
table(mdb4$Den)
# I disregarded 23992 - 2018 
table(mdb4$DenGLT, mdb4$DenCol)

table(mdb4$YEAR)
length(unique(mdb4$`ID NR...6`))

summary(mdb3)


mdb4 <- mdb3[!is.na(mdb3$DenGLT) & !is.na(mdb3$DenCol),]
table(mdb4$YEAR)
length(unique(mdb4$`ID NR...6`))
