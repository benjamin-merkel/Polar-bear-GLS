library(zoo)

MYfindHEZenith <- function(twl, tol = 0.08, range=c(250,400), ylim = c(0,75)){
  
  z <- seq(89, 99, by = 0.1)
  
  lats  <- apply(cbind(z), 1, function(x) thresholdPath(twl$Twilight, twl$Rise, zenith = x, tol=tol)$x[,2])
  sds   <- apply(lats[range[1]:range[2],], 2, sd, na.rm = T)
  
  colsT <- data.frame(sd = seq(min(sds)-0.1, max(sds)+0.1, length = 101), col = heat.colors(101))
  
  opar <- par(mfrow = c(2,1), mar=c(4,4,1,1))
  matplot(lats, col = as.character(colsT$col[cut(sds, breaks = colsT[,1], labels = F)]),
          lty = 1, type = "l", xlab = "", ylab = "Latitude", las = 1, xaxt = "n",ylim = ylim)
  lines(lats[,which.min(sds)], lwd = 3)
  abline(v = range, lty = 2, col = "cornflowerblue")
  axis(1, at = seq(1, nrow(twl), length = 5),
       labels = format(as.POSIXct(seq(twl$Twilight[1], twl$Twilight[nrow(twl)], length = 5)), "%d-%b"))
  plot(z, sds, las = 1, type = "o", pch = 16, cex = 1.3, col = as.character(colsT$col[cut(sds, breaks = colsT[,1], labels = F)]),
       xlab = "zenith", ylab = "sd in latitude (within range)")
  
  abline(v = z[which.min(sds)], lty = 3)
  mtext(paste("zenith =", z[which.min(sds)]), 3, line = -1.5, cex = 1.1)
  par(opar)
  
  return(z[which.min(sds)])
  
}

MytwilightEdit <- function (twilights, offset = 17, window = 4, outlier.mins = 45, 
          stationary.mins = 15, zlim = c(0, 64), plot = T) 
{
  day <- twilights$Twilight
  hour <- hourOffset(as.hour(twilights$Twilight), offset)
  sunr <- which(twilights$Rise)
  suns <- which(!twilights$Rise)
  
  
  fnc <- function(x) {
    ind0 <- abs(x[(window/2) + 1, 1] - x[-((window/2) + 1), 1]) > median(diff(as.numeric(day[suns]))) * ((window/2) + 1)
    if (any(ind0)) x <- x[-((window/2) + 1), ][-which(ind0), ]
    
    empt <- F
    if(!is.null(nrow(x))) {if(nrow(x)>1)  empt <- T}
     
    
    if(empt){
      if (nrow(x) < ((window/2) + 1)) {
        out <- cbind(x[(window/2), 1], 
                     FALSE, 
                     FALSE, 
                     x[(window/2), 1])
      }
      else {
        diffr <- abs(x[(window/2) + 1, 2] - median(x[-((window/2) + 1), 2])) * 60
        
        if (diffr >= outlier.mins & all(dist(x[-((window/2) +  1), 2]) * 60 <= stationary.mins)) {
          out <- cbind(x[(window/2) + 1, 1] + (median(x[c(window/2, (window/2) + 1), 2]) - x[(window/2) + 1, 2]) * 60 * 60, 
                       FALSE, 
                       TRUE, 
                       x[(window/2) + 1, 1])
        }
        if (diffr >= outlier.mins & !all(dist(x[-((window/2) + 1), 2]) * 60 <= stationary.mins)) {
          out <- cbind(x[(window/2) + 1, 1], 
                       TRUE, 
                       FALSE, 
                       x[(window/2) + 1, 1])
        }
        if (diffr < outlier.mins) 
          out <- cbind(x[(window/2) + 1, 1], 
                       FALSE, 
                       FALSE, 
                       x[(window/2) + 1, 1])
      }
      
    } else {
      out <- cbind(x[1], 
                   FALSE, 
                   FALSE, 
                   x[1])
      }
    return(out)
  }
  
  
  # xr <- cbind(day, hour)[sunr, ]
  
  sunrT <- data.frame(rollapply(cbind(day, hour)[sunr, ], width = window + 1, FUN = fnc, fill = FALSE, by.column = F))
  sunsT <- data.frame(rollapply(cbind(day, hour)[suns, ], width = window + 1, FUN = fnc, fill = FALSE, by.column = F))
  sunrT[which(sunrT[, 1] == 0), c(1, 4)] <- day[sunr][which(sunrT[, 1] == 0)]
  sunsT[which(sunsT[, 1] == 0), c(1, 4)] <- day[suns][which(sunsT[, 1] == 0)]
  sunsT$X5 <- F
  sunrT$X5 <- T
  out <- data.frame(Twilight = as.POSIXct(c(sunrT[, 1], sunsT[,1]), origin = "1970-01-01", tz = "GMT"), 
                    Rise = c(sunrT[,5],sunsT[,5]), 
                    Deleted = ifelse(c(sunrT[, 2], sunsT[, 2]) == 1, TRUE, FALSE), 
                    Edited  = ifelse(c(sunrT[, 3], sunsT[, 3]) == 1, TRUE, FALSE), 
                    Twilight0 = as.POSIXct(c(sunrT[, 4], sunsT[, 4]), origin = "1970-01-01", tz = "GMT"))
  out <- out[order(out[, 1]), ]
  rownames(out) <- 1:nrow(out)
  if (plot) {
    day0 <- out$Twilight0
    hour0 <- hourOffset(as.hour(out$Twilight0), offset)
    day <- out$Twilight
    hour <- hourOffset(as.hour(out$Twilight), offset)
    plot(day0, hour0, type = "n", xlab = "Date", ylab = "Hour")
    points(day[!out$Deleted], hour[!out$Deleted], pch = 16, cex = 0.5, col = ifelse(out$Rise[!out$Deleted], "firebrick", "cornflowerblue"))
    arrows(day0[out$Edited], hour0[out$Edited], day[out$Edited], hour[out$Edited], length = 0.1)
    points(day0[out$Deleted | out$Edited], hour0[out$Deleted | out$Edited], pch = 16, col = "grey50")
    points(day[out$Edited], hour[out$Edited], pch = 16, col = ifelse(out$Rise[out$Edited], "firebrick", "cornflowerblue"))
    points(day0[out$Deleted], hour0[out$Deleted], pch = "X")
  }
  out
}

# adjust function to work also around midnight sun conditions
MythresholdCalibration <- function (twilight, rise, lon, lat, method = "log-normal", 
                                    plot = TRUE) 
{
  if (!(method %in% c("gamma", "log-norm"))) 
    stop("Method can only be `gamma` or `log-norm`.")
  tab <- data.frame(Twilight = twilight, Rise = rise)
  sun <- solar(tab[, 1])
  z <- refracted(zenith(sun, lon, lat))
  inc = -1
  repeat {
    cat("\r",inc,"  ")
    twl_t <- twilight(tab[, 1], lon, lat, rise = tab[, 2], 
                      zenith = max(z) + inc)
    twl_dev <- ifelse(tab$Rise, as.numeric(difftime(tab[,1], twl_t, units = "mins")), 
                      as.numeric(difftime(twl_t, tab[, 1], units = "mins")))
    if (all(twl_dev >= 0, na.rm=T)) {
      break
    }
    else {
      inc <- inc + 0.01
    }
  }
  z0 <- max(z) + inc
  seq <- seq(0, max(twl_dev,na.rm=T), length = 100)
  if (method == "log-norm") {
    fitml_ng <- suppressWarnings(fitdistr(twl_dev, "log-normal"))
    lns <- dlnorm(seq, fitml_ng$estimate[1], fitml_ng$estimate[2])
  }
  if (method == "gamma") {
    fitml_ng <- suppressWarnings(fitdistr(twl_dev[!is.na(twl_dev)], "gamma"))
    lns <- dgamma(seq, fitml_ng$estimate[1], fitml_ng$estimate[2])
  }
  diffz <- as.data.frame(cbind(min = apply(cbind(tab[, 1], 
                                                 twilight(tab[, 1], lon, lat, rise = tab[, 2], zenith = z0)), 
                                           1, function(x) abs(x[1] - x[2]))/60, z = z))
  mod <- lm(z ~ min, data = diffz)
  mod2 <- lm(min ~ z, data = diffz)
  z1 <- median(z)
  if (plot) {
    opar <- par(mar = c(10, 4, 1, 1))
    if (method == "log-norm") 
      hist(twl_dev, freq = F, breaks = 26, main = "Twilight Model (log-norm)", 
           xlab = "twilight error (min)")
    if (method == "gamma") 
      hist(twl_dev, freq = F, breaks = 26, main = "Twilight Model (gamma)", 
           xlab = "twilight error (min)")
    lines(seq, lns, col = "firebrick", lwd = 3, lty = 2)
    points(predict(mod2, newdata = data.frame(z = z0)), 0, 
           pch = 21, cex = 5, bg = "white", lwd = 2)
    text(predict(mod2, newdata = data.frame(z = z0)), 0, 
         "0")
    points(z1, 0, pch = 21, cex = 5, bg = "white", 
           lwd = 2)
    text(z1, 0, "1")
    axis(1, at = seq(0, max(twl_dev,na.rm=T), 6), 
         labels = round(90 - predict(mod, newdata = data.frame(min = seq(0, max(twl_dev,na.rm=T), 6))), 1), line = 5)
    mtext("sun elevation angle (degrees)", 1, line = 8)
    if (method == "log-norm") 
      legend("topright", paste(c("0. Zenith angle (zero)", 
                                 "1. Zenith angle (median)", "log-mean", 
                                 "log-sd"), round(c(z0, z1, fitml_ng$estimate[1], 
                                                    fitml_ng$estimate[2]), 3)), bty = "n")
    if (method == "gamma") 
      legend("topright", paste(c("0. Zenith angle (zero)", 
                                 "1. Zenith angle (median)", "shape", 
                                 "scale"), round(c(z0, z1, fitml_ng$estimate[1], 
                                                   fitml_ng$estimate[2]), 3)), bty = "n")
    par(opar)
  }
  pp <- predict(mod2, newdata = data.frame(z = z0))
  pp.mode <- seq[lns==max(lns, na.rm=T)]
  pp.mode <- predict(mod, newdata = data.frame(min = pp.mode))
  
  c(a1 = z1, e0 = z0, pred.mode = pp.mode, log.mean = fitml_ng$estimate[1], log.sd = fitml_ng$estimate[2], 
    n = length(twl_dev[!is.na(twl_dev)]), n.NA = length(twl_dev[is.na(twl_dev)]))
  # list(c(a1 = z1, e0 = z0, log.mean = fitml_ng$estimate[1], log.sd = fitml_ng$estimate[2]),twl_dev)
}


MytsimagePoints <- function (date, offset, ...) {
  points(date, hourOffset(as.hour(date), offset%%24)-24, ...)
}


Seasonal_palette <- grDevices::colorRampPalette(grDevices::hsv(1 - ((1:365) + (365/4))%%365/365, s = 0.8, v = 0.8), space = "Lab")
