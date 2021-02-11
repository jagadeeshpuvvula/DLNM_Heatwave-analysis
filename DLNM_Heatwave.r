library(dlnm)
library(splines)
library(lubridate)
library(mgcv)
#=================================================================================
#coastal
dat <- read.csv ("C:/Users/jagad/Desktop/NC_Sur/NC_Dec2019/dat_octf.csv", header = T,
                 fileEncoding="UTF-8-BOM")
#piedmont
dat <- read.csv ("C:/Users/jagad/Desktop/NC_Sur/NC_Dec2019/pied_octf.csv", header = T,
                 fileEncoding="UTF-8-BOM")
dat$date <- dat$Date # only for piedmont
# ==================================================================================
dat$date <- as.Date(dat$date, format = "%m/%d/%Y")
dat$dow <- wday(as.Date(dat$date, format = "%m/%d/%Y"))
dat$month<- as.factor(month(dat$date))
dat$year<- as.factor(format(dat$date, '%Y'))
weekdays1 <- c('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday')
dat$wDay <- factor((weekdays(dat$date) %in% weekdays1), 
                       levels=c(FALSE, TRUE), labels=c('weekend', 'weekday'))

cat <- c("month", "dow", "wDay", "year")
dat[cat] <- lapply(dat[cat], factor)
dat$Rate_ER_visit <- as.numeric(dat$imp_rate)
dat$Count_ER_visit <- as.numeric(dat$imp_count)
dat$Max_temp <- as.numeric(dat$tmax)
# ==================================================================================
varknots <- equalknots(dat$tmax,fun="bs",df=5,degree=2)
lagknots <- logknots(5, 1)
cb3.temp <- crossbasis(dat$Max_temp, lag=5, 
                       argvar=list(fun="bs",knots=varknots), arglag=list(knots=lagknots))
model3 <- gam(Count_ER_visit ~ cb3.temp+dow+month+year,
              family=quasipoisson(), dat)
# ==================================================================================
#how to set the cen value
pred3.temp <- crosspred(cb3.temp, model3,coef=NULL, vcov = NULL, cen=90, by=1, 
                        from = 64, to = 98)

#Print risk ratio matrix
print(pred3.temp$allRRfit)

#3d plot - mostly useless
plot(pred3.temp, xlab="Temperature", zlab="RR", theta=200, phi=40, lphi=30,
     main="3D graph of temperature effect")

#relative risk | Temperature | lag
plot(pred3.temp, "contour", xlab="Temperature", key.title=title("RR"),
     plot.title=title("Contour plot",xlab="Temperature",ylab="Lag"))

#Relative risk and lag days
plot(pred3.temp, "slices", var=95, ci="n", col=1, ylim=c(0.95,1.25), lwd=1.5,
     main="Lag-response curves")

#what happens at different lags and temp values [change var=c() for temp | lag=c() for lag days] 
plot(pred3.temp, "slices", var=c(85,95), lag=c(0,5), col=4,
     ci.arg=list(density=40,col=grey(0.7)))

#what happens when exposed to 95 deg F (change var for change in temp)
plot(pred3.temp, "slices", var=95, ci="bars", type="p", col=2, pch=19,
     ci.level=0.95, main="Lag-response a 10-unit increase above threshold (95CI)")

# ==================================================================================
