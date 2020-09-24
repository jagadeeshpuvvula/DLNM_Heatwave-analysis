#LOAD LIBS
library(dlnm)
library(splines)

#LOAD DATA
coastal <- read.csv ("C:/Users/jagad/Desktop/NC_Sur/NC_Dec2019/coastal_dec.csv", header = T)

# DATA FORMATTING
cat <- c("tavg_95_2_coas","tavg_98_2_coas","tavg_95_3_coas","Tavg_98_3_coas","tavg_90_2_coas",
         "tavg_90_3_coas","tmax_98_3_coas","tmax_98_2_coas","tmax_90_2_coas","tmax_90_3_coas",
         "tmax_95_2_coas","tmax_95_3_coas","tmin_95_3_coas","tmin_98_2_coas","Tmin_98_3_coas",
         "tmin_90_2_coas","tmin_90_3_coas","tmin_95_2_coas", "dow", "wDay")
coastal[cat] <- lapply(coastal[cat], factor)

coastal$ADMIT_SUP <- as.numeric(coastal$ADMIT_SUP)
coastal$tmax <- as.numeric(coastal$tmax)
coastal$tmin <- as.numeric(coastal$tmin)
coastal$tavg <- as.numeric(coastal$tavg)
coastal$date <- as.Date(coastal$date, format = "%m/%d/%Y")
coastal$dow <- wday(as.Date(coastal$date, format = "%m/%d/%Y"))

weekdays1 <- c('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday')
coastal$wDay <- factor((weekdays(coastal$date) %in% weekdays1), 
                       levels=c(FALSE, TRUE), labels=c('weekend', 'weekday'))

coastal$doy <- as.numeric(format(coastal$date, "%d"))
coastal$month <- as.numeric(format(coastal$date, "%m"))
coastal$year <- as.numeric(format(coastal$date, "%Y"))
coastal$dif <- as.numeric(coastal$tmax-coastal$tmin)

# DLNM model implementation
varknots <- equalknots(coastal$tmax,fun="bs",df=5,degree=2)
lagknots <- logknots(30, 3)
cb3.temp <- crossbasis(coastal$tmax, lag=30, argvar=list(fun="bs",
                                                         knots=varknots),
                       arglag=list(knots=lagknots))

# GLM object
model3 <- glm(ADMIT_SUP ~ tmax,
              family=quasipoisson(), coastal)

#wrap in DLNM framework
pred3.temp <- crosspred(cb3.temp, model3, cen=60, by=1)

#plot results
plot(pred3.temp, xlab="Temperature", zlab="RR", theta=200, phi=40, lphi=30,
     main="3D graph of temperature effect")
####################################################################3
#####################################################################

#DLNM example 2
library(lubridate)
coastal$year <- year(as.Date(coastal$date, format = "%m/%d/%Y"))

#Cross basis function
cb2.temp <- crossbasis(coastal$tmax, lag=10,
                       argvar=list(fun="thr",thr=c(60,77)), 
                       arglag=list(fun="strata",breaks=c(2,6)),
                       group=coastal$year)
#GLM object
model2 <- glm(ADMIT_SUP ~ tmax,
              family=quasipoisson(), coastal)
pred2.o3 <- crosspred(cb2.temp, model2, at=c(60:100))

#Ref: https://cran.r-project.org/web/packages/dlnm/vignettes/dlnmTS.pdf
