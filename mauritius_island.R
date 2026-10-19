getwd()
######################################################
# Analysis of stationarity on Mauritius Island dataset
#By Franchezka Huaranca and Maëlle Anglair
######################################################

#------------
# Import data
#------------

install.packages("readxl")
library(readxl)

MauritiusTraf <- read_excel("Mauritius-Island-2008-2017.xlsx")
str(MauritiusTraf)

traf<- ts(MauritiusTraf[,2], start=c(2008,01),frequency=12)

#---------------------
# Is traf stationary??: We have 3 tools: plot of series, ACF and PACF.
#---------------------
plot.ts(traf)
# we observe several patterns of non stationarity:
# increasing variance
# trend
# regular oscillations (seasonal effects)

par(mfrow=c(2,1))
acf(ts(traf, frequency=1), main="Autocorrelation main series")
pacf(ts(traf, frequency=1), main="Partial autocorrelation main series")
par(mfrow=c(1,1))
# The ACF confirms non stationarity

#-----------------------------------------
# First transformation: log transformation
#-----------------------------------------
ltraf <- log(traf)

plot(ltraf)
acf(ts(ltraf, frequency=1), 
    main="Autocorrelation log (main series)")
pacf(ts(ltraf, frequency=1), 
     main="Partial autocorrelation log(main series)")
# No more increasing variance, but still trend and seasonality
# persistence of significant coeff in the ACF
# The log(series) is still non stat

#--------------------------------------------
# Second transformation: 1st order difference
#--------------------------------------------
dltraf <- diff(ltraf, 1)
plot(dltraf)
acf(ts(dltraf, frequency=1),
    main="Autocorrelation 1st order diff of log (main series)")
pacf(ts(dltraf, frequency=1), 
     main="Partial autocorrelation 1st order diff of log(main series)")

# We observe seasonality on the graph 
# And also on the ACF (no fast decay to 0 of the coefficients)
# the series is still non stationary

#---------------------------------------------
# Third transformation: difference of order 12
#---------------------------------------------
dltraf_12 <- diff(dltraf, 12)
plot(dltraf_12)
acf(ts(dltraf_12, frequency=1), 
    main="Autocorrelation 1st order")
pacf(ts(dltraf_12, frequency=1), 
     main="Partial autocorrelation ")

#We remove the last 2 values of the data series and we save it to a new
#time series object called traf2:

traf2<- ts(MauritiusTraf[,2], start=c(2008,01),end=c(2016,12),frequency=12)

#---------------------
# Is traf stationary??: We have 3 tools: plot of series, ACF and PACF.
#---------------------
plot.ts(traf2)
# we observe several patterns of non stationarity:
# increasing variance
# trend
# regular oscillations (seasonal effects)

par(mfrow=c(3,1))
plot.ts(traf2)
acf(ts(traf2, frequency=1), main="Autocorrelation main series")
pacf(ts(traf2, frequency=1), main="Partial autocorrelation main series")
par(mfrow=c(1,1))
# The ACF confirms non stationarity

#-----------------------------------------
# First transformation: log transformation
#-----------------------------------------
ltraf2 <- log(traf2)

plot(ltraf2)
acf(ts(ltraf2, frequency=1), 
    main="Autocorrelation log (main series)")
pacf(ts(ltraf2, frequency=1), 
     main="Partial autocorrelation log(main series)")
# No more increasing variance, but still trend and seasonality
# persistence of significant coeff in the ACF
# The log(series) is still non stat

#--------------------------------------------
# Second transformation: 1st order difference
#--------------------------------------------
dltraf2 <- diff(ltraf2, 1)
plot(dltraf2)
acf(ts(dltraf2, frequency=1),
    main="Autocorrelation 1st order diff of log (main series)")
pacf(ts(dltraf2, frequency=1), 
     main="Partial autocorrelation 1st order diff of log(main series)")

# We observe seasonality on the graph 
# And also on the ACF (no fast decay to 0 of the coefficients)
# the series is still non stationary

#---------------------------------------------
# Third transformation: difference of order 12
#---------------------------------------------
dltraf_12 <- diff(dltraf2, 12)
plot(dltraf_12)
acf(ts(dltraf_12, frequency=1), 
    main="Autocorrelation 1st order")
pacf(ts(dltraf_12, frequency=1), 
     main="Partial autocorrelation ")



#---------------------------------
# Identification of orders p and q
#---------------------------------
# d=D=1, s=12
# p=1, q=2
# P=0, Q=1 as a starting point
# Then try P=2 to identify effects at lag 24 (?)

#------------------------------------------------------------------------
# Estimation of a multiplicative SARIMA(2,1,0)(0,1,0) with seasonality 12 
#------------------------------------------------------------------------
mod1 <- arima(ltraf2, c(2,1,0), 
              seasonal=list(order=c(0,1,0), period=12), method='ML')
mod1


#AIC = -330.34

#plot of the fitted value
install.packages("TSA")
library("TSA")
fit1 <- fitted(mod1)
fit1
par(mfrow=c(1,1))
plot.ts(cbind(ltraf,fit1), plot.type='single', col=c('black','red'))

#--------------------
# Validation of mod3
#--------------------

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Significance of coefficients: pvalue of Student test on each coeff
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Student test to check for the significance of coeffs
#'''''''''''''''''''''''''''''''''''''''''''''''''''''
# Ho: coeff=0 against H1: coeff#0
# pvalue< 5%, we accept H1

mod1$coef  #value of coefficients
mod1$var.coef #variance of coefficients
tstat <- mod1$coef/sqrt(diag(mod1$var.coef)) # test statistic of the Student test
pvalue <- 2*(1-pnorm(abs(tstat))) #pvalue of the test 
tstat
pvalue
# the pvalues of AR1 and SAR1 are larger than 5%
# So we could remove these coefficients from the model

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Residuals analysis: assumption of Gaussian White Noise
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''
res1 <- mod1$residuals
plot(res1)
# We can identify one extreme value: the residuals are less than -0.2

# General overview of residuals
install.packages("forecast")
library("forecast")
checkresiduals(mod1)
boxplot(res1)

# Autocorrelation of the residuals
#'''''''''''''''''''''''''''''''''

acf(ts(res1,frequency=1))
# only one significant coeff at lag 6

# Box-Pierce test ou Ljung-Box test
#Ho: no autocorrelation against H1: autocorrelation

Box.test(res1, lag=20, type="Ljung-Box")
# Ho: all correlations up to lag 20 are =0 
# H1: at least one is different from 0
# pvalue=16% > 5%, we accept HO: no autocorrelation

# Normal Distribution assumption
#'''''''''''''''''''''''''''''''

res_norm <- res1/sqrt(mod1$sigma2) # normalized residuals
summary(res_norm)

# If the residuals follow a Normal distribution, the values of res_norm
# should lie in between -2 and 2, with 95% of chance

plot(res_norm)
abline(a=2,b=0,col="red")
abline(a=-2, b=0, col="red")

#QQplot 
# check whether the points are close to the line
qqnorm(res1)
qqline(res1)

#Shapiro test
shapiro.test(res1)
# pvalue > 5% so we accept the normality assumption

# Jarque-Bera test
install.packages('tseries')
library(tseries)
jarque.bera.test(res1)


# Constant variance assumption
#'''''''''''''''''''''''''''''
sq.res <- (res1)^2
acf(ts(sq.res, frequency=1))
# 2 highly significant coeffs: lag 1 and lag 12
# There is an issue of non constant variance

# TSA package to apply McLeod Li test
Htest <- McLeod.Li.test(mod1, plot=FALSE)
Htest
# The analysis of the 1st pvalue is enough
#pvalue < 5% so we reject the constant variance assumption

#------------
# Prediction
#-----------

# Check of the quality of the fit wrt confidence bounds
#-------------------------------------------------------

cb80 <- mod1$sigma2^.5*qnorm(0.9) # confidence bound of security 80%
plot(cbind(ltraf, fit1-cb80, fit1+cb80), plot.type='single', 
     lty=c(1,2,2), xlim=c(2008,2016))

# Proportion of points in the confidence interval
indi <- (ltraf-(fit1-cb80))>0&(fit1+cb80-ltraf)>0
prop <- 100*sum(indi)/length(indi)
prop
#prop = 80%

# if prop >= 80%, then the fit is considered good

# Validation set approach
#'''''''''''''''''''''''''

# 2008-2014 for training (~80%)
data.train <- window(traf2, start=c(2008,1), end=c(2014,12))

# 2015-2016 for testing (~20%)
data.test  <- window(traf2, start=c(2015,1), end=c(2016,12))


mod1.train <- arima(data.train, c(2,1,0), seasonal=list(order=c(0,1,0), period=12), method='ML')
pred1.test <- predict(mod1.train, n.ahead=174)


install.packages("forecast")
library("forecast")
accuracy(pred1.test$pred,data.test)

# ME     RMSE      MAE      MPE      MAPE       ACF1      Theil's U
#21298.7 25891.73 21659.08 7.149744 7.341169 0.3663249  0.477313

#The best model has the lowest values for the accuracy parameters

# plot comparing observed values and prediction of the traffic
ts.plot(traf, xlim=c(2008,2016), ylim=c(4e+05,12e+05))
ts.plot(traf, xlim=c(2008, 2018), col="black", main="Mauritius Air Traffic: Actual vs Forecast")
lines(pred1.test$pred, col="red")
lines(pred1.test$pred-1.96*pred1.test$se, col=4, lty=2)
lines(pred1.test$pred+1.96*pred1.test$se, col=4, lty=2)

# lty=line type



