# Regression with ARIMA https://robjhyndman.com/hyndsight/forecasting-weekly-data/
# https://robjhyndman.com/hyndsight/longseasonality/
#https://stats.stackexchange.com/questions/177208/seasonal-adjustment-a-la-hyndman?noredirect=1&lq=1


#regression ---------------
# Load libraries
install.packages(c( "forecast", "Matrix", "lmtest", "tseries", "uroot"))
install.packages("CADFtest")

library(stats)
library(forecast)
library(Matrix)
library(lmtest)
library(tseries)
library(uroot)
library(CADFtest)
# Model --------------

# Y ----------------------------
names(rainpcdf)

## Time Series plot of PC1 -------------
adf.test(rainpcdf$PC1)
CADFtest(rainpcdf$PC1, max.lag.y = 10)


PC1ts <- ts(rainpcdf$PC1, frequency=12, start=c(1998,1))
library(ggplot2)
PC1ts=ggplot() + geom_line(data = rainpcdf, aes(x = date,y = PC1))
PC1ts=PC1ts + labs(y = "RPC1")+ labs(x = "Date")+ theme(legend.position='none')
PC1ts+ scale_x_date(date_breaks = "1 year", date_labels =  "%Y") 

res <- ch.test(PC1ts, type = "dummy", pvalue = "raw", sid = 12)
Box.test(PC1ts, lag=12)

d <- density(PC1ts) # returns the density data
plot(d, main="") 

hist(PC1ts)

plot(PC1ts, xaxt='n', las = 2,col = "blue", cex.axis = 0.78, cex.main=0.9,  xlab="Years", ylab="P(mm)")
axis(side = 1, at = pc.global.temp$date)

qqnorm(PC1ts, main="", col='red', cex.axis = 0.78)               
qqline(PC1ts) 

acf(PC1ts,  lag.max=24)
pacf(PC1ts, lag.max=24)



## Time Series plot of PC2---------------
adf.test(rainpcdf$PC2, k = 5)
CADFtest(rainpcdf$PC2, max.lag.y = 10)

PC2ts <- ts(rainpcdf$PC2, frequency=12, start=c(1998,1))

library(ggplot2)
PC2ts=ggplot() + geom_line(data = rainpcdf, aes(x = date,y = PC2))
PC2ts=PC2ts + labs(y = "RPC2")+ labs(x = "Date")+ theme(legend.position='none')
PC2ts+ scale_x_date(date_breaks = "1 year", date_labels =  "%Y") 

summary(PC2ts)
res <- ch.test(PC2ts, type = "dummy", pvalue = "raw", sid = 12)

d <- density(PC2ts) # returns the density data
plot(d, main="") 

hist(PC2ts)

plot(PC2ts, xaxt='n', las = 2,col = "blue", cex.axis = 0.78, cex.main=0.9,  xlab="Years", ylab="P(mm)")

qqnorm(PC2ts, main="", col='red', cex.axis = 0.78)               
qqline(PC2ts) 


summary(PC2ts)

d <- density(PC2ts) # returns the density data
plot(d, main="") 

hist(PC2ts)

plot(PC2ts, xaxt='n', las = 2,col = "blue", cex.axis = 0.78, cex.main=0.9,  xlab="Years", ylab="P(mm)")

qqnorm(PC2ts, main="", col='red', cex.axis = 0.78)               
qqline(PC2ts) 

## Time Series plot of PC3---------------
adf.test(rainpcdf$PC3, k = 3)
CADFtest(rainpcdf$PC3, max.lag.y = 12)


PC3ts <- ts(rainpcdf$PC3, frequency=12, start=c(1998,1))
library(ggplot2)
PC3ts=ggplot() + geom_line(data = rainpcdf, aes(x = date,y = PC3))
PC3ts=PC3ts + labs(y = "RPC3")+ labs(x = "Date")+ theme(legend.position='none')
PC3ts+ scale_x_date(date_breaks = "1 year", date_labels =  "%Y") 

summary(PC3ts)
res <- ch.test(PC3ts, type = "dummy", pvalue = "raw", sid = 12)

d <- density(PC3ts) # returns the density data
plot(d, main="") 

hist(PC3ts)

plot(PC3ts, xaxt='n', las = 3,col = "blue", cex.axis = 0.78, cex.main=0.9,  xlab="Years", ylab="P(mm)")

qqnorm(PC3ts, main="", col='red', cex.axis = 0.78)               
qqline(PC3ts) 


# X --------------------
attach(global_temp)
names(global_temp)



############################################################
# PC1 ARIMAX ---------------------------
############################################################
# Fourier series determine k------
library(forecast)
PC1ts <- ts(rainpcdf$PC1, frequency=12, start=c(1998,1))

bestfit <- list(aicc=Inf)
for(i in 1:6)
{
  fit <- auto.arima(PC1ts, xreg=fourier(PC1ts, K=i), seasonal=FALSE)
  if(fit$aicc < bestfit$aicc)
    bestfit <- fit
  else break;
}
fc <- forecast(bestfit, xreg=fourier(PC1ts, K=4, h=24))
plot(fc)
summary(bestfit)
coeftest(bestfit)
tsdisplay(residuals(bestfit))

# Choose independent variables --------------
require(forecast)
attach(global_temp)
# X with seasonal dummies------------
x1=cbind(rainy, short_rain, nino3.4anom, t, du031999, du102011, du052016)
x10=cbind(rainy, short_rain, nino3.4anom, t)

global_temp$du031999<- as.numeric(Year == 1999 & Month ==3)
global_temp$du102011<- as.numeric(Year == 2011 & Month ==10)
global_temp$du052016<- as.numeric(Year == 2016 & Month ==5)
global_temp$du042016<- as.numeric(Year == 2016 & Month ==4)

# Fitting autoARIMAX PC1 mean model-----------------
fitpc1 <- auto.arima(PC1ts,allowdrift=FALSE, seasonal=FALSE)
summary(fitpc1)
coeftest(fitpc1)
fitpc1 =arima(PC1ts, order=c(2,0,3))
summary(fitpc1)
coeftest(fitpc1)

# Fitting autoARIMAX PC1 X10 without fourier-----------------
fitpc1x10nof<- auto.arima(PC1ts, xreg= x10, allowdrift=FALSE, seasonal=FALSE)
summary(fitpc1x10nof)
coeftest(fitpc1x10nof)
fitpc1x10nof =arima(PC1ts, order=c(2,0,3), xreg= x10)
summary(fitpc1x10nof)
coeftest(fitpc1x10nof)


# Fitting autoARIMAX PC1 X1 without fourier-----------------
fitpc1x1nof<- auto.arima(PC1ts, xreg= x1, allowdrift=FALSE, seasonal=FALSE)
summary(fitpc1x1nof)
coeftest(fitpc1x1nof)
fitpc1x1nof =arima(PC1ts, order=c(2,0,3), xreg=x1)
summary(fitpc1x1nof)
coeftest(fitpc1x1nof)

# Fitting autoARIMAX PC1 X1 with one fourier-----------------
fitpc1x1f1 <- auto.arima(PC1ts, xreg=cbind(fourier(PC1ts, K=1), x1),allowdrift=FALSE, seasonal=FALSE)
summary(fitpc1x1f1)
coeftest(fitpc1x1f1)
fitpc1x1f1 =arima(PC1ts, order=c(2,0,3), xreg=cbind(fourier(PC1ts, K=1), x1))
summary(fitpc1x1f1)
coeftest(fitpc1x1f1)


# Fitting autoARIMAX PC1 X1 with two fourier-----------------
fitpc1x1f2 <- auto.arima(PC1ts, xreg=cbind(fourier(PC1ts, K=2), x1),allowdrift=FALSE, seasonal=FALSE)
summary(fitpc1x1f2)
coeftest(fitpc1x1f2)
res1=residuals(fitpc1x1f2)
tsdisplay(residuals(fitpc1x1f2))
Box.test(residuals(fitpc1x1f2), lag=24, fitdf=3, type="Box-Pierce")
adf.test(res1, k = 6)
CADFtest(res1, max.lag.y = 12)

# Fitting ARIMAX PC1 with two fourier-----------------
fitpc1x1f2 =arima(PC1ts, order=c(2,0,3), xreg=cbind(fourier(PC1ts, K=2), x1))
summary(fitpc1x1f2)
# test significance of coefficients
coeftest(fitpc1x1f2)
#checking residuals
tsdisplay(residuals(fitpc1x1f2))
Box.test(residuals(fitpc1x1f2), lag=24, fitdf=1, type="Box-Pierce")


############################################################
# PC2 ARIMAX ---------------------------
###################################################################

# Fourier series determine k------
library(forecast)
PC2ts <- ts(rainpcdf$PC2, frequency=12, start=c(1998,1))

bestfit <- list(aicc=Inf)
for(i in 1:6)
{
  fit <- auto.arima(PC2ts, xreg=fourier(PC2ts, K=i), seasonal=FALSE)
  if(fit$aicc < bestfit$aicc)
    bestfit <- fit
  else break;
}
fc <- forecast(bestfit, xreg=fourier(PC2ts, K=4, h=24))
plot(fc)
summary(bestfit)
coeftest(bestfit)
tsdisplay(residuals(bestfit))

# Choose independent variables --------------
require(forecast)
attach(global_temp)
# X with seasonal dummies------------
x2=cbind(amm, soi, du072015)
x20=cbind(amm, soi)
global_temp$du072015<- as.numeric(Year == 2015 & Month ==7)
global_temp$du042016<- as.numeric(Year == 2016 & Month ==4)

# Fitting autoARIMAX pc2 mean model-----------------
fitpc2 <- auto.arima(PC2ts,allowdrift=FALSE, seasonal=FALSE)
summary(fitpc2)
coeftest(fitpc2)
fitpc2 =arima(PC2ts, order=c(1,0,0))
summary(fitpc2)
coeftest(fitpc2)

# Fitting autoARIMAX pc2 x20 without fourier-----------------
fitpc2x20nof<- auto.arima(PC2ts, xreg= x20, allowdrift=FALSE, seasonal=FALSE)
summary(fitpc2x20nof)
coeftest(fitpc2x20nof)
fitpc2x20nof =arima(PC2ts, order=c(1,0,0), xreg= x20)
summary(fitpc2x20nof)
coeftest(fitpc2x20nof)


# Fitting autoARIMAX pc2 x2 without fourier-----------------
fitpc2x2nof<- auto.arima(PC2ts, xreg= x2, allowdrift=FALSE, seasonal=FALSE)
summary(fitpc2x2nof)
coeftest(fitpc2x2nof)
fitpc2x2nof =arima(PC2ts, order=c(1,0,0), xreg=x2)
summary(fitpc2x2nof)
coeftest(fitpc2x2nof)

# Fitting autoARIMAX pc2 x2 with one fourier-----------------
fitpc2x2f1 <- auto.arima(PC2ts, xreg=cbind(fourier(PC2ts, K=1), x2),allowdrift=FALSE, seasonal=FALSE)
summary(fitpc2x2f1)
coeftest(fitpc2x2f1)
fitpc2x2f1 =arima(PC2ts, order=c(1,0,0), xreg=cbind(fourier(PC2ts, K=1), x2))
summary(fitpc2x2f1)
coeftest(fitpc2x2f1)


# Fitting autoARIMAX pc2 x2 with two fourier-----------------
fitpc2x2f2 <- auto.arima(PC2ts, xreg=cbind(fourier(PC2ts, K=2), x2),allowdrift=FALSE, seasonal=FALSE)
summary(fitpc2x2f2)
coeftest(fitpc2x2f2)
res1=residuals(fitpc2x2f2)
tsdisplay(residuals(fitpc2x2f2))
Box.test(residuals(fitpc2x2f2), lag=24, fitdf=3, type="Box-Pierce")
adf.test(res1, k = 6)
CADFtest(res1, max.lag.y = 12)

# Fitting ARIMAX pc2 with two fourier-----------------
fitpc2x2f2 =arima(PC2ts, order=c(1,0,0), xreg=cbind(fourier(PC2ts, K=2), x2))
summary(fitpc2x2f2)
# test significance of coefficients
coeftest(fitpc2x2f2)
#checking residuals
tsdisplay(residuals(fitpc2x2f2))
Box.test(residuals(fitpc2x2f2), lag=24, fitdf=1, type="Box-Pierce")

# Fitting autoARIMAX pc2 x2 with three fourier-----------------
fitpc2x2f3 <- auto.arima(PC2ts, xreg=cbind(fourier(PC2ts, K=3), x2),allowdrift=FALSE, seasonal=FALSE)
summary(fitpc2x2f3)
coeftest(fitpc2x2f3)
res1=residuals(fitpc2x2f3)
tsdisplay(residuals(fitpc2x2f3))
Box.test(residuals(fitpc2x2f3), lag=24, fitdf=3, type="Box-Pierce")
adf.test(res1, k = 6)
CADFtest(res1, max.lag.y = 12)

# Fitting ARIMAX pc2 with three fourier-----------------
fitpc2x2f3 =arima(PC2ts, order=c(1,0,0), xreg=cbind(fourier(PC2ts, K=3), x2))
summary(fitpc2x2f3)
# test significance of coefficients
coeftest(fitpc2x2f3)
#checking residuals
tsdisplay(residuals(fitpc2x2f3))
Box.test(residuals(fitpc2x2f3), lag=24, fitdf=1, type="Box-Pierce")


############################################################
# PC3 ARIMAX ---------------------------


############################################################
# Fourier series determine k------
library(forecast)
PC3ts <- ts(rainpcdf$PC3, frequency=12, start=c(1998,1))

bestfit <- list(aicc=Inf)
for(i in 1:6)
{
  fit <- auto.arima(PC3ts, xreg=fourier(PC3ts, K=i), seasonal=FALSE)
  if(fit$aicc < bestfit$aicc)
    bestfit <- fit
  else break;
}
fc <- forecast(bestfit, xreg=fourier(PC3ts, K=4, h=24))
plot(fc)
summary(bestfit)
coeftest(bestfit)
tsdisplay(residuals(bestfit))

# Choose independent variables --------------
require(forecast)
attach(global_temp)
# X with seasonal dummies------------
x3=cbind(nino3.4anom, du052003, du092011, du092012)
x30=cbind(nino3.4anom)

global_temp$du052003<- as.numeric(Year == 2003 & Month ==5)
global_temp$du092011<- as.numeric(Year == 2011 & Month ==9)
global_temp$du092012<- as.numeric(Year == 2012 & Month ==9)

# Fitting autoARIMAX pc3 mean model-----------------
fitpc3 <- auto.arima(PC3ts,allowdrift=FALSE, seasonal=FALSE)
summary(fitpc3)
coeftest(fitpc3)
fitpc3 =arima(PC3ts, order=c(0,0,0))
summary(fitpc3)
coeftest(fitpc3)

# Fitting autoARIMAX pc3 x30 without fourier-----------------
fitpc3x30nof<- auto.arima(PC3ts, xreg= x30, allowdrift=FALSE, seasonal=FALSE)
summary(fitpc3x30nof)
coeftest(fitpc3x30nof)
fitpc3x30nof =arima(PC3ts, order=c(0,0,0), xreg= x30)
summary(fitpc3x30nof)
coeftest(fitpc3x30nof)


# Fitting autoARIMAX pc3 x3 without fourier-----------------
fitpc3x3nof<- auto.arima(PC3ts, xreg= x3, allowdrift=FALSE, seasonal=FALSE)
summary(fitpc3x3nof)
coeftest(fitpc3x3nof)
fitpc3x3nof =arima(PC3ts, order=c(0,0,0), xreg=x3)
summary(fitpc3x3nof)
coeftest(fitpc3x3nof)

# Fitting autoARIMAX pc3 x3 with one fourier-----------------
fitpc3x3f1 <- auto.arima(PC3ts, xreg=cbind(fourier(PC3ts, K=1), x3),allowdrift=FALSE, seasonal=FALSE)
summary(fitpc3x3f1)
coeftest(fitpc3x3f1)
fitpc3x3f1 =arima(PC3ts, order=c(0,0,0), xreg=cbind(fourier(PC3ts, K=1), x3))
summary(fitpc3x3f1)
coeftest(fitpc3x3f1)


# Fitting autoARIMAX pc3 x3 with two fourier-----------------
fitpc3x3f2 <- auto.arima(PC3ts, xreg=cbind(fourier(PC3ts, K=2), x3),allowdrift=FALSE, seasonal=FALSE)
summary(fitpc3x3f2)
coeftest(fitpc3x3f2)
res1=residuals(fitpc3x3f2)
tsdisplay(residuals(fitpc3x3f2))
Box.test(residuals(fitpc3x3f2), lag=24, fitdf=3, type="Box-Pierce")
adf.test(res1, k = 6)
CADFtest(res1, max.lag.y = 12)

# Fitting ARIMAX pc3 with two fourier-----------------
fitpc3x3f2 =arima(PC3ts, order=c(0,0,0), xreg=cbind(fourier(PC3ts, K=2), x3))
summary(fitpc3x3f2)
# test significance of coefficients
coeftest(fitpc3x3f2)
#checking residuals
tsdisplay(residuals(fitpc3x3f2))
Box.test(residuals(fitpc3x3f2), lag=24, fitdf=1, type="Box-Pierce")

# Fitting autoARIMAX pc3 x3 with three fourier-----------------
fitpc3x3f3 <- auto.arima(PC3ts, xreg=cbind(fourier(PC3ts, K=3), x3),allowdrift=FALSE, seasonal=FALSE)
summary(fitpc3x3f3)
coeftest(fitpc3x3f3)
res1=residuals(fitpc3x3f3)
tsdisplay(residuals(fitpc3x3f3))
Box.test(residuals(fitpc3x3f3), lag=24, fitdf=3, type="Box-Pierce")
adf.test(res1, k = 6)
CADFtest(res1, max.lag.y = 12)

# Fitting ARIMAX pc3 with three fourier-----------------
fitpc3x3f3 =arima(PC3ts, order=c(0,0,0), xreg=cbind(fourier(PC3ts, K=3), x3))
summary(fitpc3x3f3)
# test significance of coefficients
coeftest(fitpc3x3f3)
#checking residuals
tsdisplay(residuals(fitpc3x3f3))
Box.test(residuals(fitpc3x3f3), lag=24, fitdf=1, type="Box-Pierce")


# Export results to latex----------------------
library(stargazer)
stargazer(fitpc1,fitpc1x10nof, fitpc1x1nof, fitpc1x1f1, fitpc1x1f2)
stargazer(fitpc2,fitpc2x20nof, fitpc2x2nof, fitpc2x2f1, fitpc2x2f2, fitpc2x2f3)
stargazer(fitpc3,fitpc3x30nof, fitpc3x3nof, fitpc3x3f1, fitpc3x3f2, fitpc3x3f3)
lrtest (fitpc1x10nof, fitpc1)
lrtest (fitpc1x1f2, fitpc1x1f1,fitpc1x1nof, fitpc1x10nof, fitpc1)
lrtest (fitpc2x2f3, fitpc2x2f2, fitpc2x2f1,fitpc2x2nof, fitpc2x20nof, fitpc2)
lrtest (fitpc3x3f3, fitpc3x3f2, fitpc3x3f1,fitpc3x3nof, fitpc3x30nof, fitpc3)
