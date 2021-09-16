library(tidyverse)
library(brotools)
library(lmtest)
library(anomalize)
library(timetk)
library(knitr)
library(tibbletime)
library(forecast)
library(zoo)
library(furrr)
library(purrr)
library(ggplot2)
library(tseries)
library(dplyr)
library(lubridate)
library(TSstudio)

#data handling
data1 <- read.csv("Documents/Uni Maths/Time series/TempMelbPRO.csv")
data1

head(data1)
dput(head(data1))
     
#change to monthly data
data2 <- data1 %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y"))
data2
dput(head(data2))

data_month <- data2 %>% group_by(month=floor_date(Date, "month")) %>%
  summarize(data=mean(Daily.maximum.temperatures.in.Melbourne..Australia..1981.1990))

data_month

data_month$data

#change to tim series data
y <- ts(data_month$data, frequency=12, start=c(1981,1))

#train/test plit
y_split <- ts_split(y, sample.out = 24)

y_train <- y_split$train
y_test <- y_split$test

y_test





#time series plot
ggplot(data_month, aes(month, data)) + 
  geom_line(color = "#6495ED", size = 1) + xlab("Time") + ylab("Temperature") +
  ggtitle("Maximum Tempuerature with 12 Month Moving Average") +
  geom_point(pch=20, fill="blue") +
  geom_line(aes(y=rollmean(data, 12, na.pad=TRUE)), size=1, color="red") +
  theme_minimal()

#decomposition plot
temptimeseriescomponents <- decompose(y)
plot(temptimeseriescomponents)

#stationarity
boxplot(y~cycle(y), ylab ='Temperature',
        xlab ='Time',col = 'darkgrey')

Adjusted_TS <- y - temptimeseriescomponents$seasonal
plot(Adjusted_TS)

Deseasonalised <- diff(y, lag= 12)
adf.test(diff(y, lag = 12))
# 0.1986
DesandDiff <- diff(Deseasonalised)
adf.test(diff(Deseasonalised), alternative="stationary", k=0)
#<0.01
ts.plot(Deseasonalised, col = 'darkgrey')
dev.new()

ts.plot(DesandDiff,col = 'darkgrey')
dev.new()

#ACF and PACF plots
acf(DesandDiff , lag.max= 85, main="ACF plot")
acf(DesandDiff, lag.max=85, plot=FALSE) # get the autocorrelation values
pacf(DesandDiff , lag.max= 85, main="PACF plot")
pacf(DesandDiff, lag.max=85, plot=FALSE) # get the autocorrelation values




#fitting models

#first sarima model
SARIMA1_y<-Arima(y_train, order=c(1,1,1), seasonal=c(1,1,1))
summary(SARIMA1_y)
checkresiduals(SARIMA1_y) 

#second sarima model
SARIMA2_y<-Arima(y_train, order=c(0,0,0), seasonal=c(2,1,0))
summary(SARIMA2_y)
checkresiduals(SARIMA2_y)

#ARIMA GRID SEARCH 1
forecast::auto.arima(y_train, d=1,D=1, ic="bic")
#function to loop through a range of ARIMA model to get the best order/seasonal based on BIC/AIC
#the output is ARIMA(0,1,1)(0,1,1)[12]

#ARIMA GRID SEARCH 2
SARIMA3_y <- auto.arima(y_train,d=1,D=1)
print(summary(SARIMA3_y))
#the output is ARIMA(0,1,2)(0,1,1)[12] 
#the second grid search gives us the same model, thus suggesting that this will fit the data best

#third sarima model
SARIMA3_y<-Arima(y_train, order=c(0,1,1), seasonal=c(0,1,1))
summary(SARIMA3_y)
checkresiduals(SARIMA3_y)

#third sarima model
SARIMA4_y<-Arima(y_train, order=c(0,1,2), seasonal=c(0,1,1))
summary(SARIMA4_y)
checkresiduals(SARIMA4_y)





#Forecast

plotForecastErrors <- function(forecasterrors)
{
  mybinsize <- IQR(forecasterrors)/4
  mysd   <- sd(forecasterrors)
  mymin  <- min(forecasterrors) - mysd*5
  mymax  <- max(forecasterrors) + mysd*3
  mynorm <- rnorm(10000, mean=0, sd=mysd)
  mymin2 <- min(mynorm)
  mymax2 <- max(mynorm)
  if (mymin2 < mymin) { mymin <- mymin2 }
  if (mymax2 > mymax) { mymax <- mymax2 }
  mybins <- seq(mymin, mymax, mybinsize)
  hist(forecasterrors, col="#6495ED", freq=FALSE, breaks=mybins)
  myhist <- hist(mynorm, plot=FALSE, breaks=mybins)
  points(myhist$mids, myhist$density, type="l", col="orange", lwd=2)
}

#model 1
forecast <- forecast(SARIMA1_y,h=24)
accuracy(forecast$mean, y_test)
plot(forecast,xlab = 'Time (years)', ylab = 'Maximum Temperature',
     shadecols = 'lavender', col = 'darkgrey')
forecast(SARIMA1_y)
plotForecastErrors(forecast$residuals)
acf(forecast$residuals, lag.max = 70, main= "ACF plot for ARIMA(1,1,1)(1,1,1)[12]")
#Ljung-box test
Box.test(forecast$residuals, lag=70, type="Ljung-Box")

#model 2
forecast1 <- forecast(SARIMA2_y,h=24)
accuracy(forecast1$mean, y_test)
plot(forecast1,xlab = 'Time (years)', ylab = 'Maximum Temperature',
     shadecols = 'lavender', col = 'darkgrey')
forecast(SARIMA2_y)
plotForecastErrors(forecast1$residuals)
acf(forecast1$residuals, lag.max = 70,main= "ACF plot for ARIMA(0,0,0)(2,1,0)[12]")
#Ljung-box test
Box.test(forecast1$residuals, lag=70, type="Ljung-Box")

#model 3
forecast2 <- forecast(SARIMA3_y,h=24)
accuracy(forecast2$mean, y_test)
plot(forecast(SARIMA3_y),xlab = 'Time (years)', ylab = 'Maximum Temperature',
     shadecols = 'lavender', col = 'darkgrey')
forecast(SARIMA3_y)
plotForecastErrors(forecast2$residuals)
acf(forecast2$residuals, lag.max = 70,main= "ACF plot for ARIMA(0,1,1)(0,1,1)[12]")
#Ljung-box test
Box.test(forecast2$residuals, lag=70, type="Ljung-Box")

#model 4
forecast3 <- forecast(SARIMA4_y,h=24)
accuracy(forecast3$mean, y_test)
plot(forecast(SARIMA4_y),xlab = 'Time (years)', ylab = 'Maximum Temperature',
     shadecols = 'lavender', col = 'darkgrey')
forecast(SARIMA4_y)
plotForecastErrors(forecast3$residuals)
acf(forecast3$residuals, lag.max = 70,main= "ACF plot for ARIMA(0,1,2)(0,1,1)[12]")
#Ljung-box test
Box.test(forecast3$residuals, lag=70, type="Ljung-Box")
