######
###### Exploratory Data Analysis for time series data
######

library(astsa)

## 1. Monthly price of chicken (per pound)
tsplot(chicken, ylab = "cents per pound", col= "steelblue",
       lwd = 2, main = "Monthly price of chicken")
summary(fit <- lm(chicken~time(chicken))) # regress price on time
abline(fit)           # add the fitted regression line to the plot            

y1 <- residuals(fit)
y2 <- diff(chicken)

par(mfrow = c(2, 1), mar = c(3, 2, 1, 1))
tsplot(y1, ylab = "", main = "Detrended")
tsplot(y2, ylab = "", main = "Differenced")

par(mfrow = c(3, 1), mar = c(3, 4, 4, 4))
acf(chicken, lag.max = 36, main = "ACF of chicken data")
acf(y1, lag.max = 36, main = "ACF of detrended data")
acf(y2, lag.max = 36, main = "ACF of differenced data")

## 2. Global temperature
tsplot(globtemp, type = 'o', pch = 20, ylab = ~degree~C,
       main = "Global mean land-ocean temperature\nDeparture from 1951-1980 average")

## Linear filters: Moving averages
wgts15 <- c(-3, -6, -5, 3, 21, 46, 67, 74, 67, 46, 21, 3, -5, -6, -3) / 320
gtfilt15 <- filter(globtemp, wgts15)
lines(gtfilt15, col = "hotpink", lwd = 2)

## Nonparametric smoother: 'lowess'
lines(lowess(globtemp), lwd = 2, col = "tan")
lines(lowess(globtemp, f = 0.25), lwd = 2, col = "hotpink")

## 3. SOI
tsplot(soi, ylab = "", main = "Southern Oscillation Index")

wgts <- c(0.5, rep(1, 11), 0.5) / 12
soifilt <- filter(soi, wgts)
lines(soifilt, col = "steelblue", lwd = 2)

soifilt15 <- filter(soifilt, wgts15)
lines(soifilt15, col = "tan", lwd = 4)

tsplot(soi)
lines(lowess(soi, f = 0.05), lwd = 2, col = "orange")
lines(lowess(soi), lwd = 2, col = "violet")

