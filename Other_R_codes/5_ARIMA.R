#########
######### ARIMA models for nonstationary data
#########

library(astsa)

### 1. Glacial varve series
tsplot(varve, ylab = "", main = "Glacial varve")
tsplot(log(varve), ylab = "", main = "Glacial varve")
acf(log(varve), lag.max = 50)
tsplot(diff(log(varve)), ylab = "", main = "Glacial varve")
acf(diff(log(varve)), lag.max = 50)
pacf(diff(log(varve)), lag.max = 50)

fit0 <- arima(log(varve), order = c(0, 1, 1))
fit0
tsdiag(fit0)

fit1 <- arima(log(varve), order = c(1, 1, 1))
fit1
tsdiag(fit1)

### 2. A simulated pure seasonal ARMA model
set.seed(555)
phi <- 0.75
theta <- -0.5
x <- arima.sim(model = list(ar = c(rep(0, 11), phi)), n = 200)
tsplot(x)
acf(x, lag.max = 50)
pacf(x, lag.max = 50)

ACF <- ARMAacf(ar = c(rep(0, 11), phi), lag.max = 50)
PACF <- ARMAacf(ar = c(rep(0, 11), phi), lag.max = 50, pacf = TRUE)
plot(ACF, type = "h", xlab = "lag")
abline(h = 0)
plot(PACF, type = "h", xlab = "lag")
abline(h = 0)

### 3. Price of chicken
tsplot(chicken, ylab = "", main = "Price of chicken")
chicken_lm <- lm(chicken ~ time(chicken))
abline(chicken_lm)

## Model identification
## Differencing - instead of fitting a regression line
par(mfrow = c(2, 1))
acf(diff(chicken), lag.max = 48)
pacf(diff(chicken), lag.max = 48)
acf(diff(chicken, lag = 12), lag.max = 24)
pacf(diff(chicken, lag = 12), lag.max = 24)
acf(diff(diff(chicken, lag = 12)), lag.max = 37)
pacf(diff(diff(chicken, lag = 12)), lag.max = 37)


fit0 <- arima(chicken, order = c(1, 1, 0), seasonal = list(order = c(0, 1, 0), period = 12))
fit0
tsdiag(fit0)

fit1 <- arima(chicken, order = c(1, 1, 0), seasonal = list(order = c(0, 1, 1), period = 12))
fit1
tsdiag(fit1)

fit2 <- arima(chicken, order = c(1, 1, 1), seasonal = list(order = c(1, 1, 0), period = 12))
fit2
tsdiag(fit2)

fit3 <- arima(chicken, order = c(1, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12))
fit3
tsdiag(fit3)

## Forecasting
chicken_fore <- predict(fit3, n.ahead = 30)
U <- chicken_fore$pred + chicken_fore$se
L <- chicken_fore$pred - chicken_fore$se
tsplot(chicken, xlim = c(2001.5, 2020), ylim = c(min(chicken), max(U)),
       ylab = "", main = "Chicken price")
lines(chicken_fore$pred, type = "o", pch = 20, col = 2)
xx <- c(time(U), rev(time(U))); yy <- c(L, rev(U))
polygon(xx, yy, border = 8, col = gray(0.6, alpha = 0.2))
lines(U, col = 4); lines(L, col = 4)


### 4. Airline passengers
x <- AirPassengers
tsplot(x, ylab = "", main = "Air Passengers")
frequency(x) # monthly data
lx <- log(x); dlx <- diff(lx); ddlx <- diff(dlx, 12)
plot(cbind(x, lx, dlx, ddlx), main = "")
monthplot(dlx)
monthplot(ddlx)
acf(ddlx, lag.max = 50)
pacf(ddlx, lag.max = 50)

fit0 <- arima(lx, order = c(1, 1, 1), seasonal = c(0, 1, 1)) 
fit0
tsdiag(fit0)

fit1 <- arima(lx, order = c(0, 1, 1), seasonal = c(0, 1, 1)) 
fit1
tsdiag(fit1)

## Forecasting
air_fore <- predict(fit1, n.ahead = 25)
U <- air_fore$pred + air_fore$se
L <- air_fore$pred - air_fore$se
tsplot(window(lx, start = 1957), xlim = c(1957, 1963),
       ylim = c(min(window(lx, start = 1957)), max(U)),
       ylab = "", main = "Log Air Passengers", type = 'o', pch = 16)
lines(air_fore$pred, type = "o", pch = 20, col = 2)
xx <- c(time(U), rev(time(U))); yy <- c(L, rev(U))
polygon(xx, yy, border = 8, col = gray(0.6, alpha = 0.2))
lines(U, col = 4); lines(L, col = 4)

## Back to original scale
U <- exp(U); L <- exp(L)
tsplot(window(x, start = 1957), xlim = c(1957, 1963),
       ylim = c(min(window(x, start = 1957)), max(U)),
       ylab = "", main = "Air Passengers", type = 'o', pch = 16)
lines(exp(air_fore$pred), type = "o", pch = 20, col = 2)
xx <- c(time(U), rev(time(U))); yy <- c(L, rev(U))
polygon(xx, yy, border = 8, col = gray(0.6, alpha = 0.2))
lines(U, col = 4); lines(L, col = 4)

### 5. SOI and Recruitment
par(mfrow = c(2, 1)) 
length(soi)
tsplot(soi, ylab = "", main = "Southern Oscillation Index")
tsplot(rec, ylab = "", main = "Recruitment") 

lag1.plot(soi, 12)
dev.new()
lag2.plot(soi, rec, 8)

## Note: this could benefit from a seasonal model fit, but it hasn't
##  been talked about yet - you could come back to this after the next section
dummy <- ifelse(soi < 0, 0, 1)
fish <- ts.intersect(rec, soiL6 = lag(soi, -6), dL6 = lag(dummy, -6), dframe = TRUE)
summary(fit <- lm(rec ~ soiL6 * dL6, data = fish, na.action = NULL))
attach(fish)
plot(resid(fit))
acf(resid(fit))     # indicates AR(2)
intract <- soiL6 * dL6  # interaction term
sarima(rec,2,0,0, xreg = cbind(soiL6, dL6, intract))
# sarima(rec,2,0,0,0,1,1,12, xreg = cbind(soiL6, dL6, intract))
