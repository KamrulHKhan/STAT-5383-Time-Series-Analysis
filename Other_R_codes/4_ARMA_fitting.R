######
###### Fitting ARMA models
######

library(astsa)

### 1. Recruitment series
tsplot(rec, ylab = "", main = "Recruitment")

## Model identification
par(mfrow = c(2, 1))
acf(rec, lag.max = 37)
pacf(rec, lag.max = 37)

## Model fitting
recAR2_fit <- arima(rec, order = c(2, 0, 0))
recAR2_fit

## Implications
z <- c(1, -coef(recAR2_fit)[c("ar1", "ar2")])
(a <- polyroot(z))
Mod(a) ## outside the unit disk, as it should for causality
2 * pi / Arg(a)[1]  ## (pseudo) period, in months

## Forecasting
rec_fore <- predict(recAR2_fit, n.ahead = 24)
tsplot(rec, xlim = c(1980, 1990), ylab = "", main = "Recruitment")
lines(rec_fore$pred, type = "o", pch = 20, col = 2)
U <- rec_fore$pred + rec_fore$se
L <- rec_fore$pred - rec_fore$se
lines(U, col = 4); lines(L, col = 4)
xx <- c(time(U), rev(time(U))); yy <- c(L, rev(U))
polygon(xx, yy, border = 8, col = gray(0.6, alpha = 0.2))
abline(h = coef(recAR2_fit)["intercept"], lty = "dashed", col = grey(0.8))

### 2. Price of chicken
tsplot(chicken, ylab = "", main = "Price of chicken")
chicken_lm <- lm(chicken ~ time(chicken))
abline(chicken_lm)

## Model identification
par(mfrow = c(2, 1))
acf(residuals(chicken_lm), lag.max = 37)
pacf(residuals(chicken_lm), lag.max = 37)

## Model fitting
chicken_fit <- arima(chicken, order = c(2, 0, 0), xreg = time(chicken))
chicken_fit

## Implications
z <- c(1, -coef(chicken_fit)[c("ar1", "ar2")])
(a <- polyroot(z))
2 * pi / Arg(a)[1]  ## (pseudo) period, in months

## Forecasting
new_time <- data.frame(seq(tail(time(chicken), n = 1), length = 30, by = 1/12) + 1/12)
colnames(new_time) <- "time(chicken)"
chicken_fore <- predict(chicken_fit, n.ahead = 30, newxreg = new_time)
U <- chicken_fore$pred + chicken_fore$se
L <- chicken_fore$pred - chicken_fore$se
tsplot(chicken, xlim = c(2001.5, 2020), ylim = c(min(chicken), max(U)),
       ylab = "", main = "Chicken price")
lines(chicken_fore$pred, type = "o", pch = 20, col = 2)
lines(U, col = 4); lines(L, col = 4)
xx <- c(time(U), rev(time(U))); yy <- c(L, rev(U))
polygon(xx, yy, border = 8, col = gray(0.6, alpha = 0.2))
abline(coef(chicken_fit)[c(3, 4)], lty = "dashed", col = grey(0.6))

### 3. Glacial varve
tsplot(varve, ylab = "", main = "Glacial varve")
tsplot(log(varve), ylab = "", main = "Glacial varve")
acf(log(varve), lag.max = 50)
tsplot(diff(log(varve)), ylab = "", main = "Glacial varve")
acf(diff(log(varve)), lag.max = 50)

x <- diff(log(varve))

## Model identification
par(mfrow = c(2, 1))
acf(x)
pacf(x)

## Model fitting
fit <- arima(x, order = c(0, 0, 1))
fit

### 4. US GNP, quarterly growth rate
tsplot(gnp, ylab = "", main = "US GNP")
y <- diff(log(gnp)) # growth rate
tsplot(y, ylab = "", main = "GNP growth rate")

## Model identification
par(mfrow = c(2, 1))
acf(y, lag.max = 20)
pacf(y, lag.max = 20)

## Model fitting
fitAR <- arima(y, order = c(1, 0, 0))
fitAR
fitMA <- arima(y, order = c(0, 0, 2))
fitMA

## Diagnostics
tsdiag(fitAR)
tsdiag(fitMA)

round(ARMAtoMA(ar = 0.35, ma = 0, 10), 3)
