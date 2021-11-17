######
###### Introduction to ARMA models
######

library(astsa)

### Causality and Invertibility of ARMA(p,q) models
z <- c(1, -1.5, 0.75)     ## coefficients of the polynomial
(a <- polyroot(z))        ## roots of polynomial
Mod(a)                    ## all outside the unit circle 
Arg(a)                    ## the argument of the roots, in radians
arg <- Arg(a)[1] / (2*pi)

1 / arg                   ## the period

ARMAtoMA(ar = c(1.5, -0.75), lag.max = 10)   ## first 10 psi-weights

### ACF/PACF of ARMA(p, q) models
ACF <- ARMAacf(ar = c(1.5,-0.75), lag.max = 36)
plot(ACF, type = "h", xlab = "lag")
abline(h = 0)

PACF <- ARMAacf(ar = c(1.5,-0.75), lag.max = 12, pacf = TRUE)
plot(PACF, type = "h", xlab = "lag")
abline(h = 0)

### Empirical ACF/PACF and model identification
set.seed(314)
ar2 <- arima.sim(list(order = c(2, 0, 0), ar = c(1.5, -0.75)), n = 144)
plot(ar2, axes = FALSE, xlab = "Time")
axis(2); axis(1, at = seq(0, 144, by = 12)); box()  # work the plot machine
abline(v = seq(0, 144, by = 12), lty = 2)

acf(ar2)
pacf(ar2)

### Invertibility
z <- c(1, -1, 0.25)
polyroot(z)

ACF <- ARMAacf(ma = c(-1, 0.25), lag.max = 10)
plot(ACF, type = "h", xlab = "lag")
abline(h = 0)

PACF <- ARMAacf(ma = c(-1,0.25), lag.max = 10, pacf = TRUE)
plot(PACF, type = "h", xlab = "lag")
abline(h = 0)

set.seed(2)
ma2 <- arima.sim(list(order = c(0, 0, 2), ma = c(-1, 0.25)), n = 120)

acf(ma2)
pacf(ma2)
