######
###### Maximum Likelihood Estimation for Dynamic Linear Models
######

library(dlm)
library(forecast)
library(astsa)

###
### Local linear trend
###

plot(globtemp, type = 'o', pch = 20, ylab = ~degree~C,
     main = "Global mean land-ocean temperature\nDeparture from 1951-1980 average")

mod0 <- dlmModPoly()
buildFun <- function(u) {
    mod0$V[] <- exp(u[1])
    diag(mod0$W) <- exp(u[2:3])
    mod0
}

fit <- dlmMLE(globtemp, c(-.5, .1, .6), buildFun)
fit$convergence # always check!! must be zero

modFit <- buildFun(fit$par)

## look at MLE of parameters
drop(V(modFit)) # observation variance
diag(W(modFit)) # system variance (diagonal)

## use the fitted model to forecast future temperature
modFilter <- dlmFilter(globtemp, modFit)
future_temp <- dlmForecast(modFilter, nAhead = 35)
plot(globtemp, type = 'o', pch = 20, ylab = ~degree~C, xlim = c(1880, 2050), ylim = c(min(globtemp), 2),
     main = "Global mean land-ocean temperature\nDeparture from 1951-1980 average")
lines(future_temp$f, col = "hotpink", lwd = 2, type = 'o', pch = 16, cex = 0.5)
## add 90% probability bands
lines(future_temp$f + sqrt(unlist(future_temp$Q)) * qnorm(0.05), col = "darkgreen")
lines(future_temp$f + sqrt(unlist(future_temp$Q)) * qnorm(0.95), col = "darkgreen")

###
### Local linear trend plus seasonal factors
###

wineind <- forecast::wineind / 1000
plot(wineind, type = 'o', pch = 20, ylab = expression(10^6 ~ "liters"),
     main = "Australian wine sales")
mod0 <- dlmModPoly() + dlmModSeas(12) # model 'skeleton'
buildFun <- function(u) {
    u <- exp(u)
    mod0$V[] <- u[1]
    diag(mod0$W)[1:3] <- u[2:4]
    mod0
}

fit <- dlmMLE(wineind, c(3, -1, 0, -.5), buildFun)
fit$conv
fit$value

## try random starting values, to make sure it is not a local minimum
set.seed(999)

A <- qr.Q(qr(matrix(rnorm(4^2), 4, 4))) # a random orthogonal matrix
crossprod(A) # sanity check - must be the identity matrix
A <- 3 * A
for (i in seq_len(ncol(A))) { 
    fit <- try(dlmMLE(wineind, parm = A[, i], buildFun), silent = TRUE)
    if (!inherits(fit, "try-error"))
        cat("i = ", i, ":\tconv = ", fit$conv, ",\tnegative LL = ", fit$value, "\n")
}

## use the last starting value, which seems to have converged to the global minimum
modFit <- buildFun(fit$par)

drop(V(modFit))
diag(W(modFit))

## use fitted model for forecasting
modFilter <- dlmFilter(wineind, modFit)
future_wine <- dlmForecast(modFilter, nAhead = 60)
plot(wineind, type = 'o', pch = 20, ylab = expression(10^6 ~~ "liters"),
     xlim = c(1980, 2000), ylim = c(7.9, 46.7), main = "Australian wine sales")
lines(future_wine$f, col = "hotpink", lwd = 2, type = 'o', pch = 16, cex = 0.5)
## add 90% probability bands
lines(future_wine$f + sqrt(unlist(future_wine$Q)) * qnorm(0.05), col = "darkgreen")
lines(future_wine$f + sqrt(unlist(future_wine$Q)) * qnorm(0.95), col = "darkgreen")

## smoothing
modSmooth <- dlmSmooth(modFilter)

plot(wineind, type = 'o', pch = 20, ylab = expression(10^6 ~~ "liters"),
     xlim = c(1980, 2000))# , ylim = c(7.9, 46.7), main = "Australian wine sales")
lines(modSmooth$s[, 1], col = "tomato", lwd = 2) # smoothed level
## recover posterior variance matrices, for probability bands
S <- with(modSmooth, dlmSvd2var(U.S, D.S))
level_var <- sapply(S, function(x) x[1, 1])
lines(modSmooth$s[, 1] + sqrt(level_var) * qnorm(0.05), col = "darkgreen")
lines(modSmooth$s[, 1] + sqrt(level_var) * qnorm(0.95), col = "darkgreen")
## add forecasted level ('deseasonalized')
lines(future_wine$a[, 1], col = "steelblue", lwd = 2, type = 'o', pch = 16, cex = 0.5)
future_level_var <- sapply(future_wine$R, function(x) x[1, 1])
lines(future_wine$a[, 1] + sqrt(future_level_var) * qnorm(0.05), col = "darkgreen")
lines(future_wine$a[, 1] + sqrt(future_level_var) * qnorm(0.95), col = "darkgreen")
title(main = "Deseasonalized wine sales, with forecasts for the next five years")

###
### CO2 level at Mauna Loa
###

tsplot(co2, type = 'o', pch = 20, cex = 0.5)

StructTS(co2)
mymod <- dlmModPoly(dW = c(0.11486, 0), dV = 0) +
    dlmModSeas(12, dW = c(0.09357, rep(0, 10)), dV = 1e-5)
filt <- dlmFilter(co2, mymod)
tsdiag(filt) # pretty bad...

# allow for extra variability in states
myBuild <- function(psi) {
    mod <- mymod
    psi <- exp(psi)
    mod$W[cbind(1:3, 1:3)] <- psi[1:3]
    mod$V[] <- psi[4]
    mod
}

fit <- dlmMLE(co2, rep(0.2, 4), myBuild)
fit
myFit <- myBuild(fit$par)
filt <- dlmFilter(co2, myFit)
tsdiag(filt) # still pretty bad...

## try with a 3rd order polynomial model
myBuild <- function(psi) {
    psi <- exp(psi)
    dlmModPoly(3, dV = psi[1], dW = c(0, 0, psi[2])) +
        dlmModSeas(12, dV = 0, dW = c(psi[3], rep(0, 10)))
}

fit <- dlmMLE(co2, rep(0.2, 3), myBuild)
fit
myFit <- myBuild(fit$par)
filt <- dlmFilter(co2, myFit)
tsdiag(filt) # nope...

## use trigonometric functions for the seasonal component,
##   and include a seasonal AR(1) to account for the residual
##   correlation at lag 12
myBuild <- function(psi) {
    psi[1:4] <- exp(psi[1:4])
    dlmModPoly(2, dW = c(psi[1:2]), dV = psi[3]) + 
        dlmModTrig(12, 3, dV = 0, dW = 0) +
        dlmModARMA(ar = c(rep(0, 11), ARtransPars(psi[5])), sigma2 = psi[4], dV = 0) 
    }

fit <- dlmMLE(co2, rep(-0.2, 5), myBuild)
fit
myFit <- myBuild(fit$par)
filt <- dlmFilter(co2, myFit)
tsdiag(filt) # finally! :-)






