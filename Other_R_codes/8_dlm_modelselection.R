######
###### Model Selection Criteria for Dynamic Linear Models
######

library(dlm)
library(astsa)
# library(forecast)


######
###### Nile river flow
######

#### Local level model
mod_level <- dlmModPoly(1)
build_level <- function(u) {
    u <- exp(u)
    mod_level$V[] <- u[1]
    mod_level$W[] <- u[2]
    mod_level
}

fit_level <- dlmMLE(Nile, c(.3, .7), build_level)
fit_level$conv
mod_level <- build_level(fit_level$par)

## AIC
2 * fit_level$value + 2 * length(fit_level$par)

## BIC
2 * fit_level$value + log(length(Nile)) * length(fit_level$par)

#### Local linear trend model
mod_trend <- dlmModPoly(2)
build_trend <- function(u) {
    u <- exp(u)
    mod_trend$V[] <- u[1]
    diag(mod_trend$W) <- u[2 : 3]
    mod_trend
}

## try random starting values, to make sure it is not a local minimum
set.seed(999)
n_start <- 10
out <- vector("list", n_start)
start <- matrix(rnorm(3 * n_start, sd = 4), nrow = 3) 

for (i in seq_len(n_start)) {
    fit <- try(dlmMLE(Nile, start[, i], build_trend), silent = TRUE)
    if (!inherits(fit, "try-error")) {
        cat("i = ", i, ":\tconv = ", fit$conv, ",\tnegative LL = ", fit$value, "\n")
    } else 
        cat("i = ", i, ":\tError\n")
    out[[i]] <- fit
}

goodOut <- out[sapply(out, function(x) !inherits(x, "try-error") && x[["convergence"]] == 0)]
if (length(goodOut))
    jj <- which.min(sapply(goodOut, function(x) x[["value"]]))

fit_trend <- goodOut[[jj]]
mod_trend <- build_trend(fit_trend$par)


#### Integrated random walk (IRW)
mod_IRW <- dlmModPoly(2)
build_IRW <- function(u) {
    u <- exp(u)
    mod_IRW$V[] <- u[1]
    diag(mod_IRW$W)[2] <- u[2]
    mod_IRW
}

## try random starting values, to make sure it is not a local minimum
set.seed(999)
n_start <- 10
out <- vector("list", n_start)
start <- matrix(rnorm(3 * n_start, sd = 4), nrow = 3) 

for (i in seq_len(n_start)) {
    fit <- try(dlmMLE(Nile, start[, i], build_IRW), silent = TRUE)
    if (!inherits(fit, "try-error")) {
        cat("i = ", i, ":\tconv = ", fit$conv, ",\tnegative LL = ", fit$value, "\n")
    } else 
        cat("i = ", i, ":\tError\n")
    out[[i]] <- fit
}

goodOut <- out[sapply(out, function(x) !inherits(x, "try-error") && x[["convergence"]] == 0)]
if (length(goodOut))
    jj <- which.min(sapply(goodOut, function(x) x[["value"]]))

fit_IRW <- goodOut[[jj]]
mod_IRW <- build_IRW(fit_IRW$par)

### Compare AIC/BIC
## AIC
2 * fit_level$value + 2 * length(fit_level$par)
2 * fit_trend$value + 2 * length(fit_trend$par)
2 * fit_IRW$value + 2 * length(fit_IRW$par)

## BIC
2 * fit_level$value + log(length(Nile)) * length(fit_level$par)
2 * fit_trend$value + log(length(Nile)) * length(fit_trend$par)
2 * fit_IRW$value + log(length(Nile)) * length(fit_IRW$par)

######
###### Nottingham monthly temperature
######

### Fit Fourier seasonal models
tsplot(nottem, type = 'o', pch = 20, cex = 0.75, ylab = "F", xlab = "",
       main = "Average monthly temperature in Nottingham (UK)")

modFourier <- vector(mode = "list", length = 6) # for the models
buildFourier <- vector(mode = "list", length = 6) # for the 'build' functions
for (i in seq_along(modFourier)) {
    modFourier[[i]] <- dlmModPoly(1) + dlmModTrig(s = 12, q = i) # model 'skeleton'
    buildFourier[[i]] <- function(u) {
        u <- exp(u)
        modFourier[[i]]$V[] <- u[1]
        diag(modFourier[[i]]$W)[1] <- u[2]
        diag(modFourier[[i]]$W)[-1] <- u[3]
        modFourier[[i]]
    }
    fit <- try(dlmMLE(nottem, c(3, -1, -.5), buildFourier[[i]]))
    if (!inherits(fit, "try-error")) {
        cat("i = ", i, ":\tconv = ", fit$conv, ",\tnegative LL = ", fit$value, "\n")
        modFourier[[i]] <- buildFourier[[i]](fit$par)
    }
}

for (i in seq_len(6)) 
    cat("i = ", i, ":\tMAE = ", mean(abs(tail(residuals(dlmFilter(nottem, modFourier[[i]]),
                                                        type = "raw", sd = FALSE), n = -24))), "\n")

for (i in seq_len(6)) 
    cat("i = ", i, ":\tRMSE = ", sqrt(mean(tail(residuals(dlmFilter(nottem, modFourier[[i]]),
                                                          type = "raw", sd = FALSE), n = -24)^2)), "\n")

for (i in seq_len(6)) 
    cat("i = ", i, ":\tMAPE = ", mean(tail(abs(residuals(dlmFilter(nottem, modFourier[[i]]),
                                                         type = "raw", sd = FALSE)) / nottem, n = -24)), "\n")

## best model, with 2 harmonics
i <- 2
sqrt(diag(W(modFourier[[i]]))) # standard deviations of state noises
sqrt(drop(V(modFourier[[i]]))) # standard deviations of observation error

### Use the best model to forecast temperature next 5 years
modFilter <- dlmFilter(nottem, modFourier[[i]])
future_temp <- dlmForecast(modFilter, nAhead = 60)

tsplot(window(nottem, start = 1930), type = 'o', pch = 20, cex = 0.75,
       ylab = "F", xlab = "",
       main = "Average monthly temperature in Nottingham (UK)",
       xlim = c(1929.9, 1945.08))
## point forecasts
lines(future_temp$f, col = "hotpink", lwd = 2, type = 'o', pch = 16, cex = 0.75)
## add 90% probability bands
lines(future_temp$f + sqrt(unlist(future_temp$Q)) * qnorm(0.05), col = "darkgreen")
lines(future_temp$f + sqrt(unlist(future_temp$Q)) * qnorm(0.95), col = "darkgreen")

######
###### UK quarterly gas consumption 
######

gas <- log(UKgas)
tsplot(gas, type = 'o', pch = 16, cex = 0.5, ylab = "", xlab = "",
       main = "UK quarterly gas consumption (log millions of therms)")

######
###### Annual number of sunspots
######

tsplot(sqrt(sunspot.year), type = 'o', xlab = "", ylab = "",
       main = "Annual number of sunspots (square root)", pch = 20, cex = 0.5)
spots <- sqrt(sunspot.year)

spotBuild <- function(psi) {
    psi[1:3] <- exp(psi[1:3])
    my_mod <- dlmModPoly(1) + dlmModTrig(q = 2, tau = psi[4])
    my_mod$V[] <- psi[1]
    diag(my_mod$W) <- c(psi[2], rep(psi[3], 2 * 2))
    my_mod
}

fit <- dlmMLE(spots, c(1, 0.2, 0.1, 8), spotBuild,
              lower = c(rep(-Inf, 3), 5), upper = c(rep(Inf, 3), 15))
fit$convergence
fit$par # period tau = 10.817

### Model selection
h_max <- 8 # maximum number of harmonics to be considered
modFourier <- vector(mode = "list", length = h_max) # for the models
buildFourier <- vector(mode = "list", length = h_max) # for the 'build' functions
for (i in seq_along(modFourier)) {
    buildFourier[[i]] <- function(psi) {
        psi[1:3] <- exp(psi[1:3])
        my_mod <- dlmModPoly(1) + dlmModTrig(q = i, tau = psi[4])
        my_mod$V[] <- psi[1]
        diag(my_mod$W) <- c(psi[2], rep(psi[3], 2 * i))
        my_mod
    }
    fit <- try(dlmMLE(spots, c(1, 0.2, 0.1, 8), buildFourier[[i]]))
    if (!inherits(fit, "try-error")) {
        cat("i = ", i, ":\tconv = ", fit$conv, ",\tnegative LL = ", fit$value, "\n")
        modFourier[[i]] <- buildFourier[[i]](fit$par)
    }
}

i <- 1 # Best model, according to AIC/BIC
diag(W(modFourier[[i]]))[1:2]
drop(V(modFourier[[i]]))
fit <- dlmMLE(spots, c(-1, 0.1, 0.2, 6), buildFourier[[i]])
fit

### estimating standard errors of MLEs
fit <- dlmMLE(spots, c(-1, 0.1, 0.2, 6), buildFourier[[i]], hessian = TRUE)
fit
## for log variances
avar0 <- solve(fit$hessian)
avar0 # asymptotic variance matrix
sqrt(diag(avar0)) # asymptotic standard errors MLEs
## use Delta method for an AVar of variance estimates 
avar1 <- diag(c(exp(fit$par[1:3])), nr = 3, nc = 4) %*% avar0 %*%
    t(diag(c(exp(fit$par[1:3])), nr = 3, nc = 4))
avar1
sqrt(diag(avar1)) # asymptotic standard errors of log variance estimates
