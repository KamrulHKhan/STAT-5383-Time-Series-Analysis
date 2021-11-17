######
###### Introduction to DLMs in R
######

library(dlm)
library(astsa)

### Nile river flow
plot(Nile, ylab = "Flow", type = 'o', pch = 16)

### Using base R functions
fitNile <- StructTS(Nile, "level")
fitNile
lines(fitted(fitNile), lty = "dashed", lwd = 2)
lines(tsSmooth(fitNile), lty = "dotted", lwd = 2)

modNile <- dlm(FF = 1, GG = 1, V = 15099, W = 1469, m0 = 1000, C0 = 1e6)
FF(modNile) 
V(modNile)

## Kalman Filter
filterNile <- dlmFilter(Nile, modNile)
str(filterNile, 1) # everything you are getting...

## a useful function to plot filtered means and probability bounds
lines.dlmFiltered <- function(filt_obj, prob = 0.90) {
    ## INPUTS
    ## filt_obj: an object of class 'dlmFiltered', returned by 'dlmFilter'
    ## prob:     the probability level of probability bands
    ##
    ## Adds filtered values and probability bands to the current plot
    
    y <- as.ts(filt_obj$y)
    F <- filt_obj$mod$FF

    ## recover standard deviations from SVD of covariance matrices
    sdFilt <- sapply(with(filt_obj, dlmSvd2var(U.C, D.C)[-1]), function(x) F %*% x %*% t(F))
    sdFilt <- ts(sqrt(sdFilt))
    tsp(sdFilt) <- tsp(y)
    ## compute expected noise-free observation
    expected <- ts(drop(as.matrix(filt_obj$m)[-1, ] %*% t(F)))
    tsp(expected) <- tsp(y)
    
    alpha <- 0.5 * (1 - prob) # tail probability

    ## compute probability limits and plot, together with filtering mean
    lns <- cbind(expected, qnorm(alpha, sd = sdFilt) %o% c(-1, 1) + as.vector(expected)) 
    lines(lns[, 1], col = "hotpink", lty = "longdash", lwd = 2)
    for (i in 1 : 2) lines(lns[, i + 1], col = "tomato", lty = "dotdash", lwd = 2)
    invisible()
}

## try it out
plot(Nile, ylab = "Flow", type = 'o', pch = 16)
lines(filterNile)

## Diagnostic plots
tsdiag(filterNile)
qqnorm(residuals(filterNile, sd = FALSE))
qqline(residuals(filterNile, sd = FALSE))

## CO2 data set
co2_annual <- aggregate(co2, FUN = mean)
StructTS(co2_annual, "level")

## Local level model - a terrible idea for these data
modCO2 <- dlm(FF = 1, GG = 1, V = 1e-4, W = 1.847, m0 = 1000, C0 = 1e6)
filterCO2 <- dlmFilter(co2_annual, modCO2)

tsdiag(filterCO2)
qqnorm(residuals(filterCO2, sd = FALSE))
qqline(residuals(filterCO2, sd = FALSE))

## Local trend model - definitely more appropriate
co2_mod <- dlmModPoly(dW = c(0.180, 0.008), dV = 7e-6)
co2Filter <- dlmFilter(co2_annual, co2_mod)

tsdiag(co2Filter)
qqnorm(residuals(co2Filter, sd = FALSE))
qqline(residuals(co2Filter, sd = FALSE))

## plot
plot(co2_annual, ylab = "", type = 'o', pch = 16, main = "CO2 at Mauna Loa",
     ylim = c(325, 340), xlim = c(1971, 1980))
lines(co2Filter)
points(co2Filter$f, col = "steelblue") # one-step-ahead forecasts

