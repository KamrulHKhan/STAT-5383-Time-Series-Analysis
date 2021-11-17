####
#### HMM with Normal observations
####

library(HiddenMarkov)
library(astsa) # for the weekly return data

## Weekly returns on S&P 500 
y <- 100 * ts(as.vector(sp500w), start = 2003, freq = 52)
mu <- mean(y)
plot(y, type = 'o', pch = 19, cex = 0.25, ylab = "", xlab = "", main = "Weekly returns on S&P 500")

### Fit a HMM with Normal observations
### model selection
m.max <- 5
out <- vector("list", m.max - 1)
AIC <- numeric(m.max - 1)
BIC <- numeric(m.max - 1)

allmap <- function(mod, p){
    ##    maps vector back to delta, Pi and sd
    m <- sqrt(length(p))
    mod$Pi <- vector2Pi(head(p, -m))
    mod$pm$sd <- exp(tail(p, m))
    mod$pm$mean <- rep(mu, m)
    mod$delta <- compdelta(mod$Pi) # implied stationary distribution
    mod
}

for (m in 2 : m.max) {
    mod <- dthmm(y, Pi = matrix(1/m, m, m), delta = rep(1, m) / m,
                 distn = "norm", pm = list(mean = rep(mu, m), sd = seq(10, 1, len = m)),
                 nonstat = FALSE)
    
    p <- c(Pi2vector(mod$Pi), log(mod$pm$sd)) # initial value for optimizer
    opt <- nlm(neglogLik, p, object = mod, pmap = allmap,
               print.level = 0, gradtol = 0.000001, iterlim = 250)
    
    modFit <- allmap(mod, opt$estimate)
    out[[m - 1]] <- modFit
    AIC[m - 1] <- -2 * logLik(modFit) + 2 * length(opt$estimate)
    BIC[m - 1] <- -2 * logLik(modFit) + log(length(y)) * length(opt$estimate)
}

rn <- range(c(AIC, BIC))
plot(2 : m.max, AIC, ylim = rn + c(-1, 1) * 0.1 * diff(rn), type = 'b', lty = "dashed",
     xlim = c(1, m.max + 1), ylab = "", xlab = "Number of states")
lines(2 : m.max, BIC, type = 'b', lty = "dashed")
text(x = m.max + 0.5, y = c(tail(AIC, 1), tail(BIC, 1)), c("AIC", "BIC"))

## Consider m = 3
m <- 3
modFit <- out[[m - 1]]
states <- ts(Viterbi(modFit), start = start(y), freq = frequency(y))

## look at most likely states
states <- Viterbi(modFit)
states <- ts(states, start = start(y), freq = frequency(y))

plot(y, type = 'o', pch = 19, cex = 0.25, xlab = "", ylab = "", main = "Weekly returns on S&P 500")
text(time(y), y + sign(y) * 0.25, states, col = 6 * states - 2)

### Fitted marginal distribution
dnormalmix <- function(x, mixture) {
    delta <- mixture$delta
    m <- length(delta)
    ## Calculate share of likelihood for all data for one component
    like.component <- function(x, component) 
        delta[component] * dnorm(x, mean = mixture$mean[component], sd = mixture$sd[component])
    ## Create array with likelihood shares from all components over all data
    likes <- sapply(seq_len(m), like.component, x = x)
    ## Add up contributions from components
    rowSums(likes)
}

hist(y, "FD", prob = TRUE, , col = "grey80",
     main = "Weekly returns on S&P 500", xlab = "Time"); box(); abline(h = 0)
xx <- seq(-22, 12, length = 201)
yy <- dnormalmix(xx, list(delta = modFit$delta, mean = modFit$pm$mean, sd = modFit$pm$sd))
lines(xx, yy, lwd = 2, col = "violet")


