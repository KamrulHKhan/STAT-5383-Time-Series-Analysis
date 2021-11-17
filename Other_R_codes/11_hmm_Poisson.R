
####
#### Using package HiddenMarkov
####

library(HiddenMarkov)

### Number of major earthquakes
y <- c(13, 14, 8, 10, 16, 26, 32, 27, 18, 32, 36, 24, 22, 23, 22, 18, 25, 21, 21, 14,
       8, 11, 14, 23, 18, 17, 19, 20, 22, 19, 13, 26, 13, 14, 22, 24, 21, 22, 26, 21,
       23, 24, 27, 41, 31, 27, 35, 26, 28, 36, 39, 21, 17, 22, 17, 19, 15, 34, 10, 15,
       22, 18, 15, 20, 15, 22, 19, 16, 30, 27, 29, 23, 20, 16, 21, 21, 25, 16, 18, 15,
       18, 14, 10, 15, 8, 15, 6, 11, 8, 7, 18, 16, 13, 12, 13, 20, 15, 16, 12, 18,
       15, 16, 13, 15, 16, 11, 11)
y <- ts(y, start = 1900)
plot(y, type = 'h', pch = 20, main = "Number of major earthquakes in the world")
plot(table(factor(y, levels = 0:50)), type = 'h', ylab = "Counts",
     main = "Count distribution - Poisson?"); abline(h = 0)
points(0:50, length(y) * dpois(0:50, lambda = mean(y)), pch = 19, col = "springgreen3")
mean(y)
var(y)
acf(y)
pacf(y)

### Fit a Hidden Markov Model
m <- 2 # two components
mod <- dthmm(y, Pi = matrix(1/m, m, m), delta = rep(1, m) / m,
             distn = "pois", pm = list(lambda = seq(10, 30, length = m)),
             nonstat = FALSE)

## first, approximate solution, using EM algorithm
mod1 <- BaumWelch(mod, control = bwcontrol(maxiter = 10))

## refine the solution, using direct loglikelihood maximization via nlm
allmap <- function(mod, p)
{
    ##    maps vector back to delta, Pi and lambda
    m <- sqrt(length(p))
    mod$Pi <- vector2Pi(head(p, -m))
    mod$pm$lambda <- exp(tail(p, m))
    mod$delta <- compdelta(mod$Pi) # implied stationary distribution
    return(mod)
}

p <- c(Pi2vector(mod1$Pi), log(mod1$pm$lambda)) # initial value for optimizer

out <- nlm(neglogLik, p, object = mod, pmap = allmap,
           print.level = 1, gradtol = 0.000001, iterlim = 100)

modFit <- allmap(mod, out$estimate)

## compare parameter estimates
summary(mod)
summary(mod1)
summary(modFit)

layout(matrix(c(1, 1, 2, 3), 2, 2, TRUE))
## look at most likely states
states <- Viterbi(modFit)
states <- ts(states, start = start(y))

plot(y, type = 'h', ylab = "Earthquakes", main = "Major earthquakes over time")
text(time(y), y + 0.5, states, col = 6 * states - 2)

## posterior state probabilities
probs <- with(modFit, Estep(x, Pi, delta, distn, pm))$u
probs <- ts(probs, start = start(y))
plot(probs[, 2], ylab = expression(hat(pi)[~2]*'(t|n)'),
     main = "Posterior probability of component 2")
abline(h = .5, lty = 2)

## mixture components
dpoismix <- function(x, mixture)
{
    delta <- mixture$delta # weights
    m <- length(delta) # number of components
    ## Calculate share of likelihood for all data for one component
    like.component <- function(x, component) 
        delta[component] * dpois(x, mixture$pm$lambda[component])
    ## Create array with likelihood shares from all components over all data
    likes <- sapply(seq_len(m), like.component, x = x)
    ## Add up contributions from components
    rowSums(likes)
}

plot(table(factor(y, levels = 0:45)) / length(y), type = 'h',
     xlab = "Earthquakes", ylab = "Relative frequency / Probability",
     main = "Components and mixture", ylim = c(0, 0.1)); abline(h = 0)
xx <- seq(0, 55)
yy <- dpoismix(xx, modFit[c("delta", "pm")])
lines(xx, yy, lwd = 2, col = "lawngreen")
lines(xx, dpois(xx, modFit$pm$lambda[1]), lwd = 2, col = "royalblue")
lines(xx, dpois(xx, modFit$pm$lambda[2]), lwd = 2, col = "tomato")
legend(35, 0.08, c("mixture", "component 1", "component 2"), lwd = 2,
       lty = 1, col = c("lawngreen", "royalblue", "tomato"), bty = 'n')

### model selection
m.max <- 5
out <- vector("list", m.max - 1)
AIC <- numeric(m.max - 1)
BIC <- numeric(m.max - 1)

for (m in 2:m.max) {
    mod <- dthmm(y, Pi = matrix(1/m, m, m), delta = rep(1, m) / m,
                 distn = "pois", pm = list(lambda = seq(10, 30, length = m)),
                 nonstat = FALSE)
    
    ## first, approximate solution, using EM algorithm
    mod1 <- BaumWelch(mod, control = bwcontrol(maxiter = 10))

    ## refine the solution, using direct loglikelihood maximization via nlm
    p <- c(Pi2vector(mod1$Pi), log(mod1$pm$lambda)) # initial value for optimizer
    
    opt <- nlm(neglogLik, p, object = mod, pmap = allmap,
               print.level = 0, gradtol = 0.000001, iterlim = 150)

    modFit <- allmap(mod, opt$estimate)
    out[[m - 1]] <- modFit
    AIC[m - 1] <- -2 * logLik(modFit) + 2 * length(opt$estimate)
    BIC[m - 1] <- -2 * logLik(modFit) + log(length(y)) * length(opt$estimate)
}

plot(2 : m.max, AIC, ylim = c(670, 830), type = 'b', lty = "dashed",
     xlim = c(1, m.max + 1), ylab = "", xlab = "Number of states")
lines(2 : m.max, BIC, type = 'b', lty = "dashed")
text(x = m.max + 0.5, y = c(tail(AIC, 1), tail(BIC, 1)), c("AIC", "BIC"))

## use m = 3 states for the final model
m <- 3
modFit <- out[[m - 1]]

## posterior state probabilities
a <- with(modFit, forward(x, Pi, delta, distn, pm)) # values are on log scale
b <- with(modFit, backward(x, Pi, distn, pm)) # values are on log scale
probs <- a + b
probs <- probs - apply(probs, 1, max) + 1 # to avoid overflows in exponentials
probs <- exp(probs)
probs <- probs / rowSums(probs) # rescale to probabilities
probs <- ts(probs, start = start(y))
par(mfrow = c(m, 1), oma = c(0, 0, 3, 0))
for (j in 1 : m) 
    plot(probs[, j], type = 'h', ylim = c(0, 1), xlab = "", ylab = "",
         main = paste("component", j, "( lambda =", round(modFit$pm$lambda[j], 2), ")"))
title(main = "Posterior state probabilities", outer = TRUE, cex.main = 2)

## most likely sequence of states
states <- ts(Viterbi(modFit), start = start(y))

par(mfrow = c(1, 1))
plot(y, ylab = "Earthquakes", xlab = "")
abline(h = modFit$pm$lambda, col = "lightgrey")
points(time(states), modFit$pm$lambda[states], pch = 20, col = "orangered")
