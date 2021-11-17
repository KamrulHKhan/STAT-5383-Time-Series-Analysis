######
###### Seemingly Unrelated Time Series Equations (SUTSE)
######

library(dlm)
library(astsa)

####
#### Example 1: Denmark and Spain investments
####
invest <- ts(matrix(scan("D:/Books/time_series/R_codes/denmark_spain.txt"), nc = 2, byrow = TRUE),
             start = 1960, names = c("Denmark", "Spain"))
plot(invest, type = "o", pch = 20, xlab = "", main = "Investments",
     mar = c(3, 4, 1, 2), oma = c(0, 0, 4, 0))

mod <- dlmModPoly(2)
mod$FF <- mod$FF %x% diag(2)
mod$GG <- mod$GG %x% diag(2)
mod$m0 <- rep(0, 2)
mod$C0 <- 1e7 * diag(4)

logCholeski <- function(u) {
    ## INPUT
    ## u: parameter vector of length n * (n + 1) / 2,
    ##    representing the log-Choleski parametrization of an n x n
    ##    variance matrix
    ## OUTPUT
    ## the variance matrix corresponding to u
    k <- length(u)
    n <- 0.5 * (sqrt(8 * k + 1) - 1)
    if (!all.equal(n, round(n))) error("Wrong length of the parameter vector")
    dd <- exp(tail(u, n))
    oo <- head(u, -n)
    V <- diag(dd)
    V[lower.tri(V)] <- oo
    tcrossprod(V)
}
    
build_invest <- function(psi) {
    psi <- matrix(psi, ncol = 2)
    ## W_mu <- logCholeski(psi[, 1])
    W_mu <- matrix(0, 2, 2)
    W_beta <- logCholeski(psi[, 1])
    mod$W <- bdiag(W_mu, W_beta)
    mod$V <- logCholeski(psi[, 2])
    mod
}

## try several random starting values, to be "sure" we don't find a local minimum
set.seed(1789)
n_start <- 10
out <- vector("list", n_start)
# start <- matrix(rnorm(9 * n_start, sd = 3), nrow = 9) 
start <- matrix(rexp(6 * n_start, rate = 1), nrow = 6)
start[c(1, 4), ] <- 0 # start from diagonal variance matrices

for (i in seq_len(n_start)) {
    fit <- try(dlmMLE(invest, start[, i], build_invest), silent = TRUE)
    if (!inherits(fit, "try-error")) {
        cat("i = ", i, ":\tconv = ", fit$conv, ",\tnegative LL = ", fit$value, "\n")
    } else 
        cat("i = ", i, ":\tError\n")
    out[[i]] <- fit
}

goodOut <- out[sapply(out, function(x) !inherits(x, "try-error") && x[["convergence"]] == 0)]
if (length(goodOut))
    jj <- which.min(sapply(goodOut, function(x) x[["value"]]))

fit_invest <- goodOut[[jj]]
mod_invest <- build_invest(fit_invest$par)
W(mod_invest)
V(mod_invest)

####
#### Example 2: Multivariate Capital Asset Pricing Model (CAPM)
####

###
### Univariate first...
###
capm <- read.table("D:/Books/time_series/R_codes/capm.txt", header = TRUE)
capm <- ts(capm, start = c(1978, 1), frequency = 12)
colnames(capm)
plot(capm, cex.lab = 0.7, main = "Monthly Returns", oma = c(3, 0, 5, 0))
IBM <- capm[, "IBM"] - capm[, "RKFREE"]
market <- capm[, "MARKET"] - capm[, "RKFREE"]

### static CAPM
outLM <- lm(IBM ~ market)
outLM$coef
summary(outLM)
## no obvious violations of linear model assumptions in the following plots
acf(outLM$res)
qqnorm(outLM$res); qqline(outLM$res)

### dynamic CAPM
buildCapm <- function(u) {
    dlmModReg(market, dV = exp(u[1]), dW = exp(u[2 : 3]))
}

outMLE <- dlmMLE(IBM, parm = rep(0, 3), buildCapm)
outMLE$conv
mod <- buildCapm(outMLE$par)
outS <- dlmSmooth(IBM, mod)
colnames(outS$s) <- c("alpha", "beta")
plot(dropFirst(outS$s), main = "Dynamic regression coefficients",
     mar = c(2, 4, 0, 1.1), cex.axis = 0.75)


###
### Multivariate Dynamic CAPM
###
y <- capm[, 1 : 4] - capm[, "RKFREE"]
colnames(y) <- colnames(capm)[1 : 4]
m <- NCOL(y)

### Set up the model
CAPM <- dlmModReg(market)
CAPM$FF <- CAPM$FF %x% diag(m)
CAPM$GG <- CAPM$GG %x% diag(m)
CAPM$JFF <- CAPM$JFF %x% diag(m)
CAPM$W <- CAPM$W %x% matrix(0, m, m)
CAPM$V <- CAPM$V %x% matrix(0, m, m)
CAPM$m0 <- rep(0, 2 * m)
CAPM$C0 <- diag(1e7, nr = 2 * m)

build_CAPM <- function(psi) {
    m <- NROW(CAPM$FF)
    k <- m * (m + 1) / 2
    CAPM$W[m + 1 : m, m + 1 : m] <- logCholeski(psi[1 : k])
    CAPM$V <- logCholeski(psi[k + 1 : k])
    CAPM
}

## fit model
set.seed(1789)
n_start <- 3
out <- vector("list", n_start)
start <- matrix(0, nrow = 20, ncol = n_start)
start[c(6 + 1:4, 16 + 1:4), ] <- -rexp(8 * n_start, rate = 1)

for (i in seq_len(n_start)) {
    fit <- try(dlmMLE(y, start[, i], build_CAPM, method = "SANN",
                      control = list(maxit = 5000)),
               silent = TRUE)
    if (!inherits(fit, "try-error")) {
        cat("i = ", i, ":\tconv = ", fit$conv, ",\tnegative LL = ", fit$value, "\n")
    } else 
        cat("i = ", i, ":\tError\n")
    out[[i]] <- fit
}

goodOut <- out[sapply(out, function(x) !inherits(x, "try-error") && x[["convergence"]] == 0)]
if (length(goodOut))
    jj <- which.min(sapply(goodOut, function(x) x[["value"]]))
## refine solution
fit <- dlmMLE(y, goodOut[[jj]]$par, build_CAPM, control = list(maxit = 500))
fit$convergence

CAPM <- build_CAPM(fit$par)
V(CAPM)
cov2cor(V(CAPM))
W(CAPM)
cov2cor(W(CAPM)[m + 1:m, m + 1:m])

### Smoothed betas
CAPMsmooth <- dlmSmooth(y, CAPM)

par(mar = c(3, 4, 1, 2) + 0.1, oma = c(0, 0, 0, 0), cex = 0.75)
plot(dropFirst(CAPMsmooth$s[, m + 1 : m]), lty = c("13", "6413", "431313", "B4"),
     plot.type = "s", xlab = "", ylab = "Beta")
abline(h = 1, col = "darkgrey")
legend("bottomright", legend = colnames(y), bty = "n",
       lty = c("13", "6413", "431313", "B4"), inset = 0.05)
