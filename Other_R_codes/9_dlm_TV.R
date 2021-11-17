######
###### Time-varying DLMs with package dlm
######

library(dlm)

### Nile river flow, accounting for the construction of Aswan dam, started in 1898
plot(Nile, type = 'o', pch = 20, xlab = "", main = "Nile river flow")
abline(v = 1898, lty = "dashed", col = "tomato")

####
#### Approach I: use a dummy intervention variable
####
x <- as.numeric(time(Nile) > 1898)
x
NileReg <- dlmModReg(x, dW = c(1, 0)) # constant dam effect
buildReg <- function(phi) {
    phi <- exp(phi)
    NileReg$V[] <- phi[1]
    NileReg$W[1, 1] <- phi[2]
    NileReg
}

## try random starting values, to make sure it is not a local minimum
set.seed(999)
n_start <- 10
out <- vector("list", n_start)
start <- matrix(rnorm(2 * n_start, sd = 4), nrow = 2) 

for (i in seq_len(n_start)) {
    fit <- try(dlmMLE(Nile, start[, i], buildReg), silent = TRUE)
    if (!inherits(fit, "try-error")) {
        cat("i = ", i, ":\tconv = ", fit$conv, ",\tnegative LL = ", fit$value, "\n")
    } else 
        cat("i = ", i, ":\tError\n")
    out[[i]] <- fit
}

goodOut <- out[sapply(out, function(x) !inherits(x, "try-error") && x[["convergence"]] == 0)]
if (length(goodOut))
    jj <- which.min(sapply(goodOut, function(x) x[["value"]]))

fitII <- goodOut[[jj]]
modII_fit <- buildII(fitII$par)
modII_fit


fitReg <- dlmMLE(Nile, parm = A[, 1], build = buildReg)
fitReg
NileReg <- buildReg(fitReg$par)

drop(V(NileReg))
diag(W(NileReg))

### filtering 
NileRegFilt <- dlmFilter(Nile, NileReg)
F <- cbind(1, x)
Nile_clean <- numeric(length(Nile))
Nile_clean_sd <- numeric(length(Nile))
for (j in seq_along(Nile)) {
    Nile_clean[j] <- crossprod(F[j, ], NileRegFilt$m[j+1, ])
    Nile_clean_sd[j] <- sqrt(F[j, , drop = FALSE] %*%
                             dlmSvd2var(NileRegFilt$U.C[[j+1]], NileRegFilt$D.C[j+1, ]) %*% F[j, ])
}
Nile_clean <- ts(Nile_clean, start = start(Nile))
lines(Nile_clean, col = "steelblue", lty = "longdash", lwd = 2)
lines(Nile_clean + qnorm(0.05, sd = Nile_clean_sd), col = "green", lty = "dashed", lwd = 2)
lines(Nile_clean + qnorm(0.95, sd = Nile_clean_sd), col = "green", lty = "dashed", lwd = 2)

####
#### Approach II: use a larger system variance to account for increased uncertainty
####

modII <- dlmModPoly(1)
JW(modII) <- matrix(1)
buildII <- function(psi) {
    psi <- exp(psi)
    V(modII) <- psi[1]
    X(modII) <- matrix(psi[2], nr = length(Nile), nc = 1)
    X(modII)[time(Nile) == 1899, ] <- psi[2] * psi[3]
    modII
}

## try random starting values, to make sure it is not a local minimum
set.seed(999)
n_start <- 10
out <- vector("list", n_start)
start <- matrix(rnorm(3 * n_start, sd = 4), nrow = 3) 

for (i in seq_len(n_start)) {
    fit <- try(dlmMLE(Nile, start[, i], buildII), silent = TRUE)
    if (!inherits(fit, "try-error")) {
        cat("i = ", i, ":\tconv = ", fit$conv, ",\tnegative LL = ", fit$value, "\n")
    } else 
        cat("i = ", i, ":\tError\n")
    out[[i]] <- fit
}

goodOut <- out[sapply(out, function(x) !inherits(x, "try-error") && x[["convergence"]] == 0)]
if (length(goodOut))
    jj <- which.min(sapply(goodOut, function(x) x[["value"]]))

fitII <- goodOut[[jj]]
modII_fit <- buildII(fitII$par)
modII_fit

sm <- dlmSmooth(Nile, modII_fit)
lines(sm$s)

### Standard errors
fit <- dlmMLE(Nile, start[, jj], buildII, hessian = TRUE)
H <- fit$hess
Delta <- solve(H)
sqrt(diag(Delta)) # estimated asymptotic standard errors
Dg <- diag(exp(fit$par))
Dg %*% Delta %*% t(Dg)
sqrt(diag(Dg %*% Delta %*% t(Dg))[1:2]) # estimated asymptotic standard errors of variance estimators


