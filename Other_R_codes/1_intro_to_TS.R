######
###### Examples of Time Series and introduction to statistical models
######

## pdf("intro_to_TS.pdf", width = 8, height = 8) # to output plots on a pdf file
## install.packages("astsa")
library(astsa)

## Note: tsplot is an astsa version 1.7.7+ script
## you can change tsplot to plot for an uglier graphic

## 1. Johnson & Johnson Earnings
tsplot(jj, type = "o", pch = 20, ylab = "$",
       main = "Johnson & Johnson - Quarterly Earnings per Share")

## 2. Global Temperature
tsplot(globtemp, type = "o", pch = 20, ylab = expression(degree * C),
       main = "Global Temperature Deviations")

## 3. Speech data
tsplot(speech)

## 4. Dow Jones Industrial Average
library(xts)                             
djiar <- diff(log(djia$Close))[-1]         # approximate returns
plot(djiar, main = "DJIA Returns", type = "l")  

## 5. SOI and Recruitment
par(mfrow = c(2, 1))  # set up the graphics
length(soi)
tsplot(soi, ylab = "", main = "Southern Oscillation Index")
tsplot(rec, ylab = "", main = "Recruitment") 

## 6. fMRI imaging
par(mfrow = c(2, 1), mar = c(3, 2, 1, 0) + .5, mgp = c(1.6, .6,0))  
ts.plot(fmri1[, 2:5], col = 1:4, ylab = "BOLD", xlab = "", main = "Cortex")
ts.plot(fmri1[, 6:9], col = 1:4, ylab = "BOLD", xlab = "", main = "Thalamus & Cerebellum")
mtext("Time (1 pt = 2 sec)", side = 1, line = 2)

## 7. Earthquakes and Explosions
par(mfrow = c(2, 1))
tsplot(EQ5, main = "Earthquake")
tsplot(EXP6, main = "Explosion")

## 8. White Noise and Moving Averages
w <- rnorm(500, 0, 1)  # 500 N(0,1) variates
v <- filter(w, sides = 2, rep(1/3, 3))  # moving average
par(mfrow = c(2, 1))
tsplot(w, ylim = c(-3, 3), main = "White Noise")
tsplot(v, ylim = c(-3, 3), main = "Moving Average")

## 10. Autoregression
par(mfrow = c(1, 1))
w <- rnorm(550, 0, 1)  # 50 extra to avoid startup problems
x <- filter(w, filter = c(1, -.9), method = "recursive")[-(1:50)]
tsplot(x, main = "Autoregression")

## 11. Random Walk with drift
set.seed(154) # so you can reproduce the results
w <- rnorm(200); x <- cumsum(w) # two commands in one line
wd = w + .2;    xd <- cumsum(wd)
tsplot(xd, ylim = c(-5, 55), main = "Random Walk", ylab = '')
lines(x, col = 4) 
abline(h = 0, col = 4, lty = 2)
abline(a = 0, b = .2, lty = 2)

## 12. Signal in noise
cs <- 2 * cos(2 * pi * (1:500) / 50 + .6 * pi)
w <- rnorm(500, 0, 1)
par(mfrow = c(3,1), mar = c(3, 2, 2, 1), cex.main = 1.5)   # help(par) for info
tsplot(cs, ylab = "", # help(plotmath) for more on math annotation in plots
       main = expression(x[t]==2*cos(2*pi*t/50+.6*pi)))
tsplot(cs + w, ylab= "", main = expression(x[t]==2*cos(2*pi*t/50+.6*pi)+N(0,1)))
tsplot(cs + 5 * w, ylab = "", main = expression(x[t]==2*cos(2*pi*t/50+.6*pi)+N(0,25)))

## 25. Sample ACF and scatterplots
(r <- round(acf(soi, 6, plot = FALSE)$acf[-1], 3)) # first 6 sample acf values
par(mfrow = c(1, 2))
plot(lag(soi, -1), soi); legend('topleft', legend = r[1])
plot(lag(soi, -6), soi); legend('topleft', legend = r[6])

