set.seed(7777)
library(astsa)
library(dlm)
library(DAAG)
#####################################################################################################

## to construct monthly data

data_1 = read.csv(file = paste0("Mortgage_Data.csv"), header = TRUE,row.names = NULL, sep=",")

data_temp = unname(data_1)[-1, -(1:2)]

data_full = matrix(as.numeric(as.matrix(data_temp)), nrow(data_temp),ncol(data_temp))

month = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
len = rep(31, 12)
len[2] = 28
len[c(4,6,9,11)] = 30
xx = data_1[-1, 2]
yy = rep(0,length(xx) )
for(m in 1:12)
	{
	for (i in 1:len[m])
		{
		aa = which(xx == paste0(i, "-", month[m]))
		yy[aa] = m
		}
	}

data_temp = unname(data_1)[-1, ]
year_pos = c()
for (i in 1990:2015) year_pos[i - 1989] = which(data_temp[ , 1] == i)
year_pos[length(year_pos)+1] = nrow(data_temp) + 1

ll = matrix(0, ncol = 12)
data_monthly = matrix(0, ncol = ncol(data_full))
for(i in 1:(length(year_pos)-1))
	{
	pos = year_pos[i]:(year_pos[i+1]-1)
	nn = c()
	for(j in 1:12) nn[j] = sum(yy[pos] == j)
	ll = rbind(ll, nn)
	temp_1 = data_full[pos, ]
	data_monthly = rbind(data_monthly, rowsum(temp_1, group = yy[pos])/nn)
	}
ll = ll[-1, ]
data_monthly = data_monthly[-1, ]

#tt = cbind(rep(1990:2015, each = 12), rep(1:12, 26), data_monthly)
#head(tt)

data = data_monthly[, 3:7]

Regions = c("NE", "SE", "NC", "SW", "W")

colnames(data) = Regions

interest_rate <- ts(data, start = c(1990,1), frequency = 12, names = c("NE", "SE", "NC", "SW", "W"))

interest_rate_list = list()
for (i in 1:5) interest_rate_list[[i]] = ts(data[ , i], frequency = 12, start = c(1990,1))

### data plot

pdf("data_visualization.pdf", height = 6, width = 12)
plot.ts(interest_rate, plot.type = "single", lty = c(2,4:6, 8), col = c(2:4,6,9), ylab = "Interest Rates")
legend("topright", bty="n", lty= c(2,4:6, 8), col=c(2:4,6,9),
       legend=Regions, title = "Regions")
dev.off()

#########################################################################
############################ SARIMA Model ###############################
#########################################################################

diagnosis_outputs_SARIMA = list()

IC_outputs_SARIMA = list()

Regions_full = c("North East", "South East", "North Center", "South West", "West")

pdf("SARIMA_model_forecasting_outputs.pdf", width = 6, height = 10)
par(mfrow = c(5,1))
for (region in 1:5)
	{
		log_interest_rate <- log(interest_rate_list[[region]]); 
		dldata <- diff(log_interest_rate)
		ddldata <- diff(dldata, 12)
		#par(mfrow = c(2,2))
		#acf(dldata, main = "ACF of Log Interest Rate Differences")
		#pacf(dldata, main = "PACF of Log Interest Rate Differences")
		#acf(ddldata, main = "ACF of Seasonal Difference of Log Interest Rate Differences")
		#pacf(ddldata, main = "ACF of Seasonal Difference of Log Interest Rate Difference")

		SARIMA_model <- arima(log_interest_rate, order = c(1, 1, 1), seasonal = c(0, 1, 1)) 

		diagnosis_outputs_SARIMA[[region]] = SARIMA_model
		
	### Information criterion based outputs

		AIC_SARIMA = - 2*SARIMA_model$loglik + 2*(length(SARIMA_model$coef) + 1)
		BIC_SARIMA = - 2*SARIMA_model$loglik + log(length(log_interest_rate))*(length(SARIMA_model$coef) + 1)

		IC_outputs_SARIMA[[region]] = c(SARIMA_model$loglik, AIC_SARIMA, BIC_SARIMA)
		
	## Forecasting

		pred_interest_rate_SARIMA <- predict(SARIMA_model, n.ahead = 60)
		lU <- pred_interest_rate_SARIMA$pred + pred_interest_rate_SARIMA$se
		lL <- pred_interest_rate_SARIMA$pred - pred_interest_rate_SARIMA$se

	## Back to original scale
		U <- exp(lU); L <- exp(lL)
		
		tsplot(window(interest_rate_list[[region]], start = 1990), type = 'o', pch = 16, cex = 0.75, ylab = "Interest Rate", xlab = "Time", ylim = c(0, max(interest_rate_list[[region]])), main = Regions_full[region], xlim = c(1989.9, 2020.08))
	
		lines(exp(pred_interest_rate_SARIMA$pred), type = "o", pch = 20, col = "blue")
		
		xx <- c(time(U), rev(time(U))); yy <- c(L, rev(U))
		polygon(xx, yy, border = 8, col = gray(0.6, alpha = 0.2))
		lines(U, col = "darkgreen"); lines(L, col = "darkgreen")
	}

dev.off()

pdf("SARIMA_model_diagnostic_outputs.pdf")

for (region in 1:5) tsdiag(diagnosis_outputs_SARIMA[[region]])

dev.off()

#########################################################################
################### Dynamic Linear Model (DLM) ##########################
#########################################################################

DLM_model <- function(psi) 
	{
		psi[1:4] <- exp(psi[1:4])
		dlmModPoly(2, dW = c(psi[1:2]), dV = psi[3]) + 
			dlmModTrig(12, 2, dV = 0, dW = 0) +
				dlmModARMA(ar = c(ARtransPars(psi[5])), ma = c(psi[6]),  sigma2 = psi[4], dV = 0) 
     }

diagnosis_outputs_DLM = list()

IC_outputs_DLM = list()

pdf("DLM_model_forecasting_outputs.pdf", width = 6, height = 10)

par(mfrow = c(5,1))

for (region in 1:5)
	{

	## try several random starting values, to be "sure" we don't find a local minimum		
		
		n_start <- 10
	
		start <- matrix(rnorm(6 * n_start, sd = 1), ncol = n_start) 

		out_DLM <- vector("list", n_start)

		for (i in seq_len(n_start)) 
			{
				fit <- try(dlmMLE(interest_rate_list[[region]], start[, i], DLM_model), silent = TRUE)
				if (!inherits(fit, "try-error")) {
				cat("i = ", i, ":\tconv = ", fit$conv, ",\tnegative LL = ", fit$value, "\n")
				} else 
				cat("i = ", i, ":\tError\n")
			out_DLM[[i]] <- fit
			}

		goodOut_DLM <- out_DLM[sapply(out_DLM, function(x) !inherits(x, "try-error") && x[["convergence"]] == 0)]
		if (length(goodOut_DLM)) jj <- which.min(sapply(goodOut_DLM, function(x) x[["value"]]))

		fit_DLM <- goodOut_DLM[[jj]]

		myFit_DLM <- DLM_model(fit_DLM$par)

		filt_DLM <- dlmFilter(interest_rate_list[[region]], myFit_DLM)

#		tsdiag(filt_DLM)

		diagnosis_outputs_DLM[[region]] = filt_DLM

	## AIC
	
		AIC_DLM = 2 * fit_DLM$value + 2 * length(fit_DLM$par)

	## BIC
		BIC_DLM = 2 * fit_DLM$value + log(length(interest_rate_list[[region]])) * length(fit_DLM$par)

		IC_outputs_DLM[[region]] = c(-fit_DLM$value, AIC_DLM, BIC_DLM)
		
	### Use the best model to forecast temperature next 5 years

		pred_interest_rate_DLM <- dlmForecast(filt_DLM, nAhead = 60)

		tsplot(window(interest_rate_list[[region]], start = 1990), type = 'o', pch = 20, cex = 0.75, ylab = "Interest Rate", xlab = "Time", ylim = c(0, max(interest_rate_list[[region]])), main = Regions_full[region], xlim = c(1989.9, 2020.08))
	
	## point forecasts

		lines(pred_interest_rate_DLM$f, col = "blue", lwd = 2, type = 'o', pch = 16, cex = 0.75)
		
	## add 90% probability bands
		
		lines(pred_interest_rate_DLM$f + sqrt(unlist(pred_interest_rate_DLM$Q)) * qnorm(0.05), col = "darkgreen")
		lines(pred_interest_rate_DLM$f + sqrt(unlist(pred_interest_rate_DLM$Q)) * qnorm(0.95), col = "darkgreen")
	}

dev.off()	

pdf("DLM_model_diagnostic_outputs.pdf")

for (region in 1:5) tsdiag(diagnosis_outputs_DLM[[region]])

dev.off()

#########################################################################
############################ VAR Model ##################################
#########################################################################

library(vars)

	VAR_model <- VAR(interest_rate, p = 1, season = 12, type = "both")
	
	VAR_log_likelihood <- logLik(VAR_model)

	### Information criterion based outputs

	## m*(m+1)/2 many variance-covariance parameters

	n_parameter = length(Bcoef(VAR_model)) + m*(m+1)/2
## AIC
	AIC_VAR = - 2 * c(VAR_log_likelihood) + 2 * n_parameter
## BIC
	BIC_VAR = - 2 * c(VAR_log_likelihood) + log(nrow(interest_rate)) * n_parameter

## Forecasting

	pred_interest_rate_VAR = predict(VAR_model, n.ahead = 60, ci = 0.90)

	pdf("VAR_model_forecasting_outputs.pdf", width = 6, height = 10)

	par(mfrow = c(5,1))

	for(region in 1:5)
		{
			pred_temp = (pred_interest_rate_VAR)[[1]][[region]]

#			tsplot(interest_rate_list[[region]], type = "l", xlab = "",  xlim = c(1989.9, 2020.08), ylim = c(0, max(data)), ylab = Regions[region])

			tsplot(window(interest_rate_list[[region]], start = 1990), type = 'o', pch = 20, cex = 0.75, ylab = "Interest Rate", xlab = "Time", ylim = c(0, max(interest_rate_list[[region]])), main = Regions_full[region], xlim = c(1989.9, 2020.08))

		## point forecasts

			pred_temp_1 = ts(pred_temp[ , 1], frequency = 12, start = c(2016,1))

			pred_temp_2 = ts(pred_temp[ , 2], frequency = 12, start = c(2016,1))

			pred_temp_3 = ts(pred_temp[ , 3], frequency = 12, start = c(2016,1))

			lines(pred_temp_1, col = "blue", lwd = 2, type = "o", pch = 16, cex = 0.75)

		## add 90% probability bands

			lines(pred_temp_2, col = "darkgreen", type = 'l')

			lines(pred_temp_3, col = "darkgreen", type = 'l')

		}

	dev.off()
	