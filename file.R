#Final project: Final Report time series Project time series analysis Fall22
#Elham-Nasarian
#US-heart-failure-cases(1969-2020)

library(astsa)
library(TTR)

#Read the data into R
USdeath <- read.csv("C:\\Users\\elhamn20\\Documents\\fall2022\\Applied time sries\\final project\\FINALREPORT\\USdeath.csv")
WomenHeartDeath <- USdeath$WomenHeartDeath
MenHeartDeath <- USdeath$MenHeartDeath
Year <- USdeath$Year

#plot the data
tsplot(Year, WomenHeartDeath)
tsplot(Year, MenHeartDeath)

#acf plot
acf(WomenHeartDeath, lag.max = 40)
acf(MenHeartDeath, lag.max = 40)

#pacf plot
pacf(WomenHeartDeath, lag.max = 40)
pacf(MenHeartDeath, lag.max = 40)

#Compute square root for removing trends
WomenHeartDeath.sqrt <- sqrt(WomenHeartDeath)
tsplot(Year, WomenHeartDeath.sqrt,type="l")
MenHeartDeath.sqrt <- sqrt(MenHeartDeath)
tsplot(Year, MenHeartDeath.sqrt,type="l")

#Compute log death
WomenHeartDeath.log <- log(WomenHeartDeath)
tsplot(Year, WomenHeartDeath.log)
MenHeartDeath.log <- log(MenHeartDeath)
tsplot(Year, MenHeartDeath.log)

#Estimate the trend of the logarithm transformed time series. 
ma.WomenHeartDeath.log <- filter(WomenHeartDeath.log, sides=2, filter=rep(1/2,2))
detrended.WomenHeartDeath.log <- WomenHeartDeath.log - ma.WomenHeartDeath.log
plot(Year, ma.WomenHeartDeath.log,type="l", main="Trend (log(WomenHeartDeath))")

ma.MenHeartDeath.log <- filter(MenHeartDeath.log, sides=2, filter=rep(1/2,2))
detrended.MenHeartDeath.log <- MenHeartDeath.log - ma.MenHeartDeath.log
plot(Year, ma.MenHeartDeath.log,type="l", main="Trend (log(MenHeartDeath))")

# This series has a nonstationary trend.
# Let's compute the first-order difference to remove the trend.
WomenHeartDeath.diff.log <- diff(WomenHeartDeath.log)
tsplot( Year[1:51], WomenHeartDeath.diff.log)
MenHeartDeath.diff.log <- diff(MenHeartDeath.log)
tsplot( Year[1:51], MenHeartDeath.diff.log)

# Let us take a look at the sample ACF and sample PACF
acf(WomenHeartDeath.diff.log, lag.max = 40)
pacf(WomenHeartDeath.diff.log, lag.max = 40)
acf(MenHeartDeath.diff.log, lag.max = 40)
pacf(MenHeartDeath.diff.log, lag.max = 40)

# Compute AIC
ord.max = 1
WomenHeartDeath.diff.log.model <- ar(WomenHeartDeath.diff.log,order.max=ord.max,aic=TRUE,method="mle")
WomenHeartDeath.diff.log.model$aic
ord.max = 1
MenHeartDeath.diff.log.model <- ar(MenHeartDeath.diff.log,order.max=ord.max,aic=TRUE,method="mle")
MenHeartDeath.diff.log.model$aic

# Compute BIC
n = length(WomenHeartDeath.diff.log)
WomenHeartDeath.diff.log.BIC <- WomenHeartDeath.diff.log.model$aic - 2*((0:ord.max)+1) + log(n) * ((0:ord.max)+1)
WomenHeartDeath.diff.log.BIC       # Choose model with smallest BIC
# BIC = Bayesian Information Criterion
n = length(MenHeartDeath.diff.log)
MenHeartDeath.diff.log.BIC <- MenHeartDeath.diff.log.model$aic - 2*((0:ord.max)+1) + log(n) * ((0:ord.max)+1)
MenHeartDeath.diff.log.BIC       # Choose model with smallest BIC
# BIC = Bayesian Information Criterion

# Fit an AR(0) model with maximum likelihood
WomenHeartDeath.diff.log.model0 <- ar(WomenHeartDeath.diff.log,order.max=0,aic=FALSE,method="mle")
WomenHeartDeath.diff.log.model0
MenHeartDeath.diff.log.model0 <- ar(MenHeartDeath.diff.log,order.max=0,aic=FALSE,method="mle")
MenHeartDeath.diff.log.model0

order = 0
MLE.fit = ar.mle(MenHeartDeath.diff.log, order=order, aic=FALSE) 
MLE.fit$x.mean
MLE.fit$ar
sqrt(diag(MLE.fit$asy.var.coef))


# Now let's loop through several possible models to search for the best model
n = length(MenHeartDeath.log)
max.p = 1
max.d = 1
max.q = 0
max.P = 2
max.D = 1
max.Q = 2
BIC.array =array(NA,dim=c(max.p+1,max.d+1,max.q+1,max.P+1,max.D+1,max.Q+1))
AIC.array =array(NA,dim=c(max.p+1,max.d+1,max.q+1,max.P+1,max.D+1,max.Q+1))
best.bic <- 1e8
x.ts = MenHeartDeath.log
for (p in 0:max.p) for(d in 0:max.d) for(q in 0:max.q) 
  for (P in 0:max.P) for(D in 0:max.D) for(Q in 0:max.Q) 
  {
    # This is a modification of a function originally from the book:
    # Cowpertwait, P.S.P., Metcalfe, A.V. (2009), Introductory Time 
    # Series with R, Springer.
    # Modified by M.A.R. Ferreira (2016, 2020).
    cat("p=",p,", d=",d,", q=",q,", P=",P,", D=",D,", Q=",Q,"\n")
    fit <- tryCatch(
      {  arima(x.ts, order = c(p,d,q),  
               seas = list(order = c(P,D,Q), 
                           frequency(x.ts)),method="CSS-ML")
      },
      error = function(cond){
        message("Original error message:")
        message(cond)
        # Choose a return value in case of error
        return(NA)
      }
    )
    condition = !is.na(fit)
    if(length(condition)>1)condition = all(condition) 
    if(condition) { 
      number.parameters <- length(fit$coef) + 1
      BIC.array[p+1,d+1,q+1,P+1,D+1,Q+1] = -2*fit$loglik + log(n)*number.parameters
      AIC.array[p+1,d+1,q+1,P+1,D+1,Q+1] = -2*fit$loglik + 2*number.parameters
      if (BIC.array[p+1,d+1,q+1,P+1,D+1,Q+1] < best.bic) 
      {
        best.bic <- BIC.array[p+1,d+1,q+1,P+1,D+1,Q+1]
        best.fit <- fit
        best.model <- c(p,d,q,P,D,Q) 
      }
    }
  }
best.bic
best.fit
best.model


# Let's perform model diagnostics for the ARIMA
# with the sarima function:
sarima(WomenHeartDeath.log,0,0,0,1,1,1,12)
sarima(MenHeartDeath.log,0,0,0,1,1,1,12)


# Let's perform forecasting for 5 year using the ARIMA model:
WomenHeartDeath.log.for <- sarima.for(WomenHeartDeath.log,5,0,0,0,1,1,1,12)
MenHeartDeath.log.for <- sarima.for(MenHeartDeath.log,5,0,0,0,1,1,1,12)

# Here are predictions
exp(WomenHeartDeath.log.for$pred)
exp(MenHeartDeath.log.for$pred)

# Here are 95% prediction intervals
# Lower bounds
exp(WomenHeartDeath.log.for$pred - 1.96*WomenHeartDeath.log.for$se)
exp(MenHeartDeath.log.for$pred - 1.96*MenHeartDeath.log.for$se)

# Upper bounds
exp(WomenHeartDeath.log.for$pred + 1.96*WomenHeartDeath.log.for$se)
exp(MenHeartDeath.log.for$pred + 1.96*MenHeartDeath.log.for$se)


