
# Code to fit 'Closed Population N-mixture Model' in the 'unmarked' package

# Addressing the question: "How does wildfire influence salamander abundance?"

# Code last updated on 6/14/2021 by Gabe Barrile



# install the packages that we will use:

#install.packages("unmarked")
#install.packages("reshape")

# make sure packages are installed and active
require(unmarked)
require(reshape)


# citing the unmarked package
citation("unmarked")

# read-in the salamander count data from our ten transects  

# set working directory (which will be on the Git page)
setwd("H:/WEST_video_course/3_Closed_BinNMix_unmarked")

#  read-in data from the csv
df <- read.csv("Salamander_Wildfire.csv")

# so now we have our data stored as 'df'

# check out our data
df

# how many transects did we survey?
unique(df$Transect) # 10

# order the data frame by transect
df <- with(df,df[order(Transect),])
df

# 1.
# format count data for unmarked
m <- melt(df,id.var=c("Transect","Survey"),measure.var="Count")
head(m)

y=cast(m, Transect ~ Survey)

K <- ncol(y)

C <- as.matrix(y[,2:K])

C # each transect is a row (10 rows)
# each column indicates the survey at each site (three surveys at each transect)

head(C)
head(df)


# 2.
# format covariates for abundance and detection probability
# hypothesized that percent burned area influences abundance
# hypothesized that time of day influences detection probability


# Site covariates versus Observation covariates

# Site covariates do not change during each survey at a given site
# i.e., one value for a site covariate at each transect
# e.g., ten values for percent burned area, one value for each transect 

# We are interested in the abundance at a given transect, which we are assuming
# does not change over our three surveys at that transect. Therefore, if we
# want to estimate transect abundance as a function of percent burned area,
# then the percent burned area also cannot change over our three surveys at that
# given transect


# Observation covariates can change during each survey at a given site
# Observation covariate data should match survey data (i.e., matrix of counts)
# We can have a different detection probability for each survey
# e.g., different time of day for each survey.




# Format percent burned area as a covariate on abundance (as a 'site' covariate)
burned <- unique(df[,c("Transect","Burned")])
burned <- as.matrix(burned[,"Burned"])
burned

C
burned
head(df)



# Format time of day as a covariate on detection (as an 'observation' covariate)
m<-melt(df,id.var=c("Transect","Survey"),measure.var="Time")
head(m)
y=cast(m, Transect ~ Survey)

K <- ncol(y)

time <- as.matrix(y[,2:K])

C
time
head(df)

# remove objects that we don't need
rm(m,y,K)

# Here are the data for analysis:
C      # matrix of survey data (salamander counts)
burned # percent burned area at each transect (we think it might influence abundance)
time   # time recorded during every survey (we think it might influence detection)

# Input data into an 'unmarked data frame'
umf <- unmarkedFramePCount(
  y=C,                                   # Counts matrix
  siteCovs= data.frame(burned = burned), # Site covariates
  obsCovs = list(time = time))           # Observation covariates

# look at the summary our dataframe
summary(umf)



# Fit Closed Binomial N-mixture Model

# linear model for p (detection) follows first tilde, 
# then comes linear model for abundance (may see abundance ref. to as lambda for this model)
m1 <- pcount(~time ~burned, data=umf, K=130)
# Here 'K' is the upper summation limit for the summation over the random effects
# in the integrated likelihood. In unmarked, the default choice of K is the maximum
# observed count plus 100. That should normally be adequate, but you can play with K
# to determine the sensitivity of estimates to the value set for K.
max(C)

# Model-fitting function 'pcount' stands for 'point count' as this model can be
# employed with point count data. However, applications of the binomial mixture
# model are not restricted to point count data, as we can see with our use
# of counts along transects.

# check out summary of the model
summary(m1)


# We can choose alternate models for abundance other than the Poisson distribution

# Negative Binomial 
m2 <- pcount(~time ~burned, data=umf, mixture="NB", K=130)

# Zero-inflated Poisson
m3 <- pcount(~time ~burned, data=umf, mixture="ZIP", K=130)

# Compare models using AIC
cbind(AIC.P=m1@AIC, AIC.NB=m2@AIC, AIC.ZIP=m3@AIC)


# Predictions of abundance at specified values of percent burned area, say 0, 30, and 60)
newdat <- data.frame(burned=c(0, 30, 60))
predict(m1, type="state", newdata=newdat, append = T)
# 'predict' here uses delta rule to compute SEs and 95% CIs


# Predictions of p (detection probability) for values of time (e.g., 7am, 9am, and 11am)
newdat <- data.frame(time=c(7,9,11))
predict(m1, type="det", newdata=newdat, append = T)


# Visualize covariate relationships

# For abundance, predict to a new dataframe with 
# a suitable range for % burned area values
range(burned)
newdat <- data.frame(burned=seq(0, 60, length.out = 40))
pred.lam <- predict(m1, type="state", newdata=newdat, appendData = TRUE)

# For detection, predict to a new dataframe with 
# a suitable range for time of day values
range(time)
newdat <- data.frame(time=seq(7, 12, length.out = 20))
pred.det <- predict(m1, type="det", newdata=newdat, appendData = TRUE)


# plot abundance relationship with percent burned area
min(pred.lam$lower)
max(pred.lam$upper)
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(x = pred.lam$burned, y = pred.lam$Predicted, pch=16, 
     ylab = "Salamander Abundance",
     xlab = "% Burned Area", cex.lab=1.5, cex.axis=1.2, col="darkgray", ylim=c(0,65))
box(lwd = 4, col = 'black')
lines(pred.lam$burned, pred.lam$Predicted, lwd=8, col="blue")
lines(pred.lam$burned, pred.lam$lower, lwd=4, lty=2, col="black")
lines(pred.lam$burned, pred.lam$upper, lwd=4, lty=2, col="black")

# plot detection relationship with time of day
min(pred.det$lower)
max(pred.det$upper)
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(x = pred.det$time, y = pred.det$Predicted, pch=16, 
     ylab = "Detection Probability",
     xlab = "Time of Day (24 hr clock)", cex.lab=1.5, cex.axis=1.2, 
     col="darkgray", ylim=c(0.45,1))
box(lwd = 4, col = 'black')
lines(pred.det$time, pred.det$Predicted, lwd=8, col="blue")
lines(pred.det$time, pred.det$lower, lwd=4, lty=2, col="black")
lines(pred.det$time, pred.det$upper, lwd=4, lty=2, col="black")



# Extract abundance at each transect
ranef(m1)






# code to explore on your own

plot(ranef(m1), xlim=c(0,35))


# goodness-of-fit
require(AICcmodavg)

m1.gof <- Nmix.gof.test(m1, nsim = 100) # nsim is 100 only for illustrative purposes,
# you should probably run at least 1000 bootstrap replicates in your analysis

# p-value here suggests that we fail to reject the null (CANNOT conclude that
# the observed data are statistically different from the expected values)

# if p-value was less than alpha (e.g., 0.05), by contrast, then we would conclude that
# the observed data are statistically different from the expected values (lack of fit)

# residuals 
m1.resid <- residuals(m1)


# Fit detection-naive GLM to counts and plot comparison
summary(fm.glm <- glm(c(C) ~ rep(burned, 3), family=poisson)) # p-naive  model
matplot(burned, C, xlab="% Burned Area", ylab="Counts", frame = F, cex = 1.5, 
        pch = 1, col = "black", ylim = c(0,50))
curve(exp(coef(fm.glm)[1]+coef(fm.glm)[2]*x), 0, 60, type ="l", lwd=3, add=TRUE)
lines(burned, predict(m1, type="state")[,1], col = "blue", lwd = 3)
legend(25, 45, c("'Poisson GLM' with p", "Poisson GLM without p"), 
       col=c("blue", "black"), lty = 1, lwd=3, cex = 1.2)












