
# Code to fit 'Single-Season Occupancy Model' in the 'unmarked' package

# Addressing the question: 
# "How does forest cover influence the occurrence of brown tree snakes on pacific islands?
  
# Code last updated on 8/1/2021 by Gabe Barrile



# install the packages that we will use:

#install.packages("unmarked")
#install.packages("reshape")

# make sure packages are installed and active
require(unmarked)
require(reshape)


# citing the unmarked package
citation("unmarked")

# read-in the brown tree snake detection/nondetection data from the islands that we surveyed 
# read-in data from the csv
df <- read.csv("data/BrownTreeSnake_ForestCover.csv")

# how many islands did we survey?
unique(df$Island) # 12

# order the data frame by each island
df <- with(df,df[order(Island),])
head(df, 12)

# 1.
# format detection/nondetection data for unmarked
m <- melt(df,id.var=c("Island","Survey"),measure.var="BTS")
head(m)

y=cast(m, Island ~ Survey)

K <- ncol(y)

C <- as.matrix(y[,2:K])

C # each transect is a row (12 rows)
# each column indicates the survey at each island (six surveys at each island)

# Quick aside, and let's take island 2 as an example. Because we detected brown tree snakes
# at island 2 (during surveys 2, 3, 4, and 6), we assume that brown tree snakes were
# present on island 2 during all the other surveys -- it's just that we did not detect them
# (i.e., brown tree snakes were on island 2 during surveys 1 and 5, we just
# did not detect them during those surveys)

# Also note that detection/nondetection data at each island is sort of 
# analogous to our capture histories when we format our mark-recapture data.
# However, rather than individuals, here we deal with sites



# 2.
# format covariates for occupancy and detection probability
# hypothesized that forest cover influences occupancy
# hypothesized that temperature influences detection probability


# Site covariates versus Observation covariates

# Site covariates do not change during each survey at a given site
# i.e., one value for a site covariate at each island
# e.g., twelve values for forest cover, one value for each island 

# Observation covariates can change during each survey at a given site
# Observation covariate data should match survey data 
# (i.e., matrix of detection/nondetection data)
# We can have a different detection probability for each survey
# e.g., different temperature for each survey.


# Format forest cover as a covariate on occupancy (as a 'site' covariate)
forest <- unique(df[,c("Island","ForestCover")])
forest <- as.matrix(forest[,"ForestCover"])
forest


# Format temperature as a covariate on detection (as an 'observation' covariate)
m<-melt(df,id.var=c("Island","Survey"),measure.var="Temperature")

y=cast(m, Island ~ Survey)

K <- ncol(y)

temp <- as.matrix(y[,2:K])


# remove objects that we don't need
rm(m,y,K)


# Here are the data for analysis:
C      # matrix of survey data (detection/nondetection data of brown tree snakes)
forest # percent forest cover on each island (we think it might influence occupancy)
temp   # temp (°C) recorded during every survey (we think it might influence detection)


# Input data into an 'unmarked data frame'
umf <- unmarkedFrameOccu(
  y=C,                                   # detection/nondetection data
  siteCovs= data.frame(forest = forest), # Site covariates
  obsCovs = list(temp = temp))           # Observation covariates

# look at the summary our dataframe
summary(umf)


# Fit Single Season Occupancy Model

# linear model for p (detection) follows first tilde, 
# then comes linear model for occupancy 
m1 <- occu(~ 1 ~ 1, data=umf)

m2 <- occu(~ temp ~ 1, data=umf)

m3 <- occu(~ 1 ~ forest, data=umf)

m4 <- occu(~ temp ~ forest, data=umf)

# Compare models using AIC
cbind(m01=m1@AIC, m02=m2@AIC, m03=m3@AIC, m04=m4@AIC)
# model 4 has lowest AIC 

# check out summary of the model with the most support using AIC for model selection
summary(m4)



# Goodness of fit test on detection history frequencies
require(AICcmodavg)
system.time( gof.boot <- mb.gof.test(m4, nsim = 1000) )
gof.boot
# p-value here suggests that we fail to reject the null (CANNOT conclude that
# the observed data are statistically different from the expected values)

# if p-value was less than alpha (e.g., 0.05), by contrast, then we would conclude that
# the observed data are statistically different from the expected values (lack of fit)

# in other words, the observed frequency of the brown tree snake site-level
# detection histories agrees reasonably well with that expected under the 
# AIC-best model 'm4'. We therefore conclude that this model is suitable to
# use for inference and to inspect covariate relationships.



# Predictions of occupancy at specified values of percent forest cover, say 0, 50, and 100)
newdat <- data.frame(forest=c(0, 50, 100))
predict(m4, type="state", newdata=newdat, append = T)


# Predictions of p (detection probability) for values of temperature (e.g., 20, 30, 40 ?C)
newdat <- data.frame(temp=c(20,30,40))
predict(m4, type="det", newdata=newdat, append = T)


# Visualize covariate relationships

# For occupancy, predict to a new dataframe with 
# a suitable range for % forest cover values
range(forest)
newdat <- data.frame(forest=seq(0, 95, length.out = 40))
pred.occ <- predict(m4, type="state", newdata=newdat, appendData = TRUE)

# For detection, predict to a new dataframe with 
# a suitable range for temperature values
range(temp)
newdat <- data.frame(temp=seq(21, 40, length.out = 20))
pred.det <- predict(m4, type="det", newdata=newdat, appendData = TRUE)


# plot occupancy relationship with forest cover
min(pred.occ$lower)
max(pred.occ$upper)
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(x = pred.occ$forest, y = pred.occ$Predicted, pch=16, 
     ylab = "Probability of Occupancy",
     xlab = "Forest Cover (%)", cex.lab=1.5, cex.axis=1.2, col="darkgray", ylim=c(0,1))
box(lwd = 4, col = 'black')
lines(pred.occ$forest, pred.occ$Predicted, lwd=8, col="blue")
lines(pred.occ$forest, pred.occ$lower, lwd=4, lty=2, col="black")
lines(pred.occ$forest, pred.occ$upper, lwd=4, lty=2, col="black")

# plot detection relationship with time of day
min(pred.det$lower)
max(pred.det$upper)
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(x = pred.det$temp, y = pred.det$Predicted, pch=16, 
     ylab = "Detection Probability",
     xlab = "Air Temperature (?C)", cex.lab=1.5, cex.axis=1.2, 
     col="darkgray", ylim=c(0,1))
box(lwd = 4, col = 'black')
lines(pred.det$temp, pred.det$Predicted, lwd=8, col="blue")
lines(pred.det$temp, pred.det$lower, lwd=4, lty=2, col="black")
lines(pred.det$temp, pred.det$upper, lwd=4, lty=2, col="black")



# Empirical Bayes estimates of proportion of sites occupied
re <- ranef(m4)
sum(bup(re, stat="mode"))












