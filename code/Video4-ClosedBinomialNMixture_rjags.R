
# Code to fit 'Closed Population N-mixture Model' in the 'rjags' package

# Addressing the question: "How does wildfire influence salamander abundance?"

# Code last updated on 6/16/2021 by Gabe Barrile



# install the packages that we will use:

#install.packages("rjags")
#install.packages("reshape")

# make sure packages are installed and active
require(rjags)
require(reshape)


# citing the rjags package
citation("rjags")

# read-in the salamander count data from our ten transects  
# read-in data from the csv
# df <- read.csv("data/Salamander_Wildfire.csv")

# read-in the salamander count data from our ten transects  

# set working directory (which will be on the Git page)
setwd("H:/WEST_video_course/3_Closed_BinNMix_unmarked")

#  read-in data from the csv
df <- read.csv("Salamander_Wildfire.csv")

# so now we have our data stored as 'df'

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
burned
burned <- burned$Burned

burned
head(df)


# Format time of day as a covariate on detection (as an 'observation' covariate)
m<-melt(df,id.var=c("Transect","Survey"),measure.var="Time")
head(m)

y=cast(m, Transect ~ Survey)

head(y)
head(df)

y

K <- ncol(y)

time <- as.matrix(y[,2:K])


C
time
head(df)



# remove objects that we don't need
rm(m,y,K)



# In unmarked, our data included:
C
burned
time
# For rjags, we need a few more bits of information:


# 3. We need a single integer value for the number of sites (to loop over)
# number of transects
Sites = dim(C)[1]
Sites

# 4. We need a single integer value for the number of surveys at each site (to loop over)
# number of surveys at EACH TRANSECT, not total surveys
Surveys <- dim(C)[2]
Surveys


# Okay, now we have all the requisite information:

# Here are the data for our model/analysis:

Sites # number of sites
Surveys # number of surveys at each site (3 in this case)
C # matrix of survey data (salamander counts)
burned # percent burned area at each transect (we think it might influence abundance)
time # time recorded during every survey (we think it might affect detection)




### Bayesian Analysis using JAGS


# Define model
modelstring = " model {

# Likelihood
# Biological model for true abundance
for (i in 1:Sites) {			# Loop over the number of transects
N[i] ~ dpois(lambda[i])   # State model (Poisson distribution)
log(lambda[i]) <- alpha.lam + beta1.lam * burned[i]

# The spatial variation of local salamander abundance at site i (N[i]),
# for a collection of sites, is described by a Poisson distribution with
# the mean denoted by lambda

# Mean abundance at site i, on the scale of the natural logarithm, is a
# linear function of site-specific covariate (burned) with
# intercept (alpha.lam) and slope (beta1.lam)


# Observation model for replicated counts
for (j in 1:Surveys) {			# Loop over the number of surveys
C[i,j] ~ dbin(p[i,j], N[i]) # Observation model (Binomial distribution)
lp[i,j] <- alpha.p + beta1.p * time[i,j]
p[i,j] <- exp(lp[i,j])/(1+exp(lp[i,j]))
# In place of the above two lines, you also could use logit(), as shown below:
#logit(p[i,j]) <- alpha.p + beta1.p * time[i,j]
}
}

# The observed counts (C[i,j]) at site i and during replicate survey j
# are described by a binomial distribution with sample size (N[i])
# and detection probability (p[i,j]).

# The logit transform of detection probability at site i during
# survey j is a logit-linear function of the site- and survey-
# specific covariate (time) with intercept (alpha.p) and slope (beta1.p)


# Priors
alpha.lam ~ dunif(-10, 10) # uniform distribution for intercept in abundance model
beta1.lam ~ dunif(-10, 10) # uniform distribution for slope in abundance model
alpha.p ~ dunif(-10, 10) # uniform distribution for intercept in detection model
beta1.p ~ dunif(-10, 10) # uniform distribution for slope in detection model

# Derived quantities
totalN <- sum(N[])	# Estimate total population size across all sites

}"




# The above code may seem overwhelming, but that's mostly because of the comments
# Here is the model text without comments:

# modelstring = " model {
# 
# for (i in 1:Sites) {			
# N[i] ~ dpois(lambda[i])   
# log(lambda[i]) <- alpha.lam + beta1.lam * burned[i]
# 
# for (j in 1:Surveys) {			
# C[i,j] ~ dbin(p[i,j], N[i]) 
# logit(p[i,j]) <- alpha.p + beta1.p * time[i,j]
# }
# }
# 
# alpha.lam ~ dunif(-10, 10) 
# beta1.lam ~ dunif(-10, 10) 
# alpha.p ~ dunif(-10, 10) 
# beta1.p ~ dunif(-10, 10) 
# 
# totalN <- sum(N[])
# 
# }"




# Inits function to set intial values for analysis
Nst <- apply(C, 1, max) + 1 # maximum count at each transect + 1
inits <- function(){list(N = Nst, alpha.lam=rnorm(1, 0, 1), beta1.lam=rnorm(1, 0, 1), 
                         alpha.p=rnorm(1, 0, 1), beta1.p=rnorm(1, 0, 1))}
# need initial value for each parameter that we will estimate

# adaptation
model <- jags.model(textConnection(modelstring),
                    data = list('Sites'=Sites, 'Surveys'=Surveys,
                                'C'=C, 'burned'=burned,
                                'time'=time),
                    inits = inits,
                    n.chains = 3, # number of parallel chains for the model 
                    n.adapt = 10000) # number of iterations for adaptation

# n.adapt or number of burn ins should be around 10,000 for real analysis, 
# but may take a while to run

# sample from posterior
samples <- coda.samples(model=model,
                        variable.names=c("N","totalN","alpha.lam", "beta1.lam",
                                         "alpha.p", "beta1.p"), 
                        n.iter=60000, # number of iterations to monitor
                        thin=50) # thinning interval

# n.iter should be around 60,000 for real analysis
# thin should be around 50 for real analysis
# this can make the analysis take a long time

# check out summary of samples
summary(samples)
summary(samples)$statistics
summary(samples)$quantiles


# can plot samples and diagnostics (may need to increase size of your plot window)
#plot(samples)
#gelman.plot(samples)

# Total population size across all ten transects
summary(samples)$statistics['totalN',]
# Compare with sum of maximum counts from field data
sum(apply(C, 1, max))




# Plot posterior distributions and means of intercept and slope parameters (Abundance model)
par(mfrow = c(2,1))
hist(unlist(samples[,'alpha.lam']), col = "grey", main = "alpha.lam", xlab = "")
abline(v = summary(samples)$statistics['alpha.lam',][1], lwd = 3, col = "black")
hist(unlist(samples[,'beta1.lam']), col = "grey", main = "beta1.lam", xlab = "")
abline(v = summary(samples)$statistics['beta1.lam',][1], lwd = 3, col = "black")

# Plot posterior distributions and means of intercept and slope parameters (Detection model)
par(mfrow = c(2,1))
hist(unlist(samples[,'alpha.p']), col = "grey", main = "alpha.p", xlab = "")
abline(v = summary(samples)$statistics['alpha.p',][1], lwd = 3, col = "black")
hist(unlist(samples[,'beta1.p']), col = "grey", main = "beta1.p", xlab = "")
abline(v = summary(samples)$statistics['beta1.p',][1], lwd = 3, col = "black")

# Plot posterior distribution and mean of total population size
par(mfrow = c(1,1))
hist(unlist(samples[,'totalN']), col = "grey", main = "Total N", xlab = "")
abline(v = summary(samples)$statistics['totalN',][1], lwd = 3, col = "black")


# Population size at each transect
summary(samples)$statistics[1:10,1:2]


# Plot abundance as a function of percent burned area
BinMix.pred <- exp(summary(samples)$statistics['alpha.lam',][1] + 
                     summary(samples)$statistics['beta1.lam',][1] * burned)


par(mfrow = c(1,1))
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(burned, BinMix.pred, ylab="Salamander Abundance", xlab="% Burned Area",
     main="", las=1, type = "l", col="blue", lwd=4, ylim=c(0,50))
box(lwd = 3, col = 'black')


# Plot detection as a function of time of day
post.b0 = summary(samples)$statistics['alpha.p',][1]
post.b1 = summary(samples)$statistics['beta1.p',][1]
range(time)
newx <- seq(7,12,length.out = 20)

par(mfrow = c(1,1))
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(newx, 1/(1+exp(-(post.b0 + post.b1*newx))),
     ylab="Detection Probability", xlab="Time of Day",
     main="", las=1, type = "l", col="blue", lwd=4, ylim=c(0.5,0.9))
box(lwd = 3, col = 'black')









