
# Code to fit 'Single Season Occupancy Model' in the 'rjags' package

# Addressing the question: 
# "How does forest cover influence the occurrence of brown tree snakes on pacific islands?

# Code last updated on 7/5/2021 by Gabe Barrile



# install the packages that we will use:

#install.packages("rjags")
#install.packages("reshape")

# make sure packages are installed and active
require(rjags)
require(reshape)


# citing the rjags package
citation("rjags")

# read-in the brown tree snake detection/nondetection data from the islands that we surveyed 
# read-in data from the csv
#df <- read.csv("data/BrownTreeSnake_ForestCover.csv")

# set working directory (which will be on the Git page)
setwd("H:/WEST_video_course/7_Single-Season_Occupancy_unmarked")

#  read-in data from the csv
df <- read.csv("BrownTreeSnake_ForestCover.csv")

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




# 2.
# format covariates for occupancy and detection probability
# hypothesized that forest cover influences occupancy
# hypothesized that temperature influences detection probability


# Site covariates versus Observation covariates

# Site covariates do not change during each survey at a given site
# i.e., one value for a site covariate at each island
# e.g., nine values for forest cover, one value for each island 

# Observation covariates can change during each survey at a given site
# Observation covariate data should match survey data 
# (i.e., matrix of detection/nondetection data)
# We can have a different detection probability for each survey
# e.g., different temperature for each survey.

# Format forest cover as a covariate on occupancy (as a 'site' covariate)
forest <- unique(df[,c("Island","ForestCover")])
forest
forest <- forest$ForestCover

forest
head(df)




# Format temperature as a covariate on detection (as an 'observation' covariate)
m<-melt(df,id.var=c("Island","Survey"),measure.var="Temperature")

y=cast(m, Island ~ Survey)

K <- ncol(y)

temp <- as.matrix(y[,2:K])


# remove objects that we don't need
rm(m,y,K)



# In unmarked, our data included:
C
forest
temp
# For rjags, we need a few more bits of information:


# 3. We need a single integer value for the number of sites (to loop over)
# number of islands
Sites = dim(C)[1]
Sites

# 4. We need a single integer value for the number of surveys at each site (to loop over)
# number of surveys at EACH ISLAND, not total surveys
Surveys <- dim(C)[2]
Surveys


# Okay, now we have all the requisite information:

# Here are the data for our model/analysis:

Sites # number of sites
Surveys # number of surveys at each site (6 in this case)
C      # matrix of survey data (detection/nondetection data of brown tree snakes)
forest # percent forest cover on each island (we think it might influence occupancy)
temp   # temp (°C) recorded during every survey (we think it might influence detection)




### Bayesian Analysis using JAGS


# Define model
modelstring = " model {

# Biological model for true abundance
for (i in 1:Sites) {			# Loop over the number of islands
z[i] ~ dbern(psi[i])      # State model (Bernoulli distribution)
logit(psi[i]) <- alpha.occ + beta.occ * forest[i]


# The spatial variation of brown tree snake occurrence at island i (z[i]),
# for a collection of islands, is described by a Bernoulli distribution with
# the parameter psi, which denotes the probability of occupancy

# Occupancy probability at island i, on the logit scale, is a
# linear function of site-specific covariate (forest cover, 'forest') 
# with intercept (alpha.occ) and slope (beta.occ)

# we model the effect of the covariate on the logit scale 
# link function allows relationship between predictors and response
# to be modeled as a linear equation

# Observation model for replicated detection/nondetection data
for (j in 1:Surveys) {			# Loop over the number of surveys
C[i,j] ~ dbern(p.eff[i,j])  # Observation model (Bernoulli distribution)
p.eff[i,j] <- z[i] * p[i,j]
logit(p[i,j]) <- alpha.p + beta.p * temp[i,j]

 } # j
} # i


# The observed detections/nondetections (C[i,j]) at transect i and during replicate survey j
# are another Bernoulli random variable with a success rate that is the product 
# of the actual occurrence at site i, denoted by z[i], and the detection probability p, 
# (again at site i during survey j), which we termed p.eff. 

# The logit transform of detection probability at transect i during
# survey j is a logit-linear function of the site- and survey-
# specific covariate (temperature , 'temp') with intercept (alpha.p) and slope (beta.p)

# We model the effect of the covariate on the logit scale 


# Priors
alpha.occ ~ dunif(-10, 10) # uniform distribution for intercept in abundance model
beta.occ ~ dunif(-10, 10) # uniform distribution for slope in abundance model
alpha.p ~ dunif(-30, 30) # uniform distribution for intercept in detection model
beta.p ~ dunif(-10, 10) # uniform distribution for slope in detection model
# we chose vague priors, but could use more informative priors
# if you have better prior information 

# Derived quantities
N.occ <- sum(z[])	# total number of occupied islands (among those studied)

}"




# The above code may seem overwhelming, but that's mostly because of the comments
# Here is the model text without comments:

# Define model
# modelstring = " model {
# 
# for (i in 1:Sites) {		
# z[i] ~ dbern(psi[i])   
# logit(psi[i]) <- alpha.occ + beta.occ * forest[i]
# 
# for (j in 1:Surveys) {			
# C[i,j] ~ dbern(p.eff[i,j]) 
# p.eff[i,j] <- z[i] * p[i,j]
# logit(p[i,j]) <- alpha.p + beta.p * temp[i,j]
# 
# } 
# } 
# 
# alpha.occ ~ dunif(-10, 10)
# beta.occ ~ dunif(-10, 10) 
# alpha.p ~ dunif(-30, 30) 
# beta.p ~ dunif(-10, 10) 
# 
# N.occ <- sum(z[])
# 
# }"



# Inits function to set intial values for analysis
zst <- apply(C, 1, max) 

# We will use initial values from the analysis in 'unmarked'
inits <- function(){list(z = zst, alpha.occ = -5.6, beta.occ = 0.13, 
                          alpha.p = 21.2, beta.p = -0.7)}

# need initial value for each parameter that we will estimate

# adaptation
model <- jags.model(textConnection(modelstring),
                    data = list('Sites'=Sites, 'Surveys'=Surveys,
                                'C'=C, 'forest'=forest,
                                'temp'=temp),
                    inits = inits,
                    n.chains = 3, # number of parallel chains for the model 
                    n.adapt = 10000) # number of iterations for adaptation

# n.adapt or number of burn ins will vary depending on model convergence

# sample from posterior
samples <- coda.samples(model=model,
                        variable.names=c("z","N.occ","alpha.occ", "beta.occ",
                                         "alpha.p", "beta.p"), 
                        n.iter=60000, # number of iterations to monitor
                        thin=50) # thinning interval

# n.iter and thin will vary depending on model convergence


# check out summary of samples
summary(samples)
summary(samples)$statistics
summary(samples)$quantiles


# can plot samples and diagnostics (may need to increase size of your plot window)
#plot(samples)
#gelman.plot(samples)

# Estimate for the total number of islands occupied (out of the 12 that we surveyed)
summary(samples)$statistics['N.occ',]
# Compare with sum of maximum detections from field data
sum(apply(C, 1, max))


# Plot posterior distributions and means of intercept and slope parameters (occupancy model)
par(mfrow = c(2,1))
hist(unlist(samples[,'alpha.occ']), col = "grey", 
     main = "alpha.occ (intercept in occupancy model)", xlab="")
abline(v = summary(samples)$statistics['alpha.occ',][1], lwd = 3, col = "black")
hist(unlist(samples[,'beta.occ']), col = "grey", 
     main = "beta.occ (slope in occupancy model for forest covariate)", 
     xlab = "", xlim = c(-0.1,0.5), breaks = 50)
abline(v = summary(samples)$statistics['beta.occ',][1], lwd = 3, col = "black")
abline(v = 0, lwd = 3, col = "blue")

# Plot posterior distributions and means of intercept and slope parameters (Detection model)
par(mfrow = c(2,1))
hist(unlist(samples[,'alpha.p']), col = "grey", 
     main = "alpha.p (intercept in detection model)", xlab = "")
abline(v = summary(samples)$statistics['alpha.p',][1], lwd = 3, col = "black")
hist(unlist(samples[,'beta.p']), col = "grey", 
     main = "beta.p (slope in detection model for temperature covariate)", 
     xlab = "", xlim = c(-1.1,0))
abline(v = summary(samples)$statistics['beta.p',][1], lwd = 3, col = "black")
abline(v = 0, lwd = 3, col = "blue")


# Probability of occupancy on each island
summary(samples)$statistics[6:17,1:2]


# Plot occupancy probability as a function of forest cover
post.b0o = summary(samples)$statistics['alpha.occ',][1]
post.b1o = summary(samples)$statistics['beta.occ',][1]
range(forest)
newxo <- seq(0,95,length.out = 40)

par(mfrow = c(1,1))
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(newxo, 1/(1+exp(-(post.b0o + post.b1o*newxo))),
     ylab="Occupancy Probability", xlab="Forest Cover (%)",
     main="", las=1, type = "l", col="blue", lwd=4, ylim=c(0,1))
box(lwd = 3, col = 'black')


# Plot detection as a function of average air temperature
post.b0 = summary(samples)$statistics['alpha.p',][1]
post.b1 = summary(samples)$statistics['beta.p',][1]
range(temp)
newx <- seq(21,40,length.out = 20)

par(mfrow = c(1,1))
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(newx, 1/(1+exp(-(post.b0 + post.b1*newx))),
     ylab="Detection Probability", xlab="Air Temperature (°C)",
     main="", las=1, type = "l", col="blue", lwd=4, ylim=c(0,1))
box(lwd = 3, col = 'black')









