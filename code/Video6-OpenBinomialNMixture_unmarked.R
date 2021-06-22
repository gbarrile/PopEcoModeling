
# Code to fit 'Open Population N-mixture Model' in the 'unmarked' package

# Addressing the question: 
# How do annual rainfall patterns influence 
# temporal trends in Indian roller abundance?

# Code last updated on 6/21/2021 by Gabe Barrile



# install the packages that we will use:

#install.packages("unmarked")
#install.packages("reshape")
#install.packages("lattice")
#install.packages("ggplot2")

# make sure packages are installed and active
require(unmarked)
require(reshape)
require(lattice)
require(ggplot2)


# citing the unmarked package
citation("unmarked")


# set working directory (which will be on the Git page)
setwd("H:/WEST_video_course/6_Open_BinNMix_unmarked")

#  read-in data from the csv
df <- read.csv("IndianRoller_Rainfall.csv")

# how many site did we survey?
unique(df$Site) # 14

# order the data frame by transect
df <- with(df,df[order(Site),])
head(df)


# 1.
# format count data for unmarked

# histogram of counts
hist(df$Count, col = "grey", 
     main = "Counts of Indian Rollers", xlab = "", breaks=10)
abline(v = mean(df$Count), lwd = 3, col = "black")

# Plot counts over time for each site
xyplot(Count ~ Year | Site, data=df)


m <- melt(df,id.var=c("Site","Total.Survey"),measure.var="Count")
head(m)

y=cast(m, Site ~ Total.Survey)

K <- ncol(y)


C <- as.matrix(y[,2:K])

C # each site constitutes a row (14 rows)
# each column indicates the survey at each site (12 surveys at each site)

# Repeated count data with 4 primary periods (years) 
# and 3 secondary sampling periods (surveys within a year)
head(C,1)
head(df, 12)

tail(C,1)
tail(df, 12)




# Site-specific covariates (should equal number of sites)
# Format 'Agriculture' as a site covariate)
sc <- unique(df[,c("Site","Agriculture")])
sc <- data.frame(Agriculture = sc$Agriculture)
#agri <- as.matrix(agri[,"Agriculture"])
sc

hist(sc$Agriculture, col = "grey", 
     main = "Percent agricultural land at our sites", xlab = "", breaks=10)
abline(v = mean(sc$Agriculture), lwd = 3, col = "black")


# Plot values for agriculture over time for each site
xyplot(Agriculture ~ Year | Site, data=df)





# Format 'GrassCover' as an observation covariate
hist(df$GrassCover, col = "grey", 
     main = "Grass cover during our surveys", xlab = "", breaks=10)
abline(v = mean(df$GrassCover), lwd = 3, col = "black")


# Plot grass cover over time for each site
xyplot(GrassCover ~ Year | Site, data=df)


m<-melt(df,id.var=c("Site","Total.Survey"),measure.var="GrassCover")
head(m)
y=cast(m, Site ~ Total.Survey)

K <- ncol(y)

GrassCover <- y[,2:K]
GrassCover <- as.data.frame(GrassCover)
# GrassCover <- as.matrix(y[,2:K])

GrassCover


oc <- list(
  grasscover = as.matrix(GrassCover)
)

oc



# Format 'Rain' as a yearly-site covariate
hist(df$Rain, col = "grey", 
     main = "Annual rainfall (mm)", xlab = "")
abline(v = mean(df$Rain), lwd = 3, col = "black")


# Plot annual rainfall over time for each site
xyplot(Rain ~ Year | Site, data=df)


m<-melt(df,id.var=c("Site","Total.Survey"),measure.var="Rain")
head(m)
y=cast(m, Site ~ Total.Survey)

K <- ncol(y)

Rain <- y[,2:K]
Rain
Rain <- Rain[,c(1,4,7,10)]
Rain

Rain <- as.data.frame(Rain)
# GrassCover <- as.matrix(y[,2:K])

# Yearly-site covariates
ysc <- list(
  rain = as.matrix(Rain))

ysc


# Primary periods of surveys
primaryPeriod2 <- matrix(as.integer(c(
  1,2,3,4,
  1,2,3,4,
  1,2,3,4,
  1,2,3,4,
  1,2,3,4,
  1,2,3,4,
  1,2,3,4,
  1,2,3,4,
  1,2,3,4,
  1,2,3,4,
  1,2,3,4,
  1,2,3,4,
  1,2,3,4,
  1,2,3,4)), nrow=14, ncol=4, byrow=TRUE)

primaryPeriod2


# Here are the data for analysis:
C      # matrix of survey data (Indian roller counts)
sc     # site covariate (we think 'Agriculture' might influence initial abundance)
oc     # observation covariate (we think 'grasscover' might influence detection)
ysc    # yearly site covariate (we think 'rainfall' may influence trends/changes in abundance,
       # in our case lambda which is the finite rate of increase)
primaryPeriod2

# Create the unmarkedFrame
umf <- unmarkedFramePCO(
                         y=C,               # Counts matrix
                         siteCovs=sc,       # site covariates 
                         obsCovs=oc,        # observation covariates
                         yearlySiteCovs=ysc,# yearly site covariates
                         numPrimary=4,       
                         primaryPeriod=primaryPeriod2)

# Take a look
umf
summary(umf)


# about to fit model, but here's what the parameters mean:

# lambdaformula = initial abundance

# gammaformula = recruitment rate (when dynamics is "constant", "autoreg", or "notrend") 
#                OR population growth rate (when dynamics is "trend", "ricker", or "gompertz")

# omegaformula = apparent survival probability (when dynamics is "constant", "autoreg", or "notrend") 
#                OR equilibrium abundance (i.e. carrying capacity) (when dynamics is "ricker" or "gompertz")

# pformula = detection probability

# data = An object of class unmarkedFramePCO.

# mixture = character specifying mixture: "P", "NB", or "ZIP" 
            # for the Poisson, negative binomial, and zero-inflated Poisson distributions.

# K = Integer defining upper bound of discrete integration. 
      #This should be higher than the maximum observed count and high enough 
      #that it does not affect the parameter estimates. However, the higher 
      #the value the slower the compuatation.

# dynamics = see help page

# fix = If "omega", omega is fixed at 1. If "gamma", gamma is fixed at 0.

# starts = vector of starting values

# immigration = logical specifying whether or not to include an 
                # immigration term (iota) in population dynamics.

# iotaformula = Right-hand sided formula for average number of immigrants to a site per time step


# Fit model and backtransform
(m1 <- pcountOpen(
                  lambdaformula = ~Agriculture, 
                  gammaformula = ~ rain, 
                  omegaformula = ~ 1, 
                  pformula = ~ grasscover, 
                  data = umf, 
                  mixture = "P", 
                  K=200,
                  dynamics = "trend")) 

summary(m1)





# Predictions of growth rate at specified values of rainfall, say 300, 600, and 900)
newdat <- data.frame(rain=c(300, 600, 900))
predict(m1, type="gamma", newdata=newdat, append = T)
# Lambda = 1 (population stable)
# Lambda < 1 (population decreasing)
# Lambda > 1 (population increasing)


# Predictions of initial abundance at specified values of Agriculture, say 0, 50, and 100)
newdat <- data.frame(Agriculture=c(0, 50, 100))
predict(m1, type="lambda", newdata=newdat, append = T)

# Predictions of detection probability at specified values of grasscover, say 0, 20, and 40)
newdat <- data.frame(grasscover=c(0, 20, 40))
predict(m1, type="det", newdata=newdat, append = T)


# Visualize covariate relationships

# For lambda, predict to a new dataframe with 
# a suitable range for rainfall values
X <- seq(min(ysc$rain), max(ysc$rain), length.out=40)
newdat <- data.frame(rain = X)
pred <- round(predict(m1, type = "gamma", newdata = newdat, appendData = TRUE), 3)



# plot lambda relationship with annual rainfall patterns
min(pred$lower)
max(pred$upper)
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(x = pred$rain, y = pred$Predicted, pch=16, cex=0.1, 
     ylab = expression(paste("Pop Growth Rate","  ", "(", lambda, ")")),
     xlab = "Annual Rainfall (mm)", cex.lab=2, cex.axis=1.2, col="darkgray", ylim=c(0.6,1.4))
box(lwd = 3, col = 'black')
abline(h = 1, lwd=2, lty=3, col="darkgrey")
lines(pred$rain, pred$Predicted, lwd=6, col="blue")
lines(pred$rain, pred$lower, lwd=4, lty=2, col="black")
lines(pred$rain, pred$upper, lwd=4, lty=2, col="black")


# For detection, predict to a new dataframe with 
# a suitable range for grass cover values
X <- seq(min(oc$grasscover), max(oc$grasscover), length.out=40)
newdat <- data.frame(grasscover = X)
pred <- round(predict(m1, type = "det", newdata = newdat, appendData = TRUE), 3)


# plot detection relationship with grass cover
min(pred$lower)
max(pred$upper)
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(x = pred$grasscover, y = pred$Predicted, pch=16, cex=0.1, 
     ylab = "Detection Probability",
     xlab = "% Grass Cover", cex.lab=2, cex.axis=1.2, col="darkgray", ylim=c(0.4,0.8))
box(lwd = 3, col = 'black')
lines(pred$grasscover, pred$Predicted, lwd=6, col="blue")
lines(pred$grasscover, pred$lower, lwd=4, lty=2, col="black")
lines(pred$grasscover, pred$upper, lwd=4, lty=2, col="black")


# For initial abundance, predict to a new dataframe with 
# a suitable range for Agriculture values
X <- seq(min(sc$Agriculture), max(sc$Agriculture), length.out=40)
newdat <- data.frame(Agriculture = X)
pred <- round(predict(m1, type = "lambda", newdata = newdat, appendData = TRUE), 3)


# plot initial abundance relationship with agriculture
min(pred$lower)
max(pred$upper)
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(x = pred$Agriculture, y = pred$Predicted, pch=16, cex=0.1, 
     ylab = "Initial Abundance",
     xlab = "% Agricultural Land", cex.lab=2, cex.axis=1.2, col="darkgray", ylim=c(0,200))
box(lwd = 3, col = 'black')
lines(pred$Agriculture, pred$Predicted, lwd=6, col="blue")
lines(pred$Agriculture, pred$lower, lwd=4, lty=2, col="black")
lines(pred$Agriculture, pred$upper, lwd=4, lty=2, col="black")



# Estimates of conditional abundance distribution at each site
(re <- ranef(m1))
# Best Unbiased Predictors
bup(re, stat="mean")           # Posterior mean
bup(re, stat="mode")           # Posterior mode
confint(re, level=0.9) # 90% CI

# Plots
plot(re, subset=site %in% c(1), xlim=c(0,150))
plot(re, subset=site %in% c(7), xlim=c(0,150))

plot(re, subset=site %in% c(1,7), xlim=c(0,150))


# Can customize your own plots
# Let's plot mean abundance through time at a subset of our sites
pop <- as.data.frame(bup(re, stat="mean")) 
colnames(pop)
names(pop) <- c("2015","2016","2017","2018")
pop$Site <- 1:14
pop$Site <- as.factor(as.character(pop$Site))

pop$Agriculture <- df$Agriculture[match(pop$Site, df$Site)]

pop15 <- data.frame(Site = pop$Site, Year = "2015", Abundance = pop$`2015`)
pop16 <- data.frame(Site = pop$Site, Year = "2016", Abundance = pop$`2016`)
pop17 <- data.frame(Site = pop$Site, Year = "2017", Abundance = pop$`2017`)
pop18 <- data.frame(Site = pop$Site, Year = "2018", Abundance = pop$`2018`)

pop2 <- rbind(pop15, pop16, pop17, pop18)

pop2$Agriculture <- df$Agriculture[match(pop2$Site, df$Site)]
pop2


ag1 <- pop2[pop2$Site == "3" | pop2$Site == "11" | pop2$Site == "12", ]
ag1$Agriculture <- as.factor(as.character(ag1$Agriculture))

ag1$Agriculture=relevel(ag1$Agriculture,"100")
ag1$Agriculture=relevel(ag1$Agriculture,"50")
ag1$Agriculture=relevel(ag1$Agriculture,"0")

require(ggplot2)

  ggplot(ag1, aes(x=Year, y=Abundance, group=Agriculture, color=Agriculture)) +
  geom_point(size=5) + geom_line(size=1.5) +
  ylab("Site Abundance") +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40), expand = c(0,0)) +
  theme_bw() +
  theme(#text = element_text(size = 18, family = "Times"), # controls all fonts in plot
    panel.background = element_rect(colour = "black", size=1, linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.2,"cm"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")) +
  theme(legend.title = element_text(size=17))+
  theme(legend.text = element_text(size = 17))+
  theme(legend.position="top")+
  theme(strip.text = element_text(size = 7))





























