# Code to fit 'Cormack-Jolly-Seber' model in the 'RMark' package

# Addressing the question: 
# Does brook trout survival differ in a national park versus land developed for energy extraction?"

# Code last updated on 8/6/2021 by Gabe Barrile



# install the packages that we will use:

#install.packages("RMark")
#install.packages("reshape")
#install.packages("ggplot2")

# make sure packages are installed and active
require(RMark)
require(reshape)
require(ggplot2)


# citing the RMark package
citation("RMark")

# read-in the brook trout capture-mark-recapture data 
#  read-in data from the csv
df <- read.csv("data/BrookTrout_EnergyDevelopment.csv")

# so now we have our data stored as 'df'

# how many unique individuals did we tag?
length(unique(df$Tag)) # 146 fish

# order dataframe by animal id
df <- with(df,df[order(Tag),])
head(df, 8)

# format data for RMark
input.data <- df

input.data$detect <- rep(1,nrow(input.data))

m <- melt(input.data, id.var = c("Tag","Survey"), measure.var = "detect")

y = cast(m, Tag ~ Survey)

head(y)

y[is.na(y)] = 0
k <- dim(y)[2]

#if needing to deal with repeats on 1 occasion
y[,2:k] <- (y[,2:k]>0)*1
head(y)
head(df, 8)

#function to create capture history strings
pasty<-function(x) 
{
  k<-ncol(x)
  n<-nrow(x)
  out<-array(dim=n)
  for (i in 1:n)
  {
    y<-(x[i,]>0)*1
    out[i]<-paste(y,sep="",collapse="")
  }
  return(out)
}

# create capture history data frame 
capt.hist <- data.frame(ch = pasty(y[,2:k]))
head(capt.hist)

y$ch <- pasty(y[,2:k])

head(y)
# 'ch' is character, not numeric

nchar(y$ch)
table(nchar(y$ch))

# add variables to 'y' dataframe

# Location: national park or developed for energy extraction
y$Location <- df$Location[match(y$Tag, df$Tag)]


# bt will be our dataframe that we input into RMark
bt <- data.frame(ch = y$ch, 
                 freq = 1, 
                 location = y$Location, 
                 tag = y$Tag)


# remove unneeded objects
rm(capt.hist,input.data,m,y,k,pasty)

head(df)
head(bt)

# make location a factor variable
bt$location <- as.factor(as.character(bt$location))
table(bt$location)

# are we missing any data?
table(is.na(bt)) # no missing data = good!




# process data in RMark 

###### MODEL FOR SURVIVAL ################
d.proc=process.data(bt, 
                    model="CJS", 
                    groups = c("location"),
                    begin.time = 2015) 

# create design data
d.ddl <- make.design.data(d.proc)
# NOTE: the design dataframes are just as important as your raw data!

# Let's explore the design data
# see which parameters are estimated in the model
names(d.ddl)

d.ddl$Phi # apparent survival probability
# survival is termed apparent because mortality cannot be separated from emigration
# in other words, we don't really know whether an individual died or just left our study area
# one could make assumptions about a given system and perhaps it is reasonable to assume
# that individuals have a very low probability of moving out of the study area 
# and therefore apparent survival is a close approximation of true survival. 
# just depends on the dynamics of the system

d.ddl$p # capture probability





# Specify models (p and Phi)

# define p model
d.ddl$p # capture probability
# let's allow p to vary over time
p.loc = list(formula =  ~  location)


# define Phi model
d.ddl$Phi # apparent survival 
# let's allow survival to vary between locations
Phi.loc = list(formula =  ~ location)



# fit Model 1

m1 <- mark(d.proc,
           d.ddl, 
           model.parameters = list(p   = p.loc,
                                   Phi = Phi.loc))


# look at model output

# beta coefficients
m1$results$beta
# effects of each variable on survival (Phi) or capture (p), on the logit scale. 

# real estimates
m1$results$real
# estimates of survival probability and capture probability, 
# backtransformed to be on a regular linear scale from 0 to 1.






# Specify models (p and Phi)

# define p model
d.ddl$p # capture probability
# let's allow p to vary over time
p.time = list(formula =  ~  time)


# define Phi model
d.ddl$Phi # apparent survival 
# let's allow survival to vary between locations
Phi.loc = list(formula =  ~ location)



# fit Model 2

m2 <- mark(d.proc,
           d.ddl, 
           model.parameters = list(p   = p.time,
                                   Phi = Phi.loc))


# look at model output

# beta coefficients
m2$results$beta
# effects of each variable on survival (Phi) or capture (p), on the logit scale. 

# real estimates
m2$results$real
# estimates of survival probability and capture probability, 
# backtransformed to be on a regular linear scale from 0 to 1.


# Compare models using AICc
cbind(AIC.p.loc = m1$results$AICc, AIC.p.time = m2$results$AICc)
# model in which capture probability varied by location had lower AICc


# Let's plot predicitons based on model 1

# let's plot annual survival probability in the National Park versus
# site with energy development
surv <- m1$results$real[1:2,c(1,3,4)]

surv$location <- c("Energy Development","National Park")
surv$location <- as.factor(surv$location)
surv

# plot the values
ggplot(surv, aes(x=location, y=estimate)) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=0.2, size=2, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), stat="identity", size = 7) +
  xlab(NULL) +
  ylab("Apparent survival probability \n") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw()+
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 20, vjust = 1, color = "black"), # spacing from y
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"))


# let's plot capture probability in the National Park versus
# site with energy development
cap <- m1$results$real[3:4,c(1,3,4)]

cap$location <- c("Energy Development","National Park")
cap$location <- as.factor(cap$location)
cap

# plot the values
ggplot(cap, aes(x=location, y=estimate)) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=0.2, size=2, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), stat="identity", size = 7) +
  xlab(NULL) +
  ylab("Capture probability \n") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw()+
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 20, vjust = 1, color = "black"), # spacing from y
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"))




# Work with the results: estimate difference in survival between sites with energy development
# and sites within the national park

# Estimate the difference in survival rates between the two locations and 
# put a confidence interval on that estimated difference. 

# Here, we use the output from the model and the delta method
# to develop the SE and 95% CI for expressing the uncertainty for the estimated difference. 

# For information on the delta method, refer to Appendix B in the
# Program MARK manual (i.e., Program MARK: A Gentle Introduction)

# Note: some of the code below was developed by 
# Jay Rotella, Ecology Dept., Montana State University

# store var-cov for real estimates
rr=get.real(m1, "Phi", se=TRUE, vcv=TRUE)
sigma=rr$vcv.real # store var-cov in "sigma"

# load package 'msm' that implements the delta method
# then run delta method on difference in survival
library(msm)
s.Energy = rr$estimates$estimate[1]  # first row w/ output for sites with energy development
s.Park = rr$estimates$estimate[7]    # first row w/ output for sites within the national park
Diff= s.Energy - s.Park 
seDiff=deltamethod(~x1-x2,c(s.Energy,s.Park),sigma)
lclDiff=Diff-1.96*seDiff
uclDiff=Diff+1.96*seDiff
round(c(Diff,seDiff,lclDiff,uclDiff),3)
# 95% CI does not include zero, providing statistical evidence that 
# brook trout survival is lower at sites with energy development
# compared with sites located within the national park



