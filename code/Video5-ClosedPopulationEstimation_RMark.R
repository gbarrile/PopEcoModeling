
# Code to fit 'Closed Population Estimation' in the 'RMark' package

# Addressing the question: "How does boreal toad abundance vary across ponds?"

# Code last updated on 7/22/2021 by Gabe Barrile



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

# read-in the boreal toad capture-mark-recapture data from our three ponds  
#  read-in data from the csv
df <- read.csv("data/BorealToad_CaptureRecapture.csv")

# so now we have our data stored as 'df'

# how many unique individuals did we tag?
length(unique(df$Tag)) # 137

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

# Pond
y$Pond <- df$Pond[match(y$Tag, df$Tag)]

# SVL
svl <- aggregate(SVL ~ Tag, data=df, FUN=mean)
head(svl)

# match by tag number to add SVL to y dataframe
y$SVL <- svl$SVL[match(y$Tag, svl$Tag)]


# boto will be our dataframe that we input into RMark
boto <- data.frame(ch = y$ch, freq = 1, pond = y$Pond, 
                   tag = y$Tag,
                   svl = y$SVL)


# remove unneeded objects
rm(capt.hist,input.data,m,svl,y,k,pasty)

head(df)
head(boto)

# make pond a factor variable
boto$pond <- as.factor(as.character(boto$pond))
table(boto$pond)

# are we missing any data?
table(is.na(boto)) # no missing data = good!


# let's obtain mean and standard deviation (as error bars) for svls
# create custom function for summary stats
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x) # lower limit of error bars
  ymax <- m+sd(x) # upper limit of error bars
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# plot svl data
ggplot(boto, aes(x=pond, y=svl, fill=pond)) + 
  geom_violin(width=1, position=position_dodge(1)) + 
  theme(axis.title.x = element_blank()) +
  geom_dotplot(binaxis='y', stackdir='center', 
          position=position_dodge(1), dotsize = 0.55, binwidth = 1) +
  stat_summary(fun.data=data_summary, size=1) +
  ggtitle("Violin plot of Boreal Toad SVLs at three ponds") + 
  labs(caption = "Filled dots with error bars = mean +/- standard deviation; open dots = raw data") + 
  xlab(NULL) +
  ylab("SVL (mm)") + 
  theme_bw() +
  theme( # control different aesthetics
    legend.position="none",
    panel.border = element_rect(colour = "black", fill=NA, size=1, linetype = "solid"), # place border around plot
    panel.grid.major = element_blank(), # no major gridlines
    panel.grid.minor = element_blank(), # no minor gridlines
    axis.ticks.length = unit(0.1,"cm"), # length of tick marks on axes
    plot.title = element_text(size = 18, face = "bold"), # plot title specs
    axis.title.y = element_text(size = 18, vjust = 2, color = "black"), # y-axis specs
    axis.text.x = element_text(size = 14, color = "black"), # x-axis specs
    plot.caption = element_text(hjust = 0, size=10), # bottom caption specs
    axis.text.y = element_text(size = 14, color = "black")) # y-axis labels






# process data in RMark 

###### CLOSED MODELS FOR ABUNDANCE ################
d.proc=process.data(boto, model="Closed", # NO B, NO D, NO I, NO E
                    groups = c("pond")) 

# create design data
d.ddl <- make.design.data(d.proc)
# NOTE: the design dataframes are just as important as your raw data!
# design data in this example are fairly simple, but will get more complex
# and detailed in future videos

# Let's explore the design data
# see which parameters are estimated in the model
names(d.ddl)

# capture and recapture probabilities can be the same or can be different 
# let's say you catch a rabbit in a trap for the first time 
# (that's your capture probability)
# Then let's say that rabbit subsequently avoids your traps
# In that case, it's very likely that your recapture probability of that rabbit
# is different than the probability of capturing that rabbit the first time
# This is just one example whereby capture and recapture probability can differ

# look at design data for each parameter

# capture probability
d.ddl$p

# recapture probability
d.ddl$c
# all individuals captured during survey 1 were captured for the first time,
# thus recapture probability starts at time = 2 or survey 2

# the number of individuals never captured
d.ddl$f0
# the notation 'f0' originates from the frequency (count) of animals observed 0 times). 
# The f0 parametrization is useful computationally because f0 is bounded on the interval [0,???], 
# thus forcing the logical constraint that N ??? the number of individuals marked (M). 
# In fact, MARK uses the f0 parametrization for ease of computation by using the log link function 
# to constrain f0 ??? 0, and presents the results in terms of N as a derived parameter
# in other words, the likelihood is rewritten in MARK in terms of the number of individuals never caught, f0,
# such that f0 = N ??? M 
# (M = the number of individuals marked)





# specify models for each parameter

# Model 1

# p (capture probability)
pc.=list(formula= ~ 1, share = TRUE)
# share = TRUE indicates p = c or capture probability = recapture probability
# ~ 1 or the intercept model is often referred to as the 'constant' model, 
# meaning that capture probability is constant over time and space 
# (e.g., doesn't vary across surveys or at different ponds)

# f0 (number of individuals never captured)
f0.=list(formula= ~ 1)
# this model indicates that the same number of individuals were
# never captured at all ponds
# (e.g., we failed to capture and mark 10 individuals at
# Pond 1, 10 individuals at Pond 2, and 10 at pond 3)
# it does NOT mean that the abundance is the same at all ponds,
# just that the number we failed to capture is the same


# fit Model 1
# this is not in code for github
setwd("H:/WEST_video_course/5_Closed_PopEst_RMark/models")

m1 <- mark(d.proc,
           d.ddl, 
           model.parameters = list(p = pc.,
                                   f0= f0.))


# look at model output

# beta coefficients
m1$results$beta
# these are just the intercepts in this case

# real estimates (on the scale we tend to think on)
m1$results$real

# derived parameters
m1$results$derived

# how is abundance estimated?
# f0
f0 <- m1$results$real[2,1]
f0
# how many toads did we mark at each pond?
table(boto$pond)
# Add f0 to the number marked at each pond
c(f0 + 44, f0 + 18, f0 + 75)
# compare to the estimated abundances from the derived parameters in the model
m1$results$derived$`N Population Size`[,1]
# the numbers are identical 
# In other words, the model knows how many toads we tagged at each pond:
table(boto$pond)
# then it adds the estimated number of individuals that we failed to tag
f0
# Note that we can estimate a different f0 for each pond (below)




# Model 2

# Do not share p and c
# Capture probability does not equal recapture probability
# obtain an estimate for capture probability and for recapture probability

# p
p.=list(formula= ~ 1)

# c
c.=list(formula= ~ 1)

# f0
f0.=list(formula= ~ 1)

# fit Model 2
m2 <- mark(d.proc,
           d.ddl, 
           model.parameters = list(p = p.,
                                   c = c.,
                                   f0= f0.))


# look at model output

# beta coefficients
m2$results$beta
# these are just the intercepts in this case

# real estimates
m2$results$real

# compare models with AICc
c(m1$results$AICc,m2$results$AICc)




# Model 2 again
# different way to specify, but results should be the same

# Do not share p and c
# Capture probability does not equal recapture probability
# obtain an estimate for capture probability and for recapture probability

# p
d.ddl$p
p.c.shared=list(formula=~ 1 + c, share = TRUE)

# f0
f0.=list(formula= ~ 1)

# fit Model 2
m2b <- mark(d.proc,
           d.ddl, 
           model.parameters = list(p = p.c.shared,
                                   f0= f0.))


# look at model output

# real estimates
m2b$results$real
m2$results$real

# compare models with AICc
c(m1$results$AICc,m2$results$AICc,m2b$results$AICc)










# Model 3

# fit model wherein capture probability varies over time / across surveys
# share p and c

# p
d.ddl$p
p.time =list(formula= ~ time, share = TRUE)

# Different f0 for each pond
# f0
d.ddl$f0
f0.pond =list(formula= ~ pond)

# fit Model 3
m3 <- mark(d.proc,
           d.ddl, 
           model.parameters = list(p = p.time,
                                   f0= f0.pond))


# look at model output

# beta coefficients
m3$results$beta

# real estimates 
m3$results$real

# derived parameters
m3$results$derived

# plot mean detection probability for each survey
plot(1:4,m3$results$real[1:4,1], type="b", ylim = c(0.2,0.4),
     ylab = "Detection probability", xlab = "Survey")





# Model 4

# Obtain different capture probability for each survey at each pond

# We conducted 4 surveys at each of three ponds, so 12 surveys total
# Therefore we should end up with 12 different capture probabilities,
# one for each survey at each site

# share p and c

# p
d.ddl$p
p.timepond =list(formula= ~ time * pond, share = TRUE)

# Different f0 for each pond
# f0
d.ddl$f0
f0.pond =list(formula= ~ pond)

# fit model
m4 <- mark(d.proc,
           d.ddl, 
           model.parameters = list(p = p.timepond,
                                   f0= f0.pond))


# look at model output

# beta coefficients
m4$results$beta

# real estimates 
m4$results$real

# derived parameters
m4$results$derived

# compare models with AICc
c(m1$results$AICc,m2$results$AICc,m3$results$AICc, m4$results$AICc)


# let's plot the estimated abundance at each pond,
# using estimates from Model 4 (m4)

abund <- m4$results$derived$`N Population Size`
abund$pond <- c("Pond 1", "Pond 2", "Pond 3")
abund$pond <- as.factor(abund$pond)
abund

# plot the values
ggplot(abund, aes(x=pond, y=estimate, color=pond)) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=0.2, size=2, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), stat="identity", size = 7) +
  xlab(NULL) +
  ylab("Estimated Abundance \n") +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40)) +
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



# let's plot capture probability during each survey at each pond,
# again using estimates from Model 4 (m4)

cap <- m4$results$real[1:12,c(1,3,4)]
cap$pond <- c("Pond 1", "Pond 1", "Pond 1", "Pond 1",
              "Pond 2", "Pond 2", "Pond 2", "Pond 2",
              "Pond 3", "Pond 3", "Pond 3", "Pond 3")
cap$pond <- as.factor(cap$pond)

cap$survey <- c(1,2,3,4,
                1,2,3,4,
                1,2,3,4)

cap$survey <- as.factor(as.character(cap$survey))

cap

# plot the values
ggplot(cap, aes(x=survey, y=estimate, group=pond, color=pond)) +
  geom_point(size=5) + geom_line(size=1.5) +
  xlab("Survey") +
  ylab("Capture probability \n") +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2), expand = c(0,0)) +
  theme_bw() +
  theme(#text = element_text(size = 18, family = "Times"), # controls all fonts in plot
    panel.background = element_rect(colour = "black", size=1, linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.2,"cm"),
    axis.title.x = element_text(size = 22, color = "black"),
    axis.title.y = element_text(size = 22, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")) +
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(size = 17))+
  theme(legend.position="top")+
  theme(strip.text = element_text(size = 7))







# Now let's incorporate individual covariates (e.g., SVL) into our models

# process data in RMark 

###### CLOSED MODELS FOR ABUNDANCE ################
d.proc=process.data(boto, model="Huggins") 

# Huggins models condition abundance out of the likelihood,
# thus permiting the modelling of capture probability as a function of individual covariates


# create design data
d.ddl <- make.design.data(d.proc)
# NOTE: the design dataframes are just as important as your raw data!
# design data in this example are fairly simple, but will get more complex
# and detailed in future videos

# Let's explore the design data
# see which parameters are estimated in the model
names(d.ddl)

# capture and recapture probabilities can be the same or can be different 
# let's say you catch a rabbit in a trap for the first time 
# (that's your capture probability)
# Then let's say that rabbit subsequently avoids your traps
# In that case, it's very likely that your recapture probability of that rabbit
# is different than the probability of capturing that rabbit the first time
# This is just one example whereby capture and recapture probability can differ

# look at design data for each parameter

# capture probability
d.ddl$p

# recapture probability
d.ddl$c
# all individuals captured during survey 1 were captured for the first time,
# thus recapture probability starts at time = 2 or survey 2

# the number of individuals never captured
d.ddl$f0
# In this model, the likelihood is conditioned on the number of animals detected 
# and f0 therefore drops out of the likelihood. 
# These models contain only p and c, with abundance N estimated as a derived parameter.
# As noted earlier, the primary advantage of the Huggins data type is that 
# individual covariates can be used to model p and c.


# specify models for each parameter

# Model 1

# p (capture probability)
pc.=list(formula= ~ 1, share = TRUE)
# share = TRUE indicates p = c or capture probability = recapture probability
# ~ 1 or the intercept model is often referred to as the 'constant' model, 
# meaning that capture probability is constant over time and space 
# (e.g., doesn't vary across surveys or at different ponds)


# fit Model 1
# this is not in code for github
setwd("H:/WEST_video_course/5_Closed_PopEst_RMark/models")

m1 <- mark(d.proc,
           d.ddl, 
           model.parameters = list(p = pc.))


# look at model output

# beta coefficients
m1$results$beta
# these are just the intercepts in this case

# real estimates (on the scale we tend to think on)
m1$results$real

# derived parameters
m1$results$derived

# how is abundance estimated?

# what is the probability of an individual not being captured at all during our study?
# (1 - p)(1 - p)(1 - p)(1 - p)
# what is p, our capture probability?
p <- m1$results$real$estimate[1]
p
# the probability of an individual not being captured at all during our study?
p.no <- (1 - p)*(1 - p)*(1 - p)*(1 - p)
# so what is the probability of being captured at least once?
p.once <- 1 - p.no
# then we can say that the number of individuals that we marked was equal to
# the probability of being captured at least once times the total population size (N)

# so we need to obtain the total number that we marked 
marked <- nrow(boto)

# here's the equation:      p.once * N = marked
# rearrage to solve for N:  N = marked / p.once

# So N should be:
marked / p.once
# compare to derived estimate from model 
m1$results$derived





# svl as individual covariate on capture probability

# p (capture probability)
pc.svl =list(formula= ~ svl, share = TRUE)


# fit Model 2
# this is not in code for github
setwd("H:/WEST_video_course/5_Closed_PopEst_RMark/models")

m2 <- mark(d.proc,
           d.ddl, 
           model.parameters = list(p = pc.svl))


# look at model output

# beta coefficients
m2$results$beta

# real estimates (on the scale we tend to think on)
m2$results$real

# derived parameters
m2$results$derived


# plot relationship between capture probability and SVL
range(boto$svl)
# make predictions based on model (m2)
pred.svl <- covariate.predictions(m2, data=data.frame(svl=56:88),indices=c(1))$estimates

min(pred.svl$lcl)
max(pred.svl$ucl)
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(x = pred.svl$covdata, y = pred.svl$estimate, pch=16, 
     ylab = "Capture Probability",
     xlab = "SVL (mm)", cex.lab=1.5, cex.axis=1.2, 
     col="darkgray", ylim=c(0.1,0.5))
box(lwd = 4, col = 'black')
lines(pred.svl$covdata, pred.svl$estimate, lwd=8, col="blue")
lines(pred.svl$covdata, pred.svl$lcl, lwd=4, lty=2, col="black")
lines(pred.svl$covdata, pred.svl$ucl, lwd=4, lty=2, col="black")







