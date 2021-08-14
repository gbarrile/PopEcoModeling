# Code to fit 'Robust Design Multi-State' model in the 'RMark' package

# Addressing the question: 
# How do boreal toad dispersal rates vary amoong the three ponds in our study area?"

# Code last updated on 8/9/2021 by Gabe Barrile



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

# read-in the boreal toad capture-mark-recapture data for dispersal
#  read-in data from the csv
df <- read.csv("data/BorealToad_Dispersal.csv")

# so now we have our data stored as 'df'

# how many unique individuals did we tag?
length(unique(df$Tag)) # 1129 toads

# order dataframe by animal id
df <- with(df,df[order(Tag),])
head(df, 8)

table(df$Pond2)

df$Pond2 <- df$Pond
df$Pond2 <- factor(df$Pond2)
levels(df$Pond2)[levels(df$Pond2)=="A"] <- "1"
levels(df$Pond2)[levels(df$Pond2)=="B"] <- "2"
levels(df$Pond2)[levels(df$Pond2)=="C"] <- "3"

df$Pond3 <- df$Pond2
df$Pond3 <- as.numeric(as.character(df$Pond3))

m <- melt(df, id.var=c("Tag","Survey"), measure.var="Pond3")

y=cast(m, Tag ~ Survey)

y[is.na(y)]=0

k<-dim(y)[2]

#function to create capture history strings
pasty<-function(x) 
{
  k<-ncol(x)
  n<-nrow(x)
  out<-array(dim=n)
  for (i in 1:n)
  {
    y<-(x[i,])*1
    out[i]<-paste(y,sep="",collapse="")
  }
  return(out)
}

#capture history data frame 
capt.hist<-data.frame(ch=pasty(y[,2:k]))

y$ch <- pasty(y[,2:k])

y$ch <- gsub("1","A", y$ch)
y$ch <- gsub("2","B", y$ch)
y$ch <- gsub("3","C", y$ch)



head(y)
# 'ch' is character, not numeric

table(nchar(y$ch))

# add variables to 'y' dataframe

# Location: each pond
y$Pond <- df$Pond[match(y$Tag, df$Tag)]


# bt will be our dataframe that we input into RMark
bt <- data.frame(ch = y$ch, 
                 freq = 1, 
                 pond = y$Pond, 
                 tag = y$Tag)


# remove unneeded objects
rm(capt.hist,m,y,k,pasty)

head(df)
head(bt)

# make pond a factor variable
bt$pond <- as.factor(as.character(bt$pond))
table(bt$pond) # 1129 individuals

# are we missing any data?
table(is.na(bt)) # no missing data = good!


#convert to RMark format

#specify time intervals for Robust Design
intervals <- c(
  0,0,0,0,0,0,0,0,1, # 2016
  0,0,0,0,0,1,       # 2017
  0,0,0,0,0,1,       # 2018
  0,0,0,0,1,         # 2019
  0,0)               # 2020


# process data in RMark format
toad.proc=process.data(bt, 
                       model="HCRDMS",
                       time.intervals=intervals,
                       begin.time = 2016)

# create design data
toad.ddl <- make.design.data(toad.proc)

head(toad.proc$data)
names(toad.proc$data)
str(toad.proc$data)
names(toad.proc)

rm(intervals)

names(toad.ddl)
# S = survival probability
# Psi = transition probability (transitioning between states)
# p = capture probability
# c = recapture probability


# Multi-state model

# 'States' in our models refer to individual breeding ponds and
# transition probabilities represent the probability of moving from
# one breeding site to another between seasons.

# This modeling approach assumes that no site transitions occurred
# within a breeding season

# In model parameterization, transition probability is conditional
# on survival (e.g. transition probability Psi(rst) represents the probability
# an individual in state r at time t survives and moves to state s at time t + 1).

# look at transitions in raw capture histories
transition.pairs(bt$ch)

# within years
transition.pairs(substr(bt$ch,1,9)) # 2016
transition.pairs(substr(bt$ch,10,15)) # 2017
transition.pairs(substr(bt$ch,16,21)) # 2018
transition.pairs(substr(bt$ch,22,26)) # 2019
transition.pairs(substr(bt$ch,27,29)) # 2020

# between each year
transition.pairs(substr(bt$ch,1,15)) # 2016-2017
transition.pairs(substr(bt$ch,10,21)) # 2017-2018
transition.pairs(substr(bt$ch,16,26)) # 2018-2019
transition.pairs(substr(bt$ch,22,29)) # 2019-2020



# fix capture probabilities to zero where there are surveys we didn't actually do
toad.ddl$p$fix <- NA

toad.ddl$p$fix[toad.ddl$p$session == "2017" & toad.ddl$p$time == "6" & toad.ddl$p$stratum== "A"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2018" & toad.ddl$p$time == "6" & toad.ddl$p$stratum== "A"] <- 0

toad.ddl$p$fix[toad.ddl$p$session == "2016" & toad.ddl$p$time == "5" & toad.ddl$p$stratum== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2016" & toad.ddl$p$time == "6" & toad.ddl$p$stratum== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2016" & toad.ddl$p$time == "7" & toad.ddl$p$stratum== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2016" & toad.ddl$p$time == "8" & toad.ddl$p$stratum== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2016" & toad.ddl$p$time == "9" & toad.ddl$p$stratum== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2017" & toad.ddl$p$time == "5" & toad.ddl$p$stratum== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2017" & toad.ddl$p$time == "6" & toad.ddl$p$stratum== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2018" & toad.ddl$p$time == "5" & toad.ddl$p$stratum== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2018" & toad.ddl$p$time == "6" & toad.ddl$p$stratum== "B"] <- 0

# now do the same for c
toad.ddl$c

toad.ddl$c$fix <- NA

toad.ddl$c$fix[toad.ddl$c$session == "2017" & toad.ddl$c$time == "6" & toad.ddl$c$stratum== "A"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2018" & toad.ddl$c$time == "6" & toad.ddl$c$stratum== "A"] <- 0

toad.ddl$c$fix[toad.ddl$c$session == "2016" & toad.ddl$c$time == "5" & toad.ddl$c$stratum== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2016" & toad.ddl$c$time == "6" & toad.ddl$c$stratum== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2016" & toad.ddl$c$time == "7" & toad.ddl$c$stratum== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2016" & toad.ddl$c$time == "8" & toad.ddl$c$stratum== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2016" & toad.ddl$c$time == "9" & toad.ddl$c$stratum== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2017" & toad.ddl$c$time == "5" & toad.ddl$c$stratum== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2017" & toad.ddl$c$time == "6" & toad.ddl$c$stratum== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2018" & toad.ddl$c$time == "5" & toad.ddl$c$stratum== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2018" & toad.ddl$c$time == "6" & toad.ddl$c$stratum== "B"] <- 0


# okay, we are ready to specify and fit some models

#### Parameters for the analysis are:
names(toad.ddl)


# Specify model 

# Capture probability
# allow capture probability to vary by year and in each pond
# set capture probability to equal recapture probability
head(toad.ddl$p)
p.timexpond = list(formula=~session * stratum, share = TRUE)

# Survival 
# allow survival to vary over time and at each pond
# i.e., will estimate a different survival probability between each year and at each pond
head(toad.ddl$S)
S.timexpond = list(formula =~time * stratum)

# Transition (dispersal) 
# To determine if individuals disperse from breeding ponds, we fit model that 
# allows for state transitons among all three breeding ponds and permit
# variation in movement probabilities across years. In other words, obtain
# year-specific estimates for the probability an individual moves from one breeding pond
# to each of the other two ponds. The probability of remaining at the same site is
# therefore estimated by subtraction (1 - cumulative probability of movement to
# other two sites).
head(toad.ddl$Psi)
Psi.timexpond=list(formula=~-1+stratum:tostratum*time)

# fit model

m.timexpond <- mark(toad.proc,
                    toad.ddl, 
                    model.parameters=list(S = S.timexpond,
                                          Psi = Psi.timexpond,
                                          p = p.timexpond))



# Specify model 

# Capture probability
# allow capture probability to vary by year and in each pond
# set capture probability to equal recapture probability
head(toad.ddl$p)
p.timexpond = list(formula=~session * stratum, share = TRUE)

# Survival 
# allow survival to vary over time and at each pond
# i.e., will estimate a different survival probability between each year and at each pond
head(toad.ddl$S)
S.timexpond = list(formula =~time * stratum)

# Transition (dispersal) 
# To determine if individuals disperse from breeding ponds, fit model that 
# allows for state transitons among all three breeding ponds, but does not permit
# variation in movement probabilities across years.
head(toad.ddl$Psi)
Psi.ponds=list(formula=~-1+stratum:tostratum)

# fit model

m.pond <- mark(toad.proc,
                    toad.ddl, 
                    model.parameters=list(S = S.timexpond,
                                          Psi = Psi.ponds,
                                          p = p.timexpond))



# Collect models
pema <- collect.models(type = "HCRDMS") #grabs all models in workspace
names(pema) # 2 models
temp <- pema$model.table
( pema.ModSet <- temp[ ,-c(1:5)] )
pema.ModSet

for(i in 1:length(pema)){
  sing <- pema[[i]]$results$singular
  mod <- pema[[i]]$model.name
  print(mod)
  print(sing) # last will be model.table: NULL
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load function to check for singular parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sing.check <- function(mod.set){
  out <- data.frame()
  for(i in 1:(length(mod.set)-1)){ # -1 b/c otherwise will return value for the table itself
    sing <- mod.set[[i]]$results$singular
    sing <- if(is.null(sing)) "all converged" else (sing)
    mod <- mod.set[[i]]$model.name
    temp <- data.frame(model = mod, mod.num = i, singular = sing)
    temp$singular <- as.character(temp$singular)
    out <- rbind(temp, out)
  }
  return(out)
} 


# Check for singular parameters w/ custom function for RMark model list
sing.out <- sing.check(pema)
sing.out # take a look
sing.out[sing.out$singular != "all converged", ] # mods with params that didn't converge
sing.out$mod.num[sing.out$singular != "all converged"] # just the mod numbers, which you'll need

# Create a new mod set without mods with convergence issues
bad.mods <- sing.out$mod.num[sing.out$singular != "all converged"] # just the bad mod #'s
model.list2 <- remove.mark(pema, bad.mods) # drop the bad mods via the mod #'s

# Look at AIC table
model.list2$model.table

AICtablep<- as.data.frame(model.list2$model.table)
AICtablep



# check to make sure it converged
m.timexpond$results$singular
m.pond$results$singular

# check model for beta issues
m.pond$results$beta

# check out real estimates
m.pond$results$real

# check out derived estimates
m.pond$results$derived


Psilist=get.real(m.pond,"Psi",vcv=TRUE)
Psivalues=Psilist$estimates
Psivalues <- Psivalues[Psivalues$Cohort==0,]


a0 <- TransitionMatrix(Psivalues[Psivalues$Time==0,],vcv.real=Psilist$vcv.real)
a1 <- a0[[1]]
a2 <- a0[[3]]
a3 <- a0[[4]]
a1 <- round(a1,3)
a2 <- round(a2,3)
a3 <- round(a3,3)
# est
tt0 <- t(a1)
est0 <- c(tt0)
# lcl
tt00 <- t(a2)
est00 <- c(tt00)
# ucl
tt000 <- t(a3)
est000 <- c(tt000)

dispersal <- data.frame(estimate=est0, 
                        lcl=est00, 
                        ucl=est000, 
                        Trans=c("AA","AB","AC",
                                "BA","BB","BC",
                                "CA","CB","CC"))


dispersal$Fid.Disp <- c("Fidelity","Dispersal","Dispersal",
                        "Dispersal","Fidelity","Dispersal",
                        "Dispersal","Dispersal","Fidelity")


# plot the values
plot <- ggplot(dispersal, aes(x=Trans, y=estimate, group=Fid.Disp, color=Fid.Disp)) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=0.2, size=2, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), stat="identity", size = 7) +
  xlab(NULL) +
  ylab("Transition probability \n") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw()+
  theme(
    legend.position = "none",
    #panel.border = element_blank(),
    panel.background = element_rect(colour = "black", size=1, linetype = "solid"),
    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 20, vjust = 1, color = "black"), # spacing from y
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 12, color = "black")) +
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(size = 17))+
  theme(legend.position="top")+
  theme(strip.text = element_text(size = 7))

plot <- plot + geom_vline(xintercept = 3.5, linetype=2, size=0.25, col="grey8")

plot <- plot + geom_vline(xintercept = 6.5, linetype=2, size=0.25, col="grey8")

plot



# let's plot temporal variation in survival probability at each pond over our study period
surv <- m.pond$results$real[1:12,c(1,3,4)]

surv$year <- c("2016-2017","2017-2018","2018-2019","2019-2020",
               "2016-2017","2017-2018","2018-2019","2019-2020",
               "2016-2017","2017-2018","2018-2019","2019-2020")
surv$year <- as.factor(surv$year)

surv$Pond <- c("Pond A","Pond A","Pond A","Pond A",
               "Pond B","Pond B","Pond B","Pond B",
               "Pond C","Pond C","Pond C","Pond C")
surv$Pond <- as.factor(surv$Pond)
surv

# plot the values
ggplot(surv, aes(x=year, y=estimate, group=Pond, color=Pond)) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=0.2, size=1.25, position=position_dodge(0.3)) +
  geom_line(position=position_dodge(0.3), size=1.25) +
  geom_point(position=position_dodge(0.3), stat="identity", size = 6.5) +
  xlab(NULL) +
  ylab("Survival probability \n") +
  scale_y_continuous(limits = c(0.1, 1), breaks = seq(0, 1, 0.2), expand = c(0,0)) +
  theme_bw() +
  theme(#text = element_text(size = 18, family = "Times"), # controls all fonts in plot
    panel.background = element_rect(colour = "black", size=1, linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.2,"cm"),
    #axis.title.x = element_text(size = 22, color = "black"),
    axis.title.y = element_text(size = 22, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")) +
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(size = 17))+
  theme(legend.position="top")+
  theme(strip.text = element_text(size = 7))





