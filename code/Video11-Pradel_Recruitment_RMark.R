# Code to fit 'Robust Design Pradel Recruitment' model in the 'RMark' package

# Addressing the question: 
# How does boreal toad recruitment and population growth vary over time at three ponds?"

# Code last updated on 7/16/2021 by Gabe Barrile



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

# read-in the boreal toad capture-mark-recapture data for recruitment
#  read-in data from the csv
df <- read.csv("data/BorealToad_Recruitment.csv")

# so now we have our data stored as 'df'

# how many unique individuals did we tag?
length(unique(df$Tag)) # 1199 toads

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

table(nchar(y$ch))

# add variables to 'y' dataframe

# Location: national park or developed for energy extraction
y$Pond <- df$Pond[match(y$Tag, df$Tag)]


# bt will be our dataframe that we input into RMark
bt <- data.frame(ch = y$ch, 
                 freq = 1, 
                 pond = y$Pond, 
                 tag = y$Tag)


# remove unneeded objects
rm(capt.hist,input.data,m,y,k,pasty)

head(df)
head(bt)

# make pond a factor variable
bt$pond <- as.factor(as.character(bt$pond))
table(bt$pond)

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
                       model="RDPdfHuggins",
                       time.intervals=intervals, 
                       groups = c("pond"), 
                       begin.time = 2016)

# create design data
toad.ddl <- make.design.data(toad.proc)

head(toad.proc$data)
names(toad.proc$data)
str(toad.proc$data)
names(toad.proc)

rm(intervals)

names(toad.ddl)
# Phi = apparent survival
# f = recruitment rate
# p = capture probability
# c = recapture probability


# The recruitment rate (f) represents the number of individuals added 
# to the breeding population at time t + 1 per animal in the breeding
# population at time t. 


# fix capture probabilities to zero where there are surveys we didn't actually do
toad.ddl$p$fix <- NA

toad.ddl$p$fix[toad.ddl$p$session == "2017" & toad.ddl$p$time == "6" & toad.ddl$p$pond== "A"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2018" & toad.ddl$p$time == "6" & toad.ddl$p$pond== "A"] <- 0

toad.ddl$p$fix[toad.ddl$p$session == "2016" & toad.ddl$p$time == "5" & toad.ddl$p$pond== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2016" & toad.ddl$p$time == "6" & toad.ddl$p$pond== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2016" & toad.ddl$p$time == "7" & toad.ddl$p$pond== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2016" & toad.ddl$p$time == "8" & toad.ddl$p$pond== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2016" & toad.ddl$p$time == "9" & toad.ddl$p$pond== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2017" & toad.ddl$p$time == "5" & toad.ddl$p$pond== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2017" & toad.ddl$p$time == "6" & toad.ddl$p$pond== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2018" & toad.ddl$p$time == "5" & toad.ddl$p$pond== "B"] <- 0
toad.ddl$p$fix[toad.ddl$p$session == "2018" & toad.ddl$p$time == "6" & toad.ddl$p$pond== "B"] <- 0

# now do the same for c
toad.ddl$c

toad.ddl$c$fix <- NA

toad.ddl$c$fix[toad.ddl$c$session == "2017" & toad.ddl$c$time == "6" & toad.ddl$c$pond== "A"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2018" & toad.ddl$c$time == "6" & toad.ddl$c$pond== "A"] <- 0

toad.ddl$c$fix[toad.ddl$c$session == "2016" & toad.ddl$c$time == "5" & toad.ddl$c$pond== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2016" & toad.ddl$c$time == "6" & toad.ddl$c$pond== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2016" & toad.ddl$c$time == "7" & toad.ddl$c$pond== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2016" & toad.ddl$c$time == "8" & toad.ddl$c$pond== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2016" & toad.ddl$c$time == "9" & toad.ddl$c$pond== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2017" & toad.ddl$c$time == "5" & toad.ddl$c$pond== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2017" & toad.ddl$c$time == "6" & toad.ddl$c$pond== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2018" & toad.ddl$c$time == "5" & toad.ddl$c$pond== "B"] <- 0
toad.ddl$c$fix[toad.ddl$c$session == "2018" & toad.ddl$c$time == "6" & toad.ddl$c$pond== "B"] <- 0


# okay, we are ready to specify and fit some models

#### Parameters for the analysis are:
names(toad.ddl)

head(toad.ddl$Phi)
head(toad.ddl$f)
head(toad.ddl$p)




# Specify model 

# Capture probability
# allow capture probability to vary by year and in each pond
# set capture probability to equal recapture probability
p.timexpond = list(formula=~session * pond, share = TRUE)

# Recruitment 
# allow recruitment to vary over time and at each pond
# i.e., will estimate a different recruitment rate between each year and at each pond
f.timexpond = list(formula=~time * pond)

# Survival 
# allow survival to vary over time and at each pond
# i.e., will estimate a different survival probability between each year and at each pond
Phi.timexpond = list(formula =~time * pond)

# fit model

m.timexpond <- mark(toad.proc,
                   toad.ddl, 
                   model.parameters=list(Phi = Phi.timexpond,
                                         f = f.timexpond,
                                         p = p.timexpond))



# Specify different model 

# Capture probability
# allow capture probability to vary by year and in each pond
# set capture probability to equal recapture probability
p.timexpond = list(formula=~session * pond, share = TRUE)

# Recruitment 
# allow recruitment to vary over time 
# i.e., will estimate a different recruitment rate between each year 
f.time = list(formula=~time)

# Survival 
# allow survival to vary over time
# i.e., will estimate a different survival probability between each year
Phi.time = list(formula =~time)

# fit model

m.time <- mark(toad.proc,
                   toad.ddl, 
                   model.parameters=list(Phi = Phi.time,
                                         f = f.time,
                                         p = p.timexpond))



# Specify another model 

# Capture probability
# allow capture probability to vary by year and in each pond
# set capture probability to equal recapture probability
p.timexpond = list(formula=~session * pond, share = TRUE)

# Recruitment 
# allow recruitment to vary at each pond
# i.e., will estimate a different recruitment rate at each pond
f.pond = list(formula=~pond)

# Survival 
# allow survival to vary at each pond
# i.e., will estimate a different survival probability at each pond
Phi.pond = list(formula =~pond)

# fit models

m.pond <- mark(toad.proc,
                   toad.ddl, 
                   model.parameters=list(Phi = Phi.pond,
                                         f = f.pond,
                                         p = p.timexpond))



# Specify model 

# Capture probability
# allow capture probability to vary by year and in each pond
# set capture probability to equal recapture probability
p.timexpond = list(formula=~session * pond, share = TRUE)

# Recruitment 
# allow recruitment to vary over time and at each pond
# i.e., will estimate a different recruitment rate between each year and at each pond
f.time.pond = list(formula=~time + pond)

# Survival 
# allow survival to vary over time and at each pond
# i.e., will estimate a different survival probability between each year and at each pond
Phi.time.pond = list(formula =~time + pond)

# fit model

m.time.pond <- mark(toad.proc,
                    toad.ddl, 
                    model.parameters=list(Phi = Phi.time.pond,
                                          f = f.time.pond,
                                          p = p.timexpond))





# Collect models
pema <- collect.models(type = "RDPdfHuggins") #grabs all models in workspace
names(pema) # 4 models
temp <- pema$model.table
( pema.ModSet <- temp[ ,-c(1:5)] )
pema.ModSet

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
#all converged
#AICtableS<- as.data.frame(model.list2$model.table)
AICtableS2<- as.data.frame(pema$model.table)
AICtableS2

AICc.S1B <- model.table(pema, type="RDPdfHuggins", use.lnl = TRUE)
names(AICc.S1B)
taby <- AICc.S1B[ , c(1,2,7,8,9,6,10)]
taby[1,]
class(taby)
names(taby) <- c("Survival","Recruitment","AICc","DAICc","Model Wt","K","-2log(L)")
taby$`Model Wt` <- round(taby$`Model Wt`, 2)
taby$DAICc <- round(taby$DAICc, 2)
taby$AICc <- round(taby$AICc, 2)
taby$`-2log(L)` <- round(taby$`-2log(L)`, 2)
taby


# top model is 'm.timexpond'

# beta coefficients from top model
m.timexpond$results$beta

# real estimates from top model
m.timexpond$results$real

# derived parameters and estimates from top model
m.timexpond$results$derived$`N Population Size`

m.timexpond$results$derived$`Lambda Population Change`

m.timexpond$results$derived$`log(Lambda) Population Change`




# Recall that survival + recruitment = lambda (discrete population growth)


# obtain real estimates for recruitment
# use top model
real=get.real(m.timexpond,"f",se=TRUE)
names(real)
real <- real[,c(3,5,6,11,14)]
rownames(real) <- NULL
rec <- real
rec

# obtain real estimates for apparent survival
# use top model
real=get.real(m.timexpond,"Phi",se=TRUE)
real <- real[,c(3,5,6,11,14)]
rownames(real) <- NULL
surv <- real
surv

# create data frame for plotting
all <- data.frame(pond = rec$pond,
                  year = rec$time,
                  surv = surv$estimate,
                  surv.l = surv$lcl,
                  surv.u = surv$ucl,
                  rec = rec$estimate,
                  rec.l = rec$lcl,
                  rec.u = rec$ucl,
                  lam = m.timexpond$results$derived$`Lambda Population Change`$estimate,
                  lam.l = m.timexpond$results$derived$`Lambda Population Change`$lcl,
                  lam.u = m.timexpond$results$derived$`Lambda Population Change`$ucl)


dat <- data.frame(Site = all$pond,
                  Year=all$year,
                  Variable = "Survival",
                  Value = all$surv,
                  ybegin=all$lam.l,
                  yend=all$lam.u)

dat2 <- data.frame(Site = all$pond,
                   Year=all$year,
                   Variable = "Recruitment",
                   Value = all$rec,
                   ybegin=all$lam.l,
                   yend=all$lam.u)


dat3 <- rbind(dat2, dat)


dat3$Year <- c("16-17","17-18","18-19","19-20",
               "16-17","17-18","18-19","19-20",
               "16-17","17-18","18-19","19-20",
               "16-17","17-18","18-19","19-20",
               "16-17","17-18","18-19","19-20",
               "16-17","17-18","18-19","19-20")


dat3$Site <- c("Pond A","Pond A","Pond A","Pond A",
               "Pond B","Pond B","Pond B","Pond B",
               "Pond C","Pond C","Pond C","Pond C",
               "Pond A","Pond A","Pond A","Pond A",
               "Pond B","Pond B","Pond B","Pond B",
               "Pond C","Pond C","Pond C","Pond C")



# use 'dat3' for plotting

# plot survival, recruitment, and lambda 
plot <- ggplot(dat3, aes(x = Year, y = Value, fill = Variable)) + 
  #geom_bar(stat="identity", width=.5, position = "dodge", colour="black", size = 1) +
  geom_bar(stat = "identity", color = "white", size = 0.1) +
  facet_wrap(~Site,nrow = 1)+
  scale_fill_manual(values = c("grey8","grey55")) +
  geom_errorbar(aes(ymax=ybegin , ymin=yend), width = 0.3, 
                color = "black", size = 0.5) +
  xlab(NULL) +
  ylab(expression(paste("Population Growth Rate (", lambda, ")"))) +
  # ylab(expression(paste(Delta," % Emergent Vegetation"))) +
  scale_y_continuous(limits = c(0, 2.35), breaks = seq(0, 2.2, 0.2)) +
  theme_bw() +
  theme(#text = element_text(size = 18, family = "Times"), # controls all fonts in plot
    #panel.border=element_rect(size=1, color="black"),
    panel.background = element_rect(colour = "black", size=1, linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 20, color = "black"), # spacing from y
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")) +
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(size = 18))+
  theme(strip.text = element_text(size = 18))+
  theme(axis.text.x=element_text(size = 16))+
  theme(axis.text.x=element_text(angle=45, hjust=1))


plot <- plot + theme(legend.position = c(0.15, 0.9))
plot <- plot + theme(legend.key.size = unit(1.5,"line"))
plot <- plot + theme(#legend.justification = "bottom",
  legend.spacing.x = unit(1, "pt"),
  #legend.spacing.y = unit(0, "pt"),
  legend.margin=margin(c(-5,-5,-5,-5)))

plot <- plot + geom_hline(yintercept = 1, linetype=2, size=0.25, col="grey8")

plot

# Reminder:
# in the figure, survival is apparent because we cannot separate mortality from emigration
# and for recruitment, we cannot distinguish between births (in situ rec.) and immigration



# plot abundance at each pond through time
# use estimates from top model
popsize <- m.timexpond$results$derived$`N Population Size`

popsize$Year <- rep(c("2016","2017","2018","2019","2020"), times=3)
popsize$Site <- c(rep("Pond A", times=5),rep("Pond B", times=5),
                  rep("Pond C", times=5))

popsize$Year <- as.factor(popsize$Year)
popsize$Site <- as.factor(popsize$Site)

p <- ggplot(popsize, aes(x=Year, y=estimate, group=Site, color=Site)) +
  geom_point(size=5) + geom_line(size=2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width=.1, size=0.75) +
  xlab(NULL) +
  ylab("Estimated Abundance \n") +
  scale_y_continuous(limits = c(30, 490), breaks = seq(0, 500, 50), expand = c(0,0)) +
  #scale_x_discrete(expand = expand_scale(mult = c(0.45,0.55))) +
  theme_bw() +
  theme(#text = element_text(size = 18, family = "Times"), # controls all fonts in plot
    panel.background = element_rect(colour = "black", size=1, linetype = "solid"),
    #axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
    #axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.2,"cm"),
    axis.title.y = element_text(size = 22, vjust = 1, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"))

p <- p + theme(legend.title = element_blank())

p <- p + theme(legend.position = c(0.75, 0.85))

p <- p + theme(legend.text = element_text(colour="black", size=18))

p







# Let's look at a model for survival with an environmental covariate

# there are various ways to add covariate
# for instance, we can add values directly to the design data
toad.ddl$Phi

toad.ddl$Phi$ppt <- 0

# Pond A 
toad.ddl$Phi$ppt[toad.ddl$Phi$pond == "A" & toad.ddl$Phi$time=="2016"] = 400
toad.ddl$Phi$ppt[toad.ddl$Phi$pond == "A" & toad.ddl$Phi$time=="2017"] = 350
toad.ddl$Phi$ppt[toad.ddl$Phi$pond == "A" & toad.ddl$Phi$time=="2018"] = 450
toad.ddl$Phi$ppt[toad.ddl$Phi$pond == "A" & toad.ddl$Phi$time=="2019"] = 410

# Pond B 
toad.ddl$Phi$ppt[toad.ddl$Phi$pond == "B" & toad.ddl$Phi$time=="2016"] = 500
toad.ddl$Phi$ppt[toad.ddl$Phi$pond == "B" & toad.ddl$Phi$time=="2017"] = 300
toad.ddl$Phi$ppt[toad.ddl$Phi$pond == "B" & toad.ddl$Phi$time=="2018"] = 550
toad.ddl$Phi$ppt[toad.ddl$Phi$pond == "B" & toad.ddl$Phi$time=="2019"] = 405

# Pond C
toad.ddl$Phi$ppt[toad.ddl$Phi$pond == "C" & toad.ddl$Phi$time=="2016"] = 310
toad.ddl$Phi$ppt[toad.ddl$Phi$pond == "C" & toad.ddl$Phi$time=="2017"] = 330
toad.ddl$Phi$ppt[toad.ddl$Phi$pond == "C" & toad.ddl$Phi$time=="2018"] = 290
toad.ddl$Phi$ppt[toad.ddl$Phi$pond == "C" & toad.ddl$Phi$time=="2019"] = 380


toad.ddl$Phi



# Specify model 

# Capture probability
# allow capture probability to vary by year and in each pond
# set capture probability to equal recapture probability
p.timexpond = list(formula=~session * pond, share = TRUE)

# Recruitment 
# allow recruitment to vary over time and at each pond
# i.e., will estimate a different recruitment rate between each year and at each pond
f.timexpond = list(formula=~time * pond)

# Survival 
# model survival as a function of precipitation (ppt)
Phi.ppt = list(formula =~ppt)

# fit model

m.ppt <- mark(toad.proc,
              toad.ddl, 
              model.parameters=list(Phi = Phi.ppt,
                                    f = f.timexpond,
                                    p = p.timexpond))


# beta coefficients 
m.ppt$results$beta

# real estimates from top model
m.ppt$results$real


# plot
vals=get.real(m.ppt,"Phi",se=TRUE)
real.values=vals
real <- real.values
names(real)
real <- real[,c(3,5,6,11,14)]
real

# add covariate data
real$ppt <- toad.ddl$Phi$ppt

real <- real[order(real$ppt),]

# plot
par(mar=c(7,6,5,1))
plot(real$ppt, real$estimate,type="l",
     xlab="Annual Precipitation (mm)",
     ylab="Survival Probability", ylim = c(0.25,0.9), lwd=4, col="black",
     cex.lab=1.5, cex.axis=1.25)
box(lwd=2)
lines(real$ppt,real$lcl,lty=2, lwd=2, col="blue")
lines(real$ppt,real$ucl,lty=2, lwd=2, col="blue")

























