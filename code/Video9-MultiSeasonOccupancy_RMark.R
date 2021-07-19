
# Code to fit 'Multi-Season Occupancy Model' in the 'RMark' package

# Addressing the question: 
# How does little tern breeding occurrence change over time on beaches in northern Europe?
  

# Code last updated on 7/15/2021 by Gabe Barrile



# install the packages that we will use:

#install.packages("RMark")
#install.packages("reshape")
#install.packages("ggplot2")

# make sure packages are installed and active
require(RMark)
require(reshape)
require(ggplot2)


# citing the unmarked package
citation("RMark")


# read-in the little tern detection/nondetection data from the beaches that we surveyed 
# read-in data from the csv
df <- read.csv("data/LittleTern_Breeding.csv")

# check out data
head(df)

# how many beaches?
length(unique(df$Site)) # 8 beaches

# format data into detection histories for each site (i.e., each beach)
m <- melt(df,id.var=c("Site","Survey"),measure.var="Eggs")

y = cast(m, Site ~ Survey)

head(y)

#y[is.na(y)]=0
k<-dim(y)[2]
#if needing to deal with repeats on 1 occasion
#y[,2:k]<-(y[,2:k]>0)*1
#head(y)

#function to create detection history strings
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

#detection history data frame 
det.hist <- data.frame(ch=pasty(y[,2:k]))
head(det.hist)

y$ch <- pasty(y[,2:k])

head(y)


# tern will be our dataframe with which to proceed
tern <- data.frame(ch = y$ch, freq = 1, beach = y$Site)

# make beach a factor variable
tern$beach <- as.factor(as.character(tern$beach))

# remove unneeded objects
rm(m,det.hist,y,k,pasty)

# are we missing any data?
table(is.na(tern)) # no missing data = good!

# here is our data frame
tern

#specify time intervals for Robust Design
intervals  <-  c(0,0,1, # 2015
                 0,0,1, # 2016
                 0,0,1, # 2017
                 0,0,1, # 2018
                 0,0)   # 2019

# process data in RMark format
t.proc = process.data(tern, model = "RDOccupEG", groups = "beach", 
                      time.intervals=intervals, begin.time = 2015)

# create design data
t.ddl = make.design.data(t.proc)

# look at design data
names(t.ddl)

t.ddl$Psi # initial occupancy
# probability that little terns bred on a given beach in 2015

t.ddl$Epsilon # extinction
# probability of an occupied site becoming unoccupied
# probability of a beach that little terns bred on during time t
# was not bred on during time t + 1

t.ddl$Gamma # colonization
# probability of an unoccupied site becoming occupied
# probability of a beach that little terns did not breed on during time t
# being bred on during time t + 1

t.ddl$p # detection probability




# Specify models (p, Psi, Epsilon, and Gamma)

# define p model
t.ddl$p # detection probability
# let's obtain a different detection probability in each year
p.1 = list(formula =  ~  session)
# note that, for the robust design, 'session' indicates each year
# and 'time' indicates surveys within each year for detection probability
# this can be confusing because 'time' indicates each year for the
# probability of occupancy, extinction, and colonization


# define Psi model
t.ddl$Psi # initial occupancy
# let's allow initial occupancy to vary across beaches
Psi.1 = list(formula =  ~ beach)


# define Epsilon model
t.ddl$Epsilon # extinction probability
# let's allow extinction probability to vary over time
E.1 = list(formula =  ~ time)


# define Gamma model
t.ddl$Gamma # colonization probability
# let's allow colonization probability to vary over time
G.1 = list(formula =  ~ time)




# fit model
setwd("F:/Multistate/models")
# 1
occ <- mark(t.proc, # processed data
            t.ddl,  # design data
            model.parameters=list(Psi = Psi.1,   # initial occupancy
                                  Gamma = G.1,   # colonization
                                  Epsilon = E.1, # extinciton
                                  p = p.1))      # detection


# check out model output

# beta coefficients
occ$results$beta

# real estimates
occ$results$real

# derived parameters
occ$results$derived$`psi Probability Occupied`


# let's visualize the results

# let's plot temporal variation in extinction probability over our study period
ext <- occ$results$real[9:12,c(1,3,4)]

ext$year <- c("2015-2016","2016-2017","2017-2018","2018-2019")
ext$year <- as.factor(ext$year)
ext

# plot the values
ggplot(ext, aes(x=year, y=estimate)) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=0.2, size=2, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), stat="identity", size = 7) +
  xlab(NULL) +
  ylab("Extinction probability \n") +
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


# let's plot temporal variation in colonization probability over our study period
col <- occ$results$real[13:16,c(1,3,4)]

col$year <- c("2015-2016","2016-2017","2017-2018","2018-2019")
col$year <- as.factor(col$year)
col

# plot the values
ggplot(col, aes(x=year, y=estimate)) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=0.2, size=2, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), stat="identity", size = 7) +
  xlab(NULL) +
  ylab("Colonization probability \n") +
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



# let's plot temporal variation in colonization and extinction probability 
# over our study period on one graph

ext$colext <- "Extinction"
col$colext <- "Colonization"

colext <- rbind(ext, col)

colext$colext <- as.factor(colext$colext)


# plot the values
ggplot(colext, aes(x=year, y=estimate, group=colext, color=colext)) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=0.2, size=2, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), stat="identity", size = 7) +
  xlab(NULL) +
  ylab("Probability of change in breeding \n") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw()+
  theme(
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








# let's plot temporal variation in detection probability over our study period
det <- occ$results$real[17:21,c(1,3,4)]

det$year <- c("2015","2016","2017","2018","2019")
det$year <- as.factor(det$year)
det

# plot the values
ggplot(det, aes(x=year, y=estimate)) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=0.2, size=2, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), stat="identity", size = 7) +
  xlab(NULL) +
  ylab("Detection probability \n") +
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



# let's plot temporal variation in occupancy probability at each beach
# over our study period
occup <- occ$results$derived$`psi Probability Occupied`

occup$year <- c("2015","2016","2017","2018","2019",
                "2015","2016","2017","2018","2019",
                "2015","2016","2017","2018","2019",
                "2015","2016","2017","2018","2019",
                "2015","2016","2017","2018","2019",
                "2015","2016","2017","2018","2019",
                "2015","2016","2017","2018","2019",
                "2015","2016","2017","2018","2019")
occup$year <- as.factor(occup$year)
occup

occup$Beach <- c("1","1","1","1","1",
                 "2","2","2","2","2",
                 "3","3","3","3","3",
                 "4","4","4","4","4",
                 "5","5","5","5","5",
                 "6","6","6","6","6",
                 "7","7","7","7","7",
                 "8","8","8","8","8")

occup$Beach <- as.factor(occup$Beach)
occup


# plot the values
# plot the values
ggplot(occup, aes(x=year, y=estimate, group=Beach, color=Beach)) +
  #geom_point(size=5) + geom_line(size=1.5) +
  geom_point(position=position_dodge(0.3), stat="identity", size = 5) +
  geom_line(position=position_dodge(0.3), size=1.5) +
  xlab(NULL) +
  ylab("Occupancy probability \n") +
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = seq(0, 1, 0.2), expand = c(0,0)) +
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
  theme(legend.title = element_text(size = 17))+
  theme(legend.text = element_text(size = 17))+
  theme(legend.position="top")+
  theme(strip.text = element_text(size = 7))



