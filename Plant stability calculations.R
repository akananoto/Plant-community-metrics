library(vegan)
library(plyr)
library(doBy) #for summaryby which seems to do the same thing as aggregate, as far as I can tell
library(cgwtools)#for seqle

plant<-read.csv("~/Documents/PhD docs/2012-2013/Environmental variability/Plant/all_plant.csv")
plant<-plant[,-c(11)]

#Find first and last year of data for all sites
plant.year<-summaryBy(Year~Location, data=plant, FUN=c(max,min))
write.csv(plant.year,"C:/Users/Akana/Documents/2012-2013/Environmental variability/Plant/plant_year.csv",row.names=FALSE)

#Need to figure out how many years of data there are for each experiment/plot
test<-unique(plant[,c(1,3)])
testcount<-count(test, vars=c("sitesubplot"))
plant<-merge(plant,testcount,by=c("sitesubplot"))
#plant$Abundance<-as.character(plant$Abundance)
#plant$Abundance<-as.numeric(plant$Abundance)

#Better to know how many years of consecutive data there are for each experiment/plot 
#test$sitesubplot<-paste(test$Experiment,test$Plot,test$Subplot,test$Location, sep="-")
sitenames<-sort(unique(test$sitesubplot))

output<-data.frame(cbind(Year=numeric(0),sitesubplot=character(0),length=numeric(0)))

for (i in 1:length(sitenames)){
  sub<-subset(test, test$sitesubplot==sitenames[i])
  sub<-sub[order(sub[,"Year"]),]
  
  result<-seqle(unique(sub$Year))
  int<-as.data.frame(result$length)
  int$Year<-result$values
  colnames(int)<-c("length","Year")
  
  temp<-merge(sub, int, by="Year",all=TRUE)
  
  output=rbind(output, temp) 
}

for(j in 1:nrow(output)){
  k<-j-1
  if (is.na(output$length[j])) output$length[j]<-output$length[k] 
  else output$length[j]
}

#Let's remove all plots for which there aren't at least five consecutive years of data 
cutout<-subset(output, length>4)
cutplants<-merge(cutout, plant, by=c("Year","sitesubplot"))
#Though this could still have more than one stretch for the same plot

# Calculate species richness for each experiment for each year ------------
#I think I don't need consecutive years for this, so I'm just using original plant data
#I have included SO here even though not enough consecutive years or permanent plots for a lot of the years
#Seems pretty uniform--should I be doing this for each plot for each year and then averaging those? I should probably be doing it for subplots (or whatever I called the 1m x 1m units) so that area is consistent across all locations
dummy<-plant[which(plant$Abundance>0),]
dummy$count<-1
rich<-aggregate(count~Location+Experiment+Year+Plot+Subplot, data=dummy, FUN=sum)
colnames(rich)<-c("Location","Experiment","Year","Plot","Subplot","Richness")
a<-summaryBy(Richness~Location, data=rich) #To figure out average richness in all subplots

#Do I want to calculate these things over the whole site or for each little area in each site? I should decide about that.
rich$mean<-ave(rich$Richness, rich$Location, rich$Experiment, FUN=mean)
rich$sd<-ave(rich$Richness, rich$Location, rich$Experiment, FUN=sd)
rich$cv<-rich$sd/rich$mean
#Let's start by doing for the whole Location lumping together experiments
richsum<-unique(rich[,c(1,4:7)])

averich<-aggregate(rich$Richness,by=c(rich["Location"]), FUN=mean)
colnames(averich)<-c("Location","ave.rich")
write.csv(averich, "C:/Users/Akana/Documents/2012-2013/Environmental variability/Plant/ave_sprich.csv", row.names=FALSE)

#I want to calculate mean species richness for each location, but I need to scale according to plot size (and number of plots?)
#I'm doing species accumulation curve to scale for number of plots--what can I do to scale for plot size? Maybe I could do this per area sampled? The main problem being that area is continuous
#Would there be a reason to do rarefaction? Maybe since there are different sampling efforts and some of these don't seem like they've leveled off, would be good idea
#Note: Should I be doing a different curve (e.g. mess with method, permutations?)
plant1<-plant[,c(2:10,12)]
plant.wide<-reshape(plant1, timevar="Species",idvar=c("Experiment","Plot","Subplot","Year","Location","unitAbund","Coast","freq"),direction="wide")

sites<-unique(plant.wide$Location)
output<-data.frame(cbind(sites=numeric(0),richness=numeric(0),sd=numeric(0),Location=character(0)))

for(i in 1:length(sites)){
  sub<-subset(plant.wide, Location==sites[i])
  sub<-sub[,-c(1:8)]
  sub[is.na(sub)] <- 0
  
  sp.sub<-specaccum(sub, method="random", permutations=1000)
  temp<-with(sp.sub, data.frame(sites, richness, sd))
  temp$Location<-sites[i]
  
  output<-rbind(output, temp)
}

#Plot sp richness for each number of plots for each location to see what it looks like
#Seems like SO and LPL have super-high richness. Are they including more different elevations? 
library(ggplot2)

qplot(x=output$sites, y=output$richness, 
      data=output, 
      colour=Location, 
      main="GGPLOT line plot with groups") +
  geom_line()

#Seems like some sites don't have enough sampling to have leveled out, even some with quite a lot of sampling. Could get rid of WEL maybe since it's quite unleveled... Based on this, I could choose some number of plots that will be good for all (i.e. smallest # of samples, see Gotelli and Colwell 2011 for reference) for which we can get a rarefied S, or we can use the max number since most are pretty close to leveling out?
calib<-subset(output, sites==55)
write.csv(calib, "C:/Users/Akana/Documents/2012-2013/Environmental variability/Plant/rare_sprich.csv", row.names=FALSE)

#Let's try the same thing again with Location and year
plant.wide<-reshape(plant1, timevar="Species",idvar=c("Experiment","Plot","Subplot","Year","Location","unitAbund","Coast","freq"),direction="wide")
plant.wide$Loc.year<-paste(plant.wide$Location,plant.wide$Year,sep=".")

sites<-unique(plant.wide$Loc.year)
output<-data.frame(cbind(sites=numeric(0),richness=numeric(0),sd=numeric(0),Location=character(0),year=character(0)))

for(i in 1:length(sites)){
  sub<-subset(plant.wide, Loc.year==sites[i])
  subber<-sub[,-c(1:8,163)]
  subber[is.na(subber)] <- 0
  
  sp.sub<-specaccum(subber, method="random", permutations=1000)
  temp<-with(sp.sub, data.frame(sites, richness, sd))
  temp$Loc.year<-sites[i]
  temp$Location<-sub$Location[1]
  temp$Year<-sub$Year[1]
  
  output<-rbind(output, temp)
}

#Plot sp richness for each number of plots for each location to see what it looks like
#Seems like SO and LPL have super-high richness. Are they including more different elevations? 
library(ggplot2)

qplot(x=output$sites, y=output$richness, 
      data=output, 
      colour=Location,
      main="GGPLOT line plot with groups") +
  geom_line()

#Have to do a lot fewer sites to make this work with the yearly data since fewer observations
calib<-subset(output, sites==20)

#Or figure some # of sites for which it makes sense to use sp richness as mean sp richness over time
#Something to check is whether some of these locations tend to have higher marsh plots cause that might cause higher species richness (e.g. SO, LPL)

#To rarefy with these data, use function rarefaction.sample from package rareNMtests. Only problem is that data need to be converted to presence (1) absence (0) with species in columns and sampling unit in row
library(rareNMtests)
#We already have plant.wide, so we just need to reformat it slightly
str(plant.wide)
plant.wide.r<-plant.wide[,-c(1,8:9)]

rdummy<-dummy[,-c(7,10:12)]
plant.wide.r<-reshape(rdummy, timevar="Species",idvar=c("Experiment","Plot","Subplot","Year","Location","unitAbund","sitesubplot"),direction="wide")
plant.wide.r[is.na(plant.wide.r)]<-0

sbr.q0 <- rarefaction.sample(Chiapas[,-1])
#Now need to do rarefaction for each site separately
sites<-unique(plant.wide.r$Location)
output<-data.frame(cbind(samplesize=numeric(0),rich=numeric(0),Location=character(0)))

for(i in 1:length(sites)){
  sub<-subset(plant.wide.r, Location==sites[i])
  temp<-rarefaction.sample(sub[,-c(1:7)])
  temp$Location<-sites[i]
  
  output<-rbind(output, temp)
}
colnames(output)<-c("samplesize","richness","Location")

library(ggplot2)

qplot(x=output$samplesize, y=output$richness, 
      data=output, 
      colour=Location, 
      main="GGPLOT line plot with groups") +
  geom_line()
#Oh, I just realized that specaccum if you do the randomized way is the same as sample-based rarefaction... So at what point do I take richnesses to do analysis with environment?


# Try rarefying with # of plots and plot area ------------------
#Same method as above
#Don't need consecutive years, so using original plant data
#Seems pretty uniform--should I be doing this for each plot for each year and then averaging those? I should probably be doing it for subplots (or whatever I called the 1m x 1m units) so that area is consistent across all locations
size<-read.table("~/Documents/PhD Docs/2012-2013/Environmental variability/Plot size.csv", header=T, sep=",")
plant2<-merge(plant,size, by="Location")
dummy<-plant2[which(plant2$Abundance>0),]
dummy$count<-1
rich2<-aggregate(count~Location+Experiment+Year+Plot+Subplot+Plot.size, data=dummy, FUN=sum)
colnames(rich2)<-c("Location","Experiment","Year","Plot","Subplot","Plot.size","Richness")
#scale per 1m2
rich2$Rich.scale<-rich2$Richness/rich2$Plot.size
a<-summaryBy(Rich.scale~Location, data=rich2) #To figure out average richness in all subplots

#I want to calculate mean species richness for each location, but I need to scale according to plot size and number of plots
#I'm doing species accumulation curve to scale for number of plots, but really I want to do it for area sampled, and would be better to have real area sampled instead of scaled. Can try the weighting thing in specaccum function, although all weights within sites are the same
plant3<-plant2[,c(1,3:10,12,13)]
plant2.wide<-reshape(plant3, timevar="Species",idvar=c("Experiment","Plot","Subplot","Year","Location","unitAbund","Coast","freq","Plot.size"),direction="wide")

sites<-unique(plant2.wide$Location)
output<-data.frame(cbind(sites=numeric(0),richness=numeric(0),sd=numeric(0),Location=character(0)))

for(i in 1:length(sites)){
  sub<-subset(plant2.wide, Location==sites[i])
  sub<-sub[,-c(1:8)]
  sub[is.na(sub)] <- 0
  
  sp.sub<-specaccum(sub, method="random", permutations=1000)
  temp<-with(sp.sub, data.frame(sites, richness, sd))
  temp$Location<-sites[i]
  
  output<-rbind(output, temp)
}

#But now to make area sampled, I'm going to multiply samplesize by plot area to make x-axis area sampled
out2<-merge(output,size,by="Location")
out2$SampArea<-out2$samplesize*out2$Plot.size

#Now have to find smallest area that works for all
a<-summaryBy(SampArea~Location, data=out2, FUN=max)
#Smallest area is PIE at 55m2
calib<-subset(out2, SampArea<56)
qplot(x=calib$SampArea, y=calib$richness, 
      data=calib, 
      colour=Location,
      main="GGPLOT line plot with groups") +
  geom_line()
#Actually things look surprisingly leveled out
calib2<-summaryBy(richness~Location, data=calib, FUN=max)
write.csv(calib2, "~/Documents/PhD Docs/2012-2013/Environmental variability/Plant/rare_sprich_byarea.csv", row.names=FALSE)

# Calculate CV and Tilman stability for each species over time ---------------------------------
#First, calculate mean and sd of cover for each species for each plot over all years for which we have data. 
#Note that these means and SDs are calculated for the period of record for each subplot, with no scaling for length of time. 
#Not sure what the best way to scale it would be. I guess I could choose a 5-year period from within each longer period, but that seems like kind of a waste of data. . .
#Things to figure out here: 1. Do I need to scale data in some way? What's the best way? 2. Does it make more sense to calc stats for each species across all years and then find mean for each location or to calc ave of all species for each little plot and then for whole location? 3. Is there a better way to do whole community Tilman stability? Check out code from Hao's friend.
cutplants$sitesubplot<-paste(cutplants$Experiment, cutplants$Plot, cutplants$Subplot, sep="_")
cutplants$mean<-ave(cutplants$Abundance, cutplants$Species, cutplants$Location, cutplants$sitesubplot, FUN=mean)
cutplants$sd<-ave(cutplants$Abundance, cutplants$Species, cutplants$Location, cutplants$sitesubplot, FUN=sd)

#Use that to calculate CV and Tilman stability for each species in each plot across all years
cutplants$cv<-ifelse(cutplants$mean==0, NA, cutplants$sd/cutplants$mean)
cutplants$Tilman<-ifelse(cutplants$mean==0, NA, cutplants$mean/cutplants$sd)

#To get average stability for each species in each site ("Location")
#Need everything only by Species, Location, Coast, unitAbund, avecv, aveTil
avecv<-aggregate(cv~Species+Location, data=cutplants, mean)
colnames(avecv)<-c("Species","Location","avecv")
aveTil<-aggregate(Tilman~Species+Location, data=cutplants, mean)
colnames(aveTil)<-c("Species","Location","aveTil")
cutplants<-merge(cutplants,avecv,by=c("Species","Location"))
cutplants<-merge(cutplants,aveTil,by=c("Species","Location"))
var<-cutplants[,c(1,2,7,9:11,17:18)]
var<-unique(var)

var2<-unique(var[,c(1:4,6:7)])

#To get mean stability for all species in each location 
avestab<-summaryBy(aveTil+avecv~Location+Coast+unitAbund, data=var, FUN=mean)

#To get mean stability for all species in each subplot and then for each location
avesub<-summaryBy(cv+Tilman~sitesubplot+Location+Coast+unitAbund, data=cutplants, FUN=c(mean, sd))
avestab.loc<-summaryBy(cv.mean+Tilman.mean~Location+Coast+unitAbund, data=avesub, FUN=mean)


# Bray-Curtis dissimilarities ---------------------------------------------
#Modify code from Elsa's group
#Calculating this for each subplot individually between each pair of years and then averaging over the whole site
#Not how Elsa did it, so just make sure I shouldn't rethink that strategy

#first create a lag so there is year and year+1 in each column
cutplants$yrminus<-cutplants$Year + 1
yrminus<-cutplants[, c("Experiment","sitesubplot","Location","yrminus","Species","Abundance")]
year<-cutplants[, c("Experiment","sitesubplot","Location","Year","Species","Abundance")]

#get rid of pesky rownames that mess with merge
#row.names(yrminus) <- NULL 
#row.names(year) <- NULL

#rename variables 
yrminus2 <- rename(yrminus, c(yrminus="Year", Abundance="abundprioryr"))

#merge & calculate abs diff between abund in year and prior yr
allyears<-merge(year, yrminus2, by=c("Experiment","sitesubplot","Location","Species", "Year"))
allyears$diff<-abs(allyears$Abundance-allyears$abundprioryr)
allyears$total<-allyears$Abundance+allyears$abundprioryr

#to calculate brays first sum the differences for all species in each, then divide by sum of total cover in each (if they summed to 100, divide by 200)
#What about for ones where it sums to less than 100 because of bare ground?
yrsum<-summaryBy(diff + total ~ Experiment + sitesubplot+ Location + Year, data=allyears, FUN=sum)
yrsum$braycurtis<-yrsum$diff.sum/yrsum$total.sum
yrsum[610,"braycurtis"]<-0
#Some NAs here from CAR sites that they seemingly stopped collecting data at
#Bigger problem is the one CRMS site that has 0/0=NaN. Plot has zero cover in 2009 and 2010. I just made it 0 for identical

sitesum<-summaryBy(braycurtis~Location, data=yrsum, FUN=c(mean,sd), na.rm=TRUE)
levels(sitesum$Location)[levels(sitesum$Location)=="GCE"]<-"SAP"

write.csv(sitesum, "C:/Users/Akana/Documents/2012-2013/Environmental variability/Plant/braycurtis.csv")

#Seems like the one year thing misses directional change. Let's try a larger time difference
#Note that there are fewer observations as we increase the length of the time period
cutplants$yrminus<-cutplants$Year + 2
yrminus<-cutplants[, c("Experiment","sitesubplot","Location","yrminus","Species","Abundance")]
year<-cutplants[, c("Experiment","sitesubplot","Location","Year","Species","Abundance")]

#get rid of pesky rownames that mess with merge
#row.names(yrminus) <- NULL 
#row.names(year) <- NULL

#rename variables 
yrminus2 <- rename(yrminus, c(yrminus="Year", Abundance="abundprioryr"))

#merge & calculate abs diff between abund in year and prior yr
allyears<-merge(year, yrminus2, by=c("Experiment","sitesubplot","Location","Species", "Year"))
allyears$diff<-abs(allyears$Abundance-allyears$abundprioryr)
allyears$total<-allyears$Abundance+allyears$abundprioryr

#to calculate brays first sum the differences for all species in each, then divide by sum of total cover in each (if they summed to 100, divide by 200)
#What about for ones where it sums to less than 100 because of bare ground?
yrsum<-summaryBy(diff + total ~ Experiment + sitesubplot+ Location + Year, data=allyears, FUN=sum)
yrsum$braycurtis<-yrsum$diff.sum/yrsum$total.sum
#Some NAs here from CAR sites that they seemingly stopped collecting data at
#Bigger problem is the one CRMS site that has 0/0=NaN. Plot has zero cover in 2009 and 2010. I just made it 0 for identical

sitesum2<-summaryBy(braycurtis~Location, data=yrsum, FUN=c(mean,sd), na.rm=TRUE)
levels(sitesum2$Location)[levels(sitesum2$Location)=="GCE"]<-"SAP"
colnames(sitesum2)<-c("Location","braycurtis.mean2","braycurtis.sd2")
allbray<-merge(sitesum, sitesum2, by="Location")

#3 year intervals
cutplants$yrminus<-cutplants$Year + 3
yrminus<-cutplants[, c("Experiment","sitesubplot","Location","yrminus","Species","Abundance")]
year<-cutplants[, c("Experiment","sitesubplot","Location","Year","Species","Abundance")]

yrminus2 <- rename(yrminus, c(yrminus="Year", Abundance="abundprioryr"))

allyears<-merge(year, yrminus2, by=c("Experiment","sitesubplot","Location","Species", "Year"))
allyears$diff<-abs(allyears$Abundance-allyears$abundprioryr)
allyears$total<-allyears$Abundance+allyears$abundprioryr

yrsum<-summaryBy(diff + total ~ Experiment + sitesubplot+ Location + Year, data=allyears, FUN=sum)
yrsum$braycurtis<-yrsum$diff.sum/yrsum$total.sum

sitesum3<-summaryBy(braycurtis~Location, data=yrsum, FUN=c(mean,sd), na.rm=TRUE)
levels(sitesum3$Location)[levels(sitesum3$Location)=="GCE"]<-"SAP"
colnames(sitesum3)<-c("Location","braycurtis.mean3","braycurtis.sd3")
allbray<-merge(allbray, sitesum3, by="Location")

#Let's also just try 4 year difference
cutplants$yrminus<-cutplants$Year + 4
yrminus<-cutplants[, c("Experiment","sitesubplot","Location","yrminus","Species","Abundance")]
year<-cutplants[, c("Experiment","sitesubplot","Location","Year","Species","Abundance")]

yrminus2 <- rename(yrminus, c(yrminus="Year", Abundance="abundprioryr"))

allyears<-merge(year, yrminus2, by=c("Experiment","sitesubplot","Location","Species", "Year"))
allyears$diff<-abs(allyears$Abundance-allyears$abundprioryr)
allyears$total<-allyears$Abundance+allyears$abundprioryr

yrsum<-summaryBy(diff + total ~ Experiment + sitesubplot+ Location + Year, data=allyears, FUN=sum)
yrsum$braycurtis<-yrsum$diff.sum/yrsum$total.sum

sitesum4<-summaryBy(braycurtis~Location, data=yrsum, FUN=c(mean,sd), na.rm=TRUE)
levels(sitesum4$Location)[levels(sitesum4$Location)=="GCE"]<-"SAP"
colnames(sitesum4)<-c("Location","braycurtis.mean4","braycurtis.sd4")
allbray<-merge(allbray, sitesum4, by="Location")

write.csv(allbray, "C:/Users/Akana/Documents/2012-2013/Environmental variability/Plant/allbraycurtis.csv")


# Bray-Curtis at site level only ------------------------------------------
#Calculating this at site level per Elsa's recommendation
#Since on site level, I didn't worry about whether there were data for every year in every plot. I think that makes sense...(allows me to keep sites where they switched plots every year)

#first create a lag so there is year and year+1 in each column
plant$yrminus<-plant$Year + 1
#maybe I need to sum everything up first before merging (sum all abundances of the same species in the same location and year)
plant2<-summaryBy(Abundance~Species+Location+Year+yrminus,data=plant, FUN=sum)
#This data frame still doesn't contain zeroes for years in which the species isn't there
yrminus<-plant2[, c("Location","yrminus","Species","Abundance.sum")]
year<-plant2[, c("Location","Year","Species","Abundance.sum")]

#rename variables 
yrminus2 <- rename(yrminus, c(yrminus="Year", Abundance.sum="abundprioryr"))

#merge & calculate abs diff between abund in year and prior yr
#so that species abundances are zero for years before or after species have come in or out, merge all=TRUE
allyears<-merge(year, yrminus2, by=c("Location","Species", "Year"),all=TRUE)
allyears[is.na(allyears)]<-0
allyears$diff<-abs(allyears$Abundance.sum-allyears$abundprioryr)
allyears$total<-allyears$Abundance.sum+allyears$abundprioryr

#to calculate brays first sum the differences for all species in each, then divide by sum of total cover in each (if they summed to 100, divide by 200)
#What about for ones where it sums to less than 100 because of bare ground?
yrsum<-summaryBy(diff + total ~ Location + Year, data=allyears, FUN=sum, na.rm=TRUE)
yrsum$braycurtis<-yrsum$diff.sum/yrsum$total.sum
#Some NAs here from CAR sites that they seemingly stopped collecting data at
#I just omitted those NAs since we are on a site level anyway, which I think is okay, but is it?

sitesum<-summaryBy(braycurtis~Location, data=yrsum, FUN=c(mean,sd))
sitesum$Location<-as.character(sitesum$Location)
sitesum$Location[sitesum$Location=="GCE"]<-"SAP"
sitesum$Location<-as.factor(sitesum$Location)

write.csv(sitesum, "~/Documents/PhD docs/2012-2013/Environmental variability/Plant/braycurtis_sitelevel.csv",row.names=FALSE)


# Tilman community stability - Aggregate species abundances (per Hallett et al 2014) --------
#This method only works when species abundances don't sum to 100
#Aggregate species abundances within replicate and year, then calculate mean and sd for each replicate to get stability
#Better to do this by replicate and then take the means or to calculate mean and sd based on all replicates in each site? The latter is more consistent with what I've been doing for the other measures, 
#Here also let's only do for plots with at least 5 consecutive years of data
#But actually, we probably need to do the with the moving window
unbound<-subset(plant, bounded=="no")
agg_abund<-aggregate(Abundance~Location+Year, data=unbound, FUN=sum)
stab_agg_abund<-aggregate(Abundance~Location, data=agg_abund, FUN=mean)
colnames(stab_agg_abund)<-c("Location","mean")
stab_agg_abund.sd<-aggregate(Abundance~Location, data=agg_abund, FUN=sd)
colnames(stab_agg_abund.sd)<-c("Location","sd")
stab_agg_abund<-merge(stab_agg_abund, stab_agg_abund.sd,by="Location")
stab_agg_abund$stability<-(stab_agg_abund$mean)/(stab_agg_abund$sd)
library(data.table)
out<-data.frame(cbind(Location=character(0),aggstab=numeric(0)))
sites<-unique(unbound$Location)
for(i in 1:length(sites)){
  sub<-subset(unbound,Location==sites[i])
  #Make sure ordered by year
  sub<-sub[order(sub$Year),]
  #Add index
  sub<-data.table(sub)
  sub[,index:=.GRP,by=Year]
  
  index<-unique(sub$index)
  for (j in 1:(length(index)-5)){
    subber<-subset(sub,j+5>=index&index>=j)
    
    #Calc Tilman stability
    agg_abund<-aggregate(Abundance~Year, data=subber, FUN=sum)
    agg_abund$mean<-mean(agg_abund$Abundance,trim=0)
    agg_abund$sd<-sd(agg_abund$Abundance)
    
    temp<-as.data.frame(unique(subber$Location))
    temp$stability<-unique((agg_abund$mean)/(agg_abund$sd))
    
    out<-rbind(out,temp)
  }
}
colnames(out)<-c("Site","aggstab")
out$Site<-as.character(out$Site)
out$Site[out$Site=="GCE"]<-"SAP"
out$Site<-as.factor(out$Site)

write.csv(out, "~/Documents/PhD docs/2012-2013/Environmental variability/Plant/site_aggstab_moving.csv", row.names=FALSE)

# Dominant species population stability (using Tilman for a since  --------
#First identify the dominant species in each replicate over time (highest mean abundance)
#This is above, but if it needs to be run independently:
#cutplants$mean<-ave(cutplants$Abundance, cutplants$Species, cutplants$Location, cutplants$sitesubplot, FUN=mean)
#cutplants$sd<-ave(cutplants$Abundance, cutplants$Species, cutplants$Location, cutplants$sitesubplot, FUN=sd)
#Pick out species from each plot with max mean value, then subset stabilities calculated above based on those species
cutplants$Lsitesubplot<-paste(cutplants$Location, cutplants$Experiment, cutplants$Plot, cutplants$Subplot, sep="_")
cutplants_short<-cutplants[,-c(1:4, 8:12)]
cutplants_short<-unique(cutplants_short)
domstab<-data.frame(cbind(Location = character(0), length=character(0), Species=character(0), mean = numeric(0), sd=numeric(0), Lsitesubplot=character(0)))

reps<-unique(cutplants_short$Lsitesubplot)

for (i in 1:length(reps)){
  rep<-reps[i]
  sub<-subset(cutplants_short, Lsitesubplot==rep)
  
  temp<-subset(sub, sub$mean==max(sub$mean))
  domstab<-rbind(temp, domstab)
}

#Now calculate stability of that species in each plot and find mean for each site
domstab$domstab<-domstab$mean/domstab$sd
domstab_agg<-aggregate(domstab$domstab, by=c(domstab["Location"]), FUN=mean)

#What to do about infinite values? (i.e. where sd=0)
#Could try calculating stabilities based on mean cover over all the plots, but that will only work if the dominant species is the same in all. Alternatively, could omit those values or could substitute with something else small



# Proper Tilman stability -------------------------------------------------
#Calculate using code from Adam for each plot and then average for site, I guess
#sum of biomass per plot
til<-read.csv("C:/Users/Akana/Documents/2012-2013/Environmental variability/Plant/all_plant_long.csv")
til$Lsitesubplot<-paste(til$Location, til$sitesubplot, sep="-")
#til<-subset(til, Location!="SO")
til<-til[,c("Year","Lsitesubplot","Species","Abundance")]

#Here we should also make sure to get rid of plots for which there isn't consecutive data 
tilcheck<-unique(til[c("Year","Lsitesubplot")])

sitenames<-sort(unique(tilcheck$Lsitesubplot))

output<-data.frame(cbind(Year=numeric(0),Lsitesubplot=character(0),length=numeric(0)))

for (i in 1:length(sitenames)){
  sub<-subset(tilcheck, tilcheck$Lsitesubplot==sitenames[i])
  sub<-sub[order(sub[,"Year"]),]
  
  result<-seqle(unique(sub$Year))
  int<-as.data.frame(result$length)
  int$Year<-result$values
  colnames(int)<-c("length","Year")
  
  temp<-merge(sub, int, by="Year",all=TRUE)
  
  output=rbind(output, temp) 
}

for(j in 1:nrow(output)){
  k<-j-1
  if (is.na(output$length[j])) output$length[j]<-output$length[k] 
  else output$length[j]
}

#Let's remove all plots for which there aren't at least five consecutive years of data 
til<-merge(output, til, by=c("Year","Lsitesubplot"))
cuttil<-subset(til, length>4)

write.csv(cuttil, "C:/Users/Akana/Documents/2012-2013/Environmental variability/Plant/raw_abund_4+_years.csv", row.names=FALSE)

temp<-as.data.frame(unique(cuttil$Lsitesubplot))
temp$plot<-seq(1:410)
colnames(temp)<-c("Lsitesubplot","plot")
cuttil<-merge(cuttil,temp, by="Lsitesubplot")

sum_abund<-rowSums(tapply(cuttil$Abundance, list(cuttil$plot, cuttil$Species), mean))
#sum of variance for each species in each plot
sum_var<-rowSums(tapply(cuttil$Abundance, list(cuttil$plot, cuttil$Species), var))

#NAs are a problem for this code. Why do I have them?
find_covar_sum<-function(x) {
  nspecies<-length(unique(x$Species))
  nsamples<-length(unique(x$Year))
  
  comparisons<-outer(1:nspecies, 1:nspecies, ">")
  
  #Make a matrix with species biomass for the plot in each column, and rows ordered by sampling event
  sp_mat<-matrix(nrow=nsamples, ncol=nspecies, data=x$Abundance[order(paste(x$Species, x$Year))])
  
  return(var(sp_mat)[comparisons]) #Record only unique, non-diagonal values 
}

#Run function for each plot
nsites<-length(unique(cuttil$plot))
sum_covar<-numeric(nsites)
n<-1
for(i in sort(unique(cuttil$plot))) {
  x<-cuttil[cuttil$plot==i,]
  sum_covar[n]<-sum(find_covar_sum(x))
  names(sum_covar)[n]<-i
  n<-n+1
}

#nsites<-length(unique(cuttil$plot))
#sum_covar2<-data.frame(cbind(plot=numeric(0),sum_covar=numeric(0),nrow=numeric(0),datalength=numeric(0),ratio=numeric(0)))
#for(i in sort(unique(cuttil$plot))) {
  #x<-cuttil[cuttil$plot==i,]
  #sum_covar2[i,1]<-i
  #sum_covar2[i,2]<-sum(find_covar_sum(x))
  #sum_covar2[i,3]<-length(unique(x$Year))#same as nsamples
  #sum_covar2[i,4]<-length(x$Abundance)
  #sum_covar2[i,5]<-length(x$Abundance)/length(unique(x$Year))
#}

#Formula: by plot, sum(biomass of species)/sqrt(sum(var(biomass of species))+sum(cov(biomass of species)))
#Sensu Lehman and Tilman, 2000. American Naturalist, Vol. 156, No. 5 (November 2000), pp. 534-552
Lehman_Tilman_Stability<-sum_abund/sqrt(sum_var+sum_covar)

Lehman_Tilman_Stability
#Small problem. In some plots, there is literally no change. What to do in those cases? Answer is infinite, but that seems problematic for a mean...
#I feel like I can sub in with something, or I could omit, but not changing does mean something...

temp$tilman_stab<-Lehman_Tilman_Stability
temp$tilstab_corr<-ifelse(temp$tilman_stab==Inf, 50, temp$tilman_stab)
temp$Location<-substr(temp$Lsitesubplot,1,3)
meantilstab<-aggregate(cbind(temp$tilman_stab,temp$tilstab_corr), by=c(temp["Location"]), FUN=mean)
meantilstab$Location[meantilstab$Location=="CRM"]<-"CRMS"
meantilstab$Location[meantilstab$Location=="GCE"]<-"SAP"
colnames(meantilstab)<-c("Site","tilman_stab_mean","tilstab_corr_mean")
write.csv(meantilstab, "C:/Users/Akana/Documents/2012-2013/Environmental variability/Plant/tilstab.csv", row.names=FALSE)



# Proper Tilman stability at site level -----------------------------------
#If we do at site level, then we can include PDB and SO
#Best way to do? Means across all plots for each site and year and then calculate stability on those?
#First, mean abundance by species for each site and year
til<-read.csv("C:/Users/Akana/Documents/2012-2013/Environmental variability/Plant/all_plant_long.csv")
til.site<-summaryBy(Abundance~Species+unitAbund+Coast+Location+Year, data=til, FUN=mean)

temp<-as.data.frame(unique(til.site$Location))
temp$plot<-seq(1:11)
colnames(temp)<-c("Location","plot")
til.site<-merge(til.site,temp, by="Location")

sum_abund<-rowSums(tapply(til.site$Abundance.mean, list(til.site$plot, til.site$Species), mean))
#sum of variance for each species in each plot
sum_var<-rowSums(tapply(til.site$Abundance.mean, list(til.site$plot, til.site$Species), var))

#NAs are a problem for this code. Why do I have them?
find_covar_sum<-function(x) {
  nspecies<-length(unique(x$Species))
  nsamples<-length(unique(x$Year))
  
  comparisons<-outer(1:nspecies, 1:nspecies, ">")
  
  #Make a matrix with species biomass for the plot in each column, and rows ordered by sampling event
  sp_mat<-matrix(nrow=nsamples, ncol=nspecies, data=x$Abundance.mean[order(paste(x$Species, x$Year))])
  
  return(var(sp_mat)[comparisons]) #Record only unique, non-diagonal values 
}

#Run function for each plot
nsites<-length(unique(til.site$plot))
sum_covar<-numeric(nsites)
n<-1
for(i in sort(unique(til.site$plot))) {
  x<-til.site[til.site$plot==i,]
  sum_covar[n]<-sum(find_covar_sum(x))
  names(sum_covar)[n]<-i
  n<-n+1
}

#Error message here, I think just for plot 1, but I can't really figure out why and it seems to calculate okay

#nsites<-length(unique(cuttil$plot))
#sum_covar2<-data.frame(cbind(plot=numeric(0),sum_covar=numeric(0),nrow=numeric(0),datalength=numeric(0),ratio=numeric(0)))
#for(i in sort(unique(cuttil$plot))) {
#x<-cuttil[cuttil$plot==i,]
#sum_covar2[i,1]<-i
#sum_covar2[i,2]<-sum(find_covar_sum(x))
#sum_covar2[i,3]<-length(unique(x$Year))#same as nsamples
#sum_covar2[i,4]<-length(x$Abundance)
#sum_covar2[i,5]<-length(x$Abundance)/length(unique(x$Year))
#}

#Formula: by plot, sum(biomass of species)/sqrt(sum(var(biomass of species))+sum(cov(biomass of species)))
#Sensu Lehman and Tilman, 2000. American Naturalist, Vol. 156, No. 5 (November 2000), pp. 534-552
Lehman_Tilman_Stability<-sum_abund/sqrt(sum_var+sum_covar)

Lehman_Tilman_Stability

temp$site.tilman_stab<-Lehman_Tilman_Stability
temp$Location<-as.character(temp$Location)
temp$Location[temp$Location=="GCE"]<-c("SAP")
temp$Location<-as.factor(temp$Location)
tilstab<-temp[,c(1,3)]
colnames(tilstab)<-c("Site","site.tilstab")
write.csv(tilstab, "C:/Users/Akana/Documents/2012-2013/Environmental variability/Plant/site_tilstab.csv", row.names=FALSE)


# Proper Tilman stability and covariance ratio, site-level with moving five-year window --------

#If we do at site level, then we can include PDB and SO
#Best way to do? Means across all plots for each site and year and then calculate stability on those?
#First, mean abundance by species for each site and year
til<-read.csv("~/Documents/PhD docs/2012-2013/Environmental variability/Plant/all_plant_long.csv")
til.site<-summaryBy(Abundance~Species+unitAbund+Coast+Location+Year, data=til, FUN=mean)
library(data.table)

#Function we need to do Tilman stability
find_covar_sum<-function(x) {
  nspecies<-length(unique(x$Species))
  nsamples<-length(unique(x$Year))
  
  comparisons<-outer(1:nspecies, 1:nspecies, ">")
  
  #Make a matrix with species biomass for the plot in each column, and rows ordered by sampling event
  sp_mat<-matrix(nrow=nsamples, ncol=nspecies, data=x$Abundance.mean[order(paste(x$Species, x$Year))])
  
  return(var(sp_mat)[comparisons]) #Record only unique, non-diagonal values 
}

out<-data.frame(cbind(Location=character(0),tilstab=numeric(0),sum_abund=numeric(0),sum_var=numeric(0),sum_covar=numeric(0)))
sites<-unique(til.site$Location)
for(i in 1:length(sites)){
  sub<-subset(til.site,Location==sites[i])
  #Make sure ordered by year
  sub<-sub[order(sub$Year),]
  #Add index
  sub<-data.table(sub)
  sub[,index:=.GRP,by=Year]
  
  index<-unique(sub$index)
  for (j in 1:(length(index)-5)){
    subber<-subset(sub,j+5>=index&index>=j)
    
    #Calc Tilman stability
    sum_abund<-sum(tapply(subber$Abundance.mean, subber$Species, mean))
    sum_var<-sum(tapply(subber$Abundance.mean, subber$Species, var))
    sum_covar<-sum(find_covar_sum(subber))
    Lehman_Tilman_Stability<-sum_abund/sqrt(sum_var+sum_covar)
    
    temp<-as.data.frame(unique(subber$Location))
    temp$tilstab<-Lehman_Tilman_Stability
    temp$sum_abund<-sum_abund
    temp$sum_var<-sum_var
    temp$sum_covar<-sum_covar
    
    out<-rbind(out,temp)
  }
}

#Formula: by plot, sum(biomass of species)/sqrt(sum(var(biomass of species))+sum(cov(biomass of species)))
#Sensu Lehman and Tilman, 2000. American Naturalist, Vol. 156, No. 5 (November 2000), pp. 534-552

colnames(out)<-c("Site","tilstab","sum_abund","sum_var","sum_covar")
out$Site<-as.character(out$Site)
out$Site[out$Site=="GCE"]<-"SAP"
out$Site<-as.factor(out$Site)
negcov<-out[,c(1,3:5)]
write.csv(out, "~/Documents/PhD docs/2012-2013/Environmental variability/Plant/site_tilstab_moving.csv", row.names=FALSE)

#Covariance ratio from Hallett et al 2014 (but really from an earlier paper) is var(community)/sum(var(pop)) where var(comm)=sum(var(pop))+2*sum(covar)
#Can do for each five year window
negcov$var.ratio<-(negcov$sum_var+2*negcov$sum_covar)/negcov$sum_var
summary(negcov)
write.csv(negcov, "~/Documents/PhD docs/2012-2013/Environmental variability/Plant/site_covar_moving.csv", row.names=FALSE)

# Bray-curtis and Tilman look very different ------------------------------
#Let's see which makes more sense by doing nMDS of each year for plots to see if ones moving around more do look different
#Need data to be in vegan format, presumably
#Can't do nMDS for this giant number of plots, so average by year for each site
#Do a comparison where I choose just a few random plots from each site just to make sure things look the same-ish
#Then make subset for each site, run nMDS and plot on same plot
cutplants2<-cutplants[,-c(2:3,12)]
cutveg<-reshape(cutplants2, timevar="Species",idvar=c("Experiment","Plot","Subplot","Year","Location","unitAbund","Coast","freq"),direction="wide")
cutvegav<-aggregate(cutveg[,c(9:72)], by=list(Location=cutveg$Location, Year=cutveg$Year), FUN=mean, na.rm=T)

CARveg<-subset(cutvegav, Location=="CAR")
rownames(CARveg)<-CARveg$Year
CARveg<-CARveg[,-c(1:2)]
CRMSveg<-subset(cutvegav, Location=="CRMS")
rownames(CRMSveg)<-CRMSveg$Year
CRMSveg<-CRMSveg[,-c(1:2)]
GCEveg<-subset(cutvegav, Location=="GCE")
rownames(GCEveg)<-GCEveg$Year
GCEveg<-GCEveg[,-c(1:2)]
LPLveg<-subset(cutvegav, Location=="LPL")
rownames(LPLveg)<-LPLveg$Year
LPLveg<-LPLveg[,-c(1:2)]
NARveg<-subset(cutvegav, Location=="NAR")
rownames(NARveg)<-NARveg$Year
NARveg<-NARveg[,-c(1:2)]
PAPveg<-subset(cutvegav, Location=="PAP")
rownames(PAPveg)<-PAPveg$Year
PAPveg<-PAPveg[,-c(1:2)]
PIEveg<-subset(cutvegav, Location=="PIE")
rownames(PIEveg)<-PIEveg$Year
PIEveg<-PIEveg[,-c(1:2)]
SFBveg<-subset(cutvegav, Location=="SFB")
rownames(SFBveg)<-SFBveg$Year
SFBveg<-SFBveg[,-c(1:2)]
VCRveg<-subset(cutvegav, Location=="VCR")
rownames(VCRveg)<-VCRveg$Year
VCRveg<-VCRveg[,-c(1:2)]

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

CARveg[is.nan(CARveg)] <- 0
CRMSveg[is.nan(CRMSveg)] <- 0
GCEveg[is.nan(GCEveg)] <- 0
LPLveg[is.nan(LPLveg)] <- 0
NARveg[is.nan(NARveg)] <- 0
PAPveg[is.nan(PAPveg)] <- 0
PIEveg[is.nan(PIEveg)] <- 0
SFBveg[is.nan(SFBveg)] <- 0
VCRveg[is.nan(VCRveg)] <- 0

CARn<-metaMDS(CARveg,k=2,trymax=100)
stressplot(CARn)
ordiplot(CARn, type='n', ylim=c(-0.4,0.6),xlim=c(-1,1))
orditorp(CARn, display="sites", col="black")
LPLn<-metaMDS(LPLveg,k=2,trymax=100)
stressplot(LPLn)
ordiplot(LPLn, type='n')
orditorp(LPLn, display="sites", col="red")
SFBn<-metaMDS(SFBveg,k=2,trymax=100)
stressplot(SFBn)
ordiplot(SFBn, type='n')
orditorp(SFBn, display="sites", col="blue")
VCRn<-metaMDS(VCRveg,k=2,trymax=100)
stressplot(VCRn)
ordiplot(VCRn, type='n')
orditorp(VCRn, display="sites", col="green")
GCEn<-metaMDS(GCEveg,k=2,trymax=100)
stressplot(GCEn)
ordiplot(GCEn, type='n')
orditorp(GCEn, display="sites", col="pink")
CRMSn<-metaMDS(CRMSveg,k=2,trymax=100)
stressplot(CRMSn)
ordiplot(CRMSn, type='n')
orditorp(CRMSn, display="sites", col="black")
#This doesn't seem consistent with either Tilman stabilities or Bray-Curtis distances, but this probably doesn't really take into account the covariances that Tilman stability does. Let's try a different way to visualize


# Plot total abundances vs. species abundances over time ------------------
#For each site, I need a plot of total abundance and plots of each species' abundance over time to see if communities with greater Tilman stability have more covarying or less variance than others
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

cutvegav[is.nan(cutvegav)] <- 0
cutvegav$Total<-rowSums(cutvegav[,c(3:66)])

sites<-unique(cutvegav$Location)

for (i in 1:length(sites)){
  sub<-subset(cutvegav, cutvegav$Location==sites[i])
   
  plot(Total~Year, data=sub, ylim=c(0,max(sub$Total)+5), main=sites[i])
  lines(Total~Year, data=sub)
  
  color<-rainbow(63)
  
  for (j in 3:66) {
    if (sum(sub[,c(j)])>0){
      lines(sub[,c(j)]~sub$Year, col=color[j])
    }
  }
}

#This makes me feel pretty good about Tilman stabilities! Things seem pretty consistent with the possible exception of SFB. Might be worth looking at stability of total biomass, too

tilstab



# Variance of total abundance  -----------------------------------
tot<-read.csv("~/Documents/PhD docs/2012-2013/Environmental variability/Plant/raw_abund_4+_years.csv")
tot<-aggregate(tot$Abundance, by=c(tot["Lsitesubplot"],tot["Year"]), FUN=sum)
tot<-aggregate(tot$x, by=tot["Lsitesubplot"], FUN=var)
colnames(tot)<-c("Lsitesubplot","Variance")
tot$Lsitesubplot<-as.character(tot$Lsitesubplot)
tot$Site<-substr(tot$Lsitesubplot,1,3)
total.var<-aggregate(tot$Variance, by=tot["Site"],FUN=mean)
total.var$Site[total.var$Site=="CRM"]<-c("CRMS")
total.var$Site[total.var$Site=="GCE"]<-c("SAP")
write.csv(total.var,"C:/Users/Akana/Documents/2012-2013/Environmental variability/Plant/total_variance.csv", row.names=FALSE)



# Species turnover among consecutive years at site level ------------------
#Like Bray Curtis, at site-level so use all sites that have data in consecutive years even if not at same plot
#Also like Bray Curtis, need this year and prev year
#Since at site-level and just presence-absence at site, sum abund for each species at each site
plantless<-summaryBy(Abundance~Species+Location+Year,data=plant, FUN=sum)
plantless$yrminus<-plantless$Year + 1
yrminus<-plantless[, c("Location","yrminus","Species","Abundance.sum")]
year<-plantless[, c("Location","Year","Species","Abundance.sum")]

#rename variables 
yrminus2 <- rename(yrminus, c(yrminus="Year", Abundance.sum="abundprioryr"))

#merge & calculate abs diff between abund in year and prior yr
allyears<-merge(year, yrminus2, by=c("Location","Species", "Year"),all=TRUE)
#With merge (and all=TRUE), every species has a row if it or prior year abundance isn't zero. Replace those NAs with 0
allyears[is.na(allyears)]<-0

allyears$gained<-ifelse(allyears$Abundance.sum>0 & allyears$abundprioryr==0, 1, 0)
allyears$lost<-ifelse(allyears$Abundance.sum==0 & allyears$abundprioryr>0,1,0)
allyears$total<-ifelse(allyears$Abundance.sum>0 | allyears$abundprioryr>0,1,0)

#For each Year, sp turnover=(#sp gained+#lost)/total so first need to sum all of those counters by year and site
spturn<-summaryBy(gained+lost+total~Location+Year, data=allyears, FUN=sum,na.rm=TRUE)
spturn$turnover<-(spturn$gained.sum+spturn$lost.sum)/spturn$total.sum
turnover<-summaryBy(turnover~Location, data=spturn, FUN=c(mean,var))
turnover$Location<-as.character(turnover$Location)
turnover$Location[turnover$Location=="GCE"]<-"SAP"
turnover$Location<-as.factor(turnover$Location)

write.csv(turnover,"~/Documents/PhD docs/2012-2013/Environmental variability/Plant/turnover_sitelevel.csv",row.names=FALSE)


# Synchrony using Loreau and de Mazancourt method -------------------------
unbound<-subset(plant, bounded=="no")
library(data.table)
out<-data.frame(cbind(Location=character(0),synchrony=numeric(0)))
sites<-unique(unbound$Location)
for(i in 1:length(sites)){
  sub<-subset(unbound,Location==sites[i])
  #Make sure ordered by year
  sub<-sub[order(sub$Year),]
  #Add index
  sub<-data.table(sub)
  sub[,index:=.GRP,by=Year]
  
  index<-unique(sub$index)
  for (j in 1:(length(index)-5)){
    subber<-subset(sub,j+5>=index&index>=j)
    
    agg_abund<-aggregate(Abundance~Year, data=subber, FUN=sum)
    agg_abund$var<-var(agg_abund$Abundance)  
    a<-aggregate(Abundance~Species+Year, data=subber, FUN=sum)
    b<-aggregate(Abundance~Species, data=a, FUN=sd)
    agg_abund$sp<-sum(b$Abundance, na.rm=TRUE)
    
    temp<-as.data.frame(unique(subber$Location))
    temp$synchrony<-unique((agg_abund$var)/(agg_abund$sp^2))
    
    out<-rbind(out,temp)
  }
}
colnames(out)<-c("Site","synchrony")
out$Site<-as.character(out$Site)
out$Site[out$Site=="GCE"]<-"SAP"
out$Site<-as.factor(out$Site)
d<-aggregate(synchrony~Site, data=out, FUN=mean)

write.csv(d, "~/Documents/PhD docs/2012-2013/Environmental variability/Plant/synchrony_moving.csv", row.names=FALSE)
