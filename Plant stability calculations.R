library(vegan)
library(ggplot2)
library(plyr)
library(doBy)
library(cgwtools)#for seqle

#See attached csv for what's in data file. I can't load the whole thing cause it has data from other people (but see paper for references to where data are from)
plant<-read.csv("~/Documents/PhD docs/2012-2013/Environmental variability/Plant/all_plant.csv")


# Rarefied species richness taking into account plots sampled and area ------------------
#Don't need consecutive years, so using original plant data
size<-read.table("~/Documents/PhD Docs/2012-2013/Environmental variability/Plot size.csv", header=T, sep=",")
plant2<-merge(plant,size, by="Location")

#I want to calculate mean species richness for each location, but I need to account for plot size and number of plots
#Species accumulation accounts for number of plots
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

#Now to make area sampled, I'm going to multiply samplesize by plot area to make x-axis area sampled
out2<-merge(output,size,by="Location")
out2$SampArea<-out2$samplesize*out2$Plot.size

#Now find smallest area that works for all
a<-summaryBy(SampArea~Location, data=out2, FUN=max)
#Smallest area is PIE at 55m2
calib<-subset(out2, SampArea<56)
#Look at curve to make sure sites are all leveled out
qplot(x=calib$SampArea, y=calib$richness, 
      data=calib, 
      colour=Location,
      main="GGPLOT line plot with groups") +
  geom_line()
calib2<-summaryBy(richness~Location, data=calib, FUN=max)
write.csv(calib2, "~/Documents/PhD Docs/2012-2013/Environmental variability/Plant/rare_sprich_byarea.csv", row.names=FALSE)


#Caluclating covariance ratio (measure of negative covariance)---------------------
#First, mean abundance by species for each site and year
til<-read.csv("~/Documents/PhD docs/2012-2013/Environmental variability/Plant/all_plant_long.csv")
til.site<-summaryBy(Abundance~Species+unitAbund+Coast+Location+Year, data=til, FUN=mean)
library(data.table)

#Function we need to find sum of covariances from Adam Clark (adamtclark on Github)
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
    
    sum_abund<-sum(tapply(subber$Abundance.mean, subber$Species, mean))
    sum_var<-sum(tapply(subber$Abundance.mean, subber$Species, var))
    sum_covar<-sum(find_covar_sum(subber))
    
    temp<-as.data.frame(unique(subber$Location))
    temp$sum_abund<-sum_abund
    temp$sum_var<-sum_var
    temp$sum_covar<-sum_covar
    
    out<-rbind(out,temp)
  }
}

colnames(out)<-c("Site","sum_abund","sum_var","sum_covar")
out$Site<-as.character(out$Site)
out$Site[out$Site=="GCE"]<-"SAP"
out$Site<-as.factor(out$Site)
negcov<-out

#Use those to calculate covariance ratio from Hallett et al 2014 (but really from an earlier paper) is var(community)/sum(var(pop)) where var(comm)=sum(var(pop))+2*sum(covar)
#Can do for each five year window
negcov$var.ratio<-(negcov$sum_var+2*negcov$sum_covar)/negcov$sum_var
summary(negcov)
write.csv(negcov, "~/Documents/PhD docs/2012-2013/Environmental variability/Plant/site_covar_moving.csv", row.names=FALSE)


# Species turnover among consecutive years at site level ------------------
#Do this at site-level so use all sites that have data in consecutive years even if not at same plot
#Since at site-level and just presence-absence at site, sum abund for each species at each site
#Need this year and prev year, so can create dummy variable
plantless<-summaryBy(Abundance~Species+Location+Year,data=plant, FUN=sum)
plantless$yrminus<-plantless$Year + 1
yrminus<-plantless[, c("Location","yrminus","Species","Abundance.sum")]
year<-plantless[, c("Location","Year","Species","Abundance.sum")]
#rename variables 
yrminus2 <- rename(yrminus, c(yrminus="Year", Abundance.sum="abundprioryr"))

#merge & calculate abs diff between abund in year and prior yr
allyears<-merge(year, yrminus2, by=c("Location","Species", "Year"),all=TRUE)
#Every species has a row if it or prior year abundance isn't zero. Replace those NAs with 0
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


#Aggregate species abundances (per Lehman and Tilman 2000) --------
#Here also only for plots with at least 5 consecutive years of data 
#Also use a moving window to make sure sites with more years of data aren't totally different just because of that
#This method only works when species abundances don't sum to 100, so remove sites that were limited to 100% cover
unbound<-subset(plant, bounded=="no")

#Aggregate species abundances within replicate and year, then calculate mean and sd for each replicate to get stability (mean/sd from Lehman and Tilman 2000)
library(data.table)
out<-data.frame(cbind(Location=character(0),aggstab=numeric(0)))
sites<-unique(unbound$Location)
for(i in 1:length(sites)){
  sub<-subset(unbound,Location==sites[i])
  #Order by year
  sub<-sub[order(sub$Year),]
  #Add index
  sub<-data.table(sub)
  sub[,index:=.GRP,by=Year]
  
  index<-unique(sub$index)
  #For five year periods
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
