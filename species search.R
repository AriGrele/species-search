library(rgbif);library(stringr);library(ggplot2)

key='50c9509d-22c7-4a22-a47d-8c48425ef4a7'                                      #Key for Inaturalist searches
r=0.25                                                                          #radius (degrees) for search area

start=occ_search(limit=5000,                                                    #collect all Merope tuber observations
           datasetKey=key,
           taxonKey='4987990')

lat=start$data$decimalLatitude                                                  #lat/long of observations
long=start$data$decimalLongitude

out=as.data.frame(matrix(ncol=ncol(start$data)+1))|>                            #output dataframe
    setNames(c('search',names(start$data)))
plot(long,lat,type='n')                                                         #plot locations
for(i in 1:length(lat)){                                                        #for each location
  print(i)
  pos=c(long[i],lat[i])                                                         # ordered pair of long/lat 
  circle=seq(0,2*3.14159,length.out=10)                                         # 10 angles for a circle
  x=cos(circle)*r+pos[1];x=c(x,x[1])                                            # x and y coord on circle
  y=sin(circle)*r+pos[2];y=c(y,y[1])
  lines(x,y)                                                                    #plot regions of interest
  
  data=occ_search(limit=300,                                                    #search for most recent n taxa within region of interest. Every 300 results requires additional search, 100,000 max limit
                  datasetKey=key,
                  geometry=paste0('POLYGON((',paste(paste(x,y),collapse=','),'))'))
  
  d=data$data[,which(names(data$data)%in%names(out))]                           #get colnames to match up with output
  d$search=i
  o=out[which(names(out)%in%names(d))]
  
  if(!is.null(data$data)){out=rbind(o,d)}                                       #bind with output if any results found (should always be >=1 result, but just in case)
  Sys.sleep(2)}                                                                 #sleep 2 seconds to make sure rate-limit doesn't get triggered

out$og_lat=sapply(out$search,\(x)lat[x])                                        #lat/long of original observation
out$og_long=sapply(out$search,\(x)long[x])
out$distance=sqrt((out$decimalLongitude-out$og_long)^2+                         #distance between observation and other taxon
             (out$decimalLatitude-out$og_lat)^2)

o=tapply(out$species,list(out$search,out$species),length)                       #table of taxa abundances for each location
for(i in 1:ncol(o)){o[,i][is.na(o[,i])]=0}                                      #NA => 0

pca=prcomp(o)
#View(pca$x)
ggplot(pca$x,aes(x=PC1,y=PC2,color=PC3))+
  geom_point()+
  theme_classic()


a=aggregate(out$species,list(out$species),length)
a=a[order(-a$x),]
means=data.frame('taxon'=NA,'distance'=NA,'meters'=NA,'stdev'=NA,'count'=NA)
for(i in 1:nrow(a[a$x>10,])){
  group=out$distance[out$species==a$Group.1[a$x>10][i]]
  dist=mean(group,na.rm=T)
  means=rbind(means,data.frame('taxon'=a$Group.1[a$x>10][i],
                               'distance'=dist,
                               'meters'=dist*111*1000,
                               'stdev'=sd(group*111*1000,na.rm=T),
                               'count'=a$x[a$x>10][i]))}
head(means[order(means$distance),])
