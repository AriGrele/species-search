library(rgbif);library(stringr);library(ggplot2);library(ggforce)

clamp=function(x,l,u){return(ifelse(x<l,l,ifelse(x>u,u,x)))}

community=function(pos,radius,taxon,limit=300){                                 #search for communities of specified taxon within circular region
  circle=seq(0,2*3.14159,length.out=10)                                         #10 angles for a circle            
  x=clamp(cos(circle)*radius+pos[1],-180,180)                                   #x and y coord on circle
  x=c(x,x[1])                                           
  y=clamp(sin(circle)*radius+pos[2],-90,90)
  y=c(y,y[1])
  
  data=occ_search(limit=30000,                                                  # search for most recent n taxa within region of interest. Every 300 results requires an additional search, 100,000 max limit
                  datasetKey='50c9509d-22c7-4a22-a47d-8c48425ef4a7',            #Key for Inaturalist searches
                  taxonKey=taxon,                                                   
                  geometry=paste0('POLYGON((',paste(paste(x,y),collapse=','),'))'))
  return(data)}

commerge=function(data,output,i=NA){                                            #get colnames to match up with output by removing columns unique to only one dataframe
  d=data$data[,which(names(data$data)%in%names(output))]                        
  d$search=i
  o=output[which(names(output)%in%names(d))]
  if(!is.null(data$data)){return(rbind(o,d))}                                   #bind with output if any results found (should always be >=1 result, but just in case)}
  return(out)}

#############################data scraping#######################################

radius=1/111                                                                    #radius (degrees) for search area
start=community(c(0,0),180,'4987990',1000)                                      #collect all Merope tuber observations
lat =start$data$decimalLatitude                                                 #lat/long of observations
long=start$data$decimalLongitude
out=as.data.frame(matrix(ncol=ncol(start$data)+1))|>                            #output dataframe
    setNames(c('search',names(start$data)))

for(i in 1:length(lat)){                                                        #for each location:
  data=community(c(long[i],lat[i]),radius,6,30000)                              #scrape data for certain radius
  out=commerge(data,out,i)                                                      #merge with output
  cat(i,data$meta$count,'\n')
  Sys.sleep(2)}                                                                 #sleep 2 seconds to make sure rate-limit doesn't get triggered

pos=c(-76.4313087,42.4637713)  
data=community(pos,radius,6,30000)
out=commerge(data,out,'test7')

out$og_lat=sapply(out$search,\(x)lat[as.numeric(x)])                            #lat/long of original observation
out$og_long=sapply(out$search,\(x)long[as.numeric(x)])
out$distance=sqrt((out$decimalLongitude-out$og_long)^2+                         #distance between test observation and other taxon
             (out$decimalLatitude-out$og_lat)^2)

#################################stats###########################################
taxon='species'
o=tapply(out$genus,list(out$search,out[,taxon]),length)                         #table of taxa abundances for each location
for(i in 1:ncol(o)){o[,i][is.na(o[,i])]=0}                                      #NA => 0

pca=prcomp(o)
cluster=(nrow(pca$x)-1)*sum(apply(pca$x,2,var))
for(i in 2:15){cluster[i]=sum(kmeans(pca$x,centers=i,nstart=25,iter.max=1000)$withinss)}
plot(1:15,cluster,type="b",xlab="Number of Clusters",
     ylab="Within groups sum of squares")
k=kmeans(pca$x,6,nstart=25,iter.max=1000)

#clust=names(sort(table(k$clust)))
#View(pca$rotation[k$clust==clust[1],])

first=as.data.frame(pca$x[,1:2])
r=first[nrow(first),]
for(i in 1:nrow(first)){
  first$dist[i]=sqrt((r[,1]-first[i,1])^2+(r[,2]-first[i,2])^2)}
first$cluster=k$cluster

fa=with(first,aggregate(cluster,list(cluster),length))
large=fa$Group.1[fa$x==max(fa$x)]
hab_percent=max(fa$x)/sum(fa$x)

drange=max(first$dist/20)
first=first[order(first$dist),]
first=first[first$dist<=drange,]

loc_percent=sum(first$cluster==large)/nrow(first)

hab_percent;loc_percent

data=as.data.frame(pca$x)
data$col=factor(k$clust,levels=fa$Group.1[order(fa$x)])
ggplot(data)+
  geom_point(aes(x=PC1,y=PC2,color=col),size=3)+
  theme_classic()+
  scale_color_manual(values=pals::ocean.haline(7)[-7])+
  geom_circle(aes(x0=first[1,1],y0=first[1,2],r=drange),color='black')+
  geom_point(aes(x=first[1,1],y=first[1,2]),size=3)

a=aggregate(out$search,list(out[,taxon]),length)
a$nsite=aggregate(out$search,list(out[,taxon]),\(x)length(unique(x)))$x
a=a[order(-a$x),]

out=out[!is.na(out[,taxon]),]

means=as.data.frame(matrix(ncol=6))|>
      setNames(c('taxon','distance','meters','stdev','count','nsite'))
n=10
for(i in 1:nrow(a[a$nsite>n,])){
  group=out$distance[out[,taxon]==a$Group.1[a$nsite>n][i]]
  dist=mean(group,na.rm=T)
  means=rbind(means,data.frame('taxon'=a$Group.1[a$nsite>n][i],
                               'distance'=dist,
                               'meters'=dist*111*1000,
                               'stdev'=sd(group*111*1000,na.rm=T),
                               'count'=a$x[a$nsite>n][i],
                               'nsite'=a$nsite[a$nsite>n][i]))}
head(means[order(-means$nsite),],20)



