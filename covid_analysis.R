library(dplyr)
library(ggmap)
library(ggplot2)
library(raster)
library(rgeos)
library(maptools)
library(rgdal)
library(inlmisc)
library(fda)
library(funFEM)
library(forecast)
library(funHDDC)
library(mclust)
library(cvTools)
library(clue)
source("functions.R")

###################################################
######### Read Data and Pre-processing ############ 
###################################################
covid=read.csv("data/covid19.csv")

head(covid)

covid=covid[,-1]
dim(covid)

date=covid[,1]
### convert yymmddhh to yymmdd
gr <- sort( as.POSIXct(unique(date), format = "%Y.%m.%d.") )
class(gr)
gr[1:5]

covid_data=  covid[,-c(1:4)]

head(covid_data)
covid_data=apply(covid_data, 2, rev)
head(covid_data)
###set the length of time-series to be 500
covid_data=covid_data[1:500,]
gr=gr[1:500]

dim(covid_data)

data=covid_data[,-26]
eng.name=c(  "Jongno-gu",  "Jung-gu" ,  "Yongsan-gu", " Seongdong-gu" ,  "Gwangjin-gu",
             " Dongdaemun-gu ", " Jungnang-gu ",
             " Seongbuk-gu" , " Gangbuk-gu ",  "Dobong-gu" ,  "Nowon-gu" ,  "Eunpyeong-gu",   "Seodaemun-gu",
             "Mapo-gu ","  Yangcheon-gu",  "Gangseo-gu",  "Guro-gu",  "Geumcheon-gu",
             "Yeongdeungpo-gu","  Dongjak-gu" ,  "Gwanak-gu ", " Seocho-gu" ,"Gangnam-gu", 
             " Songpa-gu" , " Gangdong-gu")


length(which(data==0) )/length(data)


for(j in c(1:ncol(data))){
  print( c(  eng.name[j], sum(data[,j]) ,  length(which(data[,j]==0) ) ,
             round( length(which(data[,j]==0) ) /500  ,3 ), max(data[,j])   ))
}

### EDA plot of Number of new COVID-19; total
par(mfrow=c(1,1))
plot(apply(data,1,sum), col=1, type='b' , ylim=c(0,600),  xaxt='n', xlab='', ylab='Number of new COVID-19')
axis(side=1, at=seq(1, length(gr), length.out=7), labels=gr[seq(1, length(gr), length.out=7)] , cex.axis=1) 


### EDA plot of Number of new COVID-19 from 3 different districts
par(mfrow=c(1,3))
for(j in c(1,16,2)){
  plot(data[,1], col=1, type='n' , ylim=c(0,80),  xaxt='n', xlab='', ylab='Number of new COVID-19', main=eng.name[j])
  axis(side=1, at=seq(1, length(gr), length.out=7), labels=gr[seq(1, length(gr), length.out=7)] , cex.axis=1) 
  lines( data[,j] , col=1, type='b')
  print(length(which(data[,j]==0) )/length(data[,j]))
}
dev.off()

### EDA histogram of COVID-19 New cases
hist(data, main='Histogram of COVID-19 New cases', xlab='')



###################################################
############### Number of clusters ################# 
###################################################
K=3

t1=5
r1=0.1
Up_covid = prepare_up(data,FUN=sqBound_es,t1,r1) 
hclusCut <- function(data, k) list(cluster = cutree(hclust(dist(data,method = "manhattan" )), k=k))
gap.covid <- clusGap(t(log(Up_covid)), FUN = hclusCut, K.max = 10, B = 100)
plot(gap.covid)

################################################################
########################## Do Clustering  ######################
################################################################

## various t value ####
t1_index=c(10,30,100)
yhat1=matrix(nrow=ncol(data), ncol=length(t1_index))

for(k in 1:length(t1_index)){
  t1=t1_index[k]
  r1=0.1
  Up_covid = prepare_up(data,FUN=sqBound_es,t1,r1)
  
  a1= cluster_zits(K=K,Up=log(Up_covid), FUN=l1_mean ,t1,r1,is.median=1)
  yhat=vector()
  for(l in 1:K ){
    yhat[a1$clusters[[l]]]=l
  }
  
  yhat1[,k]=yhat
  
}


save(yhat1, file='./R_result/covid_yhat.RData')
#load('./R_result/covid_yhat.RData')

################################################################
################## Comparision Clustering  ######################
################################################################


basis <- create.fourier.basis(c(0, nrow(data)), nbasis=5)
fdobj <- smooth.basis(1:nrow(data), data,basis)$fd
cl.sinu = funFEM(fdobj, K=K, init="kmeans")
elements = list()
for(c in 1:K){
  # get the point of cluster c
  elements[[c]] <- which(cl.sinu$cls == c)
}
yhat_funfem=vector()
for(l in 1:K ){
  yhat_funfem[elements[[l]]]=l
}

res.uni <- funHDDC(fdobj,K=K,model="AkBkQkDk",init="kmeans",threshold=0.2)
yhat_funhddc= res.uni$class


### DTW cluster
pc.l2 <-tsclust(t(data), k = K,  distance = "dtw", centroid = "pam",seed = 3247, trace = FALSE,  control = partitional_control(nrep = 5L))
yhat_dtw=cl_class_ids(cl_medoid(pc.l2))



#############################################
############### Result  ###############
#############################################
for(k in 1:length(t1_index)){
  print(table(yhat1[,k]))
}


cbind(eng.name,yhat1[,3] )



j=2
x=yhat1[,j]
#table(yhat_funfem , x)
x2=ifelse(x == 1, 2,  (ifelse(x == 2,1,3)  ) )  ## change 1 & 2
x2=ifelse(x == 1, 3,  (ifelse(x == 3,1,2)  ) )   ## change 1 &3    
x2=ifelse(x == 2, 3,  (ifelse(x == 3,2,1)  ) )   ## change 2 &3  
#table(yhat_funfem , x3)
yhat1[,j]=x2


par(mfrow=c(3,K))
for(t1 in c(1:3)){
  for(k in 1:K){
    
    plot(  apply( data[, which(yhat1[,k]==2 )] ,1, mean), col=1, type='n' ,main=paste('Thickness=',t1_index[t1]), xaxt='n', ylim=range(data), ylab='', xlab='')
    axis(side=1, at=seq(1, length(gr), length.out=7), labels=gr[seq(1, length(gr), length.out=7)] , cex.axis=1) 
    
    temp=data[, which(yhat1[,t1]==k )] 
    if(length(temp)==500){
      lines( temp, col=GetTolColors(5, "muted",  gray = TRUE)[c(2,3,5)][k],type='l') 
    }
    if(length(temp)>500){
      for(j in 1:ncol(temp)){
        lines( temp[,j], col=GetTolColors(5, "muted",  gray = TRUE)[c(2,3,5)][k],type='l') 
      }
    }
  }
  
}


j=3

tb2=matrix(nrow=4, ncol=K)
for(k in 1:K){
  
  if(table(yhat1[,j])[k]>1 ){
    summ= c(  table(yhat1[,j])[k],
              mean(data[, which(yhat1[,j]==k)]),
              #range(  apply(data[, which(yhat1[,j]==k)],2, mean) ),
              sqrt( var(  apply(data[, which(yhat1[,j]==k)],2, mean) )),
              length(which( data[, which(yhat1[,j]==k)]==0)) / (length(which(yhat1[,j]==k))*500 )*100
    )
    tb2[,k]=(round(summ,3))
  }
  if(table(yhat1[,j])[k]==1 ){
    summ= c(  table(yhat1[,j])[k],
              mean(data[, which(yhat1[,j]==k)]),
              #range(  apply(data[, which(yhat1[,j]==k)],2, mean) ),
              NA,
              length(which( data[, which(yhat1[,j]==k)]==0)) / (length(which(yhat1[,j]==k))*500 )*100
    )
    tb2[,k]=(round(summ,3))
  }
  
}
tb2

###############  cluster validation  ############### 
dim(data)
clust_val<-matrix(nrow=6, ncol=4)
df=t(data)
dim(df)
library(fpc)

for(j in 1:3){
  clust_stats <- cluster.stats(dist(df),  yhat1[,j])
  clust_val[j,]=round(c( clust_stats $dunn, clust_stats$avg.silwidth, clust_stats $entropy , clust_stats $ch),3)
}


clust_stats <- cluster.stats(dist(df),  yhat_funfem)
clust_val[4,]=round(c( clust_stats $dunn, clust_stats$avg.silwidth, clust_stats $entropy , clust_stats $ch ),3)


clust_stats <- cluster.stats(dist(df),  yhat_funhddc)
clust_val[5,]=round(c( clust_stats $dunn, clust_stats$avg.silwidth, clust_stats $entropy , clust_stats $ch ),3)

clust_stats <- cluster.stats(dist(df),  yhat_dtw)
clust_val[6,]=round(c( clust_stats $dunn, clust_stats$avg.silwidth, clust_stats $entropy , clust_stats $ch ),3)


colnames(clust_val)=c('Dunn','silwidth','entropy','CH')
rownames(clust_val)=c('tau10','tau30','tau100','FEM','HDDC', 'DTW')
clust_val


## plot of comparision methods
x=yhat_dtw
table(yhat_funfem , x)
x2=ifelse(x == 1, 2,  (ifelse(x == 2,1,3)  ) )  ## change 1 & 2
x2=ifelse(x == 1, 3,  (ifelse(x == 3,1,2)  ) )   ## change 1 &3    
x3=ifelse(x2 == 2, 3,  (ifelse(x2 == 3,2,1)  ) )   ## change 2 &3  
table(yhat_funfem , x2)
yhat_dtw=x2

par(mfrow=c(3,K))
for(k in 1:K){
  plot(  apply( data[, which(yhat_funfem==1 )] ,1, mean), col=1, type='n' ,main='FunFEM', xaxt='n', ylim=range(data), ylab='', xlab='')
  axis(side=1, at=seq(1, length(gr), length.out=7), labels=gr[seq(1, length(gr), length.out=7)] , cex.axis=1) 
  
  temp=data[, which(yhat_funfem==k )] 
  if(length(temp)==500){
    lines( temp, col=GetTolColors(5, "muted",  gray = TRUE)[c(2,3,5)][k]) 
  }
  if(length(temp)>500){
    for(j in 1:ncol(temp)){
      lines( temp[,j], col=GetTolColors(5, "muted",  gray = TRUE)[c(2,3,5)][k])
    }
  }
}
for(k in 1:K){
  plot(  apply( data[, which(yhat_funhddc==1 )] ,1, mean), col=1, type='n' ,main='FunHDDC', xaxt='n', ylim=range(data), ylab='', xlab='')
  axis(side=1, at=seq(1, length(gr), length.out=7), labels=gr[seq(1, length(gr), length.out=7)] , cex.axis=1) 
  
  temp=data[, which(yhat_funhddc==k )] 
  if(length(temp)==500){
    lines( temp, col=GetTolColors(5, "muted",  gray = TRUE)[c(2,3,5)][k]) 
  }
  if(length(temp)>500){
    for(j in 1:ncol(temp)){
      lines( temp[,j], col=GetTolColors(5, "muted",  gray = TRUE)[c(2,3,5)][k])
    }
  }
}
for(k in 1:K){
  plot(  apply( data[, which(yhat_dtw==1 )] ,1, mean), col=1, type='n' ,main='DTW', xaxt='n', ylim=range(data), ylab='', xlab='')
  axis(side=1, at=seq(1, length(gr), length.out=7), labels=gr[seq(1, length(gr), length.out=7)] , cex.axis=1) 
  
  temp=data[, which(yhat_dtw==k )] 
  if(length(temp)==500){
    lines( temp, col=GetTolColors(5, "muted",  gray = TRUE)[c(2,3,5)][k]) 
  }
  if(length(temp)>500){
    for(j in 1:ncol(temp)){
      lines( temp[,j], col=GetTolColors(5, "muted",  gray = TRUE)[c(2,3,5)][k])
    }
  }
}


# plot 
load('data/P_merge.RData')


gu=unique(P_merge$id)
num_gu<-vector()
for(j in 1:length(gu)){
  num_gu[j]=length( which(P_merge$id==gu[j]  ) )
}


pal=  GetTolColors(5, "muted",  gray = TRUE)[c(2,3,5)][rep(yhat1[,3], times=num_gu )]

ggplot() + geom_polygon(data = P_merge, aes(x=long, y=lat, group=group), fill = pal, color='black')+
  ylab("Latitude")+ggtitle("Thickness 100")+
  xlab("Longitude")+
  theme_bw() +
  scale_shape_manual(values=c(1, 8))+
  theme(legend.text=element_text(size=11,face = "bold"),
        axis.title = element_text(size=13,face = "bold"),
        axis.text=element_text(size=11,face = "bold"),
        legend.title =element_text(size=13,face = "bold"),
        #legend.position = c(0.25, 0.9)
        legend.position = "bottom", plot.title = element_text( size=14, face="bold"),
  )
#dev.off()
