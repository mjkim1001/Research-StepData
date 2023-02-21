source("./step_readData.R")

dim(step) ## input: 1440 * 21394

#step=t(step)

### 79 individual
length(step.list)

num.day<-vector()
for(i in 1: length(step.list) ){
    num.day[i]=ncol(step.list[[i]])
}
num.day
which.max(num.day)
range(num.day)


##thick pen 
K=6
t1_index=c(10,20,30,50,100)
result_step<-matrix(nrow=length(t1_index), ncol=ncol(step))

for(k in 1:length(t1_index)){
  t1=t1_index[k]
  r1=0.2
  Up_step = prepare_up(step,FUN=sqBound_es,t1,r1) ## Ensemble pen result
  
  a1= cluster_zits(K=K,Up=log(Up_step), FUN=l1_mean ,t1,r1,is.median=T)
  yhat1=vector()
  for(j in 1:K){
    yhat1[a1$e[[j]]]=j
  }
  result_step[k,]=yhat1

}

table( result_step[k,])
dim(result_step)


save(result_step, file='step_result.RData')
########################################################
############## cv code for thickness ####################
########################################################
data=step

t1_index=c(10,20,30,50,100)
cv_error<-vector()
for(k in 1:length(t1_index)){
  
  t1=t1_index[k]
  r1=0.2
  
  K_cv=5
  folds <- cvFolds(ncol(data), K=K_cv)
  
  temp_error_si<-vector()
  for(i in 1:K_cv){
    train <- data[,folds$subsets[folds$which != i] ] #Set the training set
    validation <- data[,folds$subsets[folds$which == i] ] #Set the validation set
    Up_step = prepare_up(train,FUN=sqBound_es,t1,r1) ## Ensemble pen result
    
    a1= cluster_zits(K=K,Up=log(Up_step), FUN=l1_mean ,t1,r1,is.median=T)
    for(l in 1:K){
      a1$pivot[,l]= apply( log(Up_step)[ ,a1$e[[l]]],1, median)}
    
    Up_test = prepare_up(validation,FUN=sqBound_es,t1,r1) ## Ensemble pen result
    
    pmeasure=matrix(nrow=ncol(Up_test) , ncol=K)
    for(clstr in 1:K){
      for(j in 1:ncol(Up_test)){
        pmeasure[j,clstr] = l1_mean(a1$pivot[,clstr],log(Up_test[,j]))
      }
    }
    
    yhat=apply(pmeasure,1, which.min) 
    #assign each validation set to the pivots derived using train set
    
    si2 <- silhouette((yhat), dist(t(validation), "euclidean"))
    #for each fold, the silhouette values are saved in temp_error_si 
    temp_error_si[i]= summary(si2)$avg.width
    
  }
  
  cv_error[k]=mean(temp_error_si)
  
}


opt_k=which.max(cv_error)
opt_k 
t1=t1_index[opt_k] #optimal thickness
Up_step = prepare_up(data,FUN=sqBound_es,t1,r1) ## Ensemble pen result

a1= cluster_zits(K=K,Up=log(Up_step), FUN=l1_mean ,t1,r1,is.median=T)
yhat=vector()
for(l in 1:K ){
  yhat[a1$e[[l]]]=l
}


cv_result=yhat

########################  comparison  ######################## 
################## Functional Clustering  ######################
################################################################

basis <- create.fourier.basis(c(0, nrow(data)), nbasis=11)
fdobj <- smooth.basis(1:nrow(data), data,basis)$fd
cl.sinu = funFEM(fdobj, K=K, init="kmeans", model="all")
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
  
yhat_hddc= res.uni$class

  
### DTW cluster
sam_id=sample(ncol(data), 150)
pc.l2 <-tsclust(t(data[,sam_id]), k = K,  distance = "dtw_basic", centroid = "pam",seed = 1, trace = FALSE,  control = partitional_control(nrep = 1L))
yhat_dtw=cl_class_ids((pc.l2))

###############  Plots for comparison methods  ############### 
tp = seq(as.POSIXct("00:00:00",format="%H:%M"),
           as.POSIXct("23:59:00",format="%H:%M"), by = "1 min")

par(mfrow=c(1,2))
for(k in c(2,5)){
  
  aa=(table(result_step[k,]))

  
  plot( tp, apply(step[,which(result_step[k,]==1)],1, mean), type='n', ylab='', xlab='time', col=1, 
        main=paste('Thickness - ',t1_index[k] ),  lwd=2, ylim=c(0,23)) 
#axis(side=1, at=seq(1, length(tp), length.out=7), labels=tp[seq(1, length(tp), length.out=7)])
  for(j in 1:6){ 
    lines( tp, apply(step[,which(result_step[k,]==order(aa, decreasing=T)[j])],1, mean), col=j,type='l',lwd=2 )
  }
  
  leg_tex<-vector()
  for(j in order(aa, decreasing=T)){
    leg_tex=c(leg_tex,aa[j]) 
    }
  
  legend('topright', paste('Group', c(1:length(aa)) ) , col=c(1:6), lty=1, lwd=2)
  
}

save(yhat_funfem, yhat_hddc,file='step_comparison_result.RData' )

png(file="./figure/step_com2.png",
width=1400, height=600)
par(mfrow=c(1,2))
aa=(table(yhat_funfem))
plot( tp, apply(step[,which(yhat_funfem ==1)],1, mean), type='n', ylab='', xlab='time', col=1, 
        main='FunFEM',  lwd=2, ylim=c(0,30)) 
#axis(side=1, at=seq(1, length(tp), length.out=7), labels=tp[seq(1, length(tp), length.out=7)])
for(j in 1:6){ 
  lines( tp, apply(step[,which(yhat_funfem ==order(aa, decreasing=T)[j])],1, mean), col=j,type='l',lwd=2 )
}
  
leg_tex<-vector()
for(j in order(aa, decreasing=T)){
  leg_tex=c(leg_tex,aa[j]) 
}
legend('topright', paste('Group', c(1:length(aa)) ) , col=c(1:6), lty=1, lwd=2)
 
aa=(table(yhat_hddc))
plot( tp, apply(step[,which(yhat_hddc ==1)],1, mean), type='n', ylab='', xlab='time', col=1, 
        main='FunHDDC',  lwd=2, ylim=c(0,30)) 
#axis(side=1, at=seq(1, length(tp), length.out=7), labels=tp[seq(1, length(tp), length.out=7)])
for(j in 1:6){ 
  lines( tp, apply(step[,which(yhat_hddc ==order(aa, decreasing=T)[j])],1, mean), col=j,type='l',lwd=2 )
}
  
leg_tex<-vector()
for(j in order(aa, decreasing=T)){
  leg_tex=c(leg_tex,aa[j]) 
}
legend('topright', paste('Group', c(1:length(aa)) ) , col=c(1:6), lty=1, lwd=2)

aa=(table(yhat_dtw))
plot( tp, apply(step[,sam_id][,which(yhat_dtw ==1)],1, mean), type='n', ylab='', xlab='time', col=1, 
        main='DTW',  lwd=2, ylim=c(0,30)) 
#axis(side=1, at=seq(1, length(tp), length.out=7), labels=tp[seq(1, length(tp), length.out=7)])
for(j in order(aa, decreasing=T)){ 
  lines( tp, apply(step[,sam_id][,which(yhat_dtw ==j)],1, mean), col=j,type='l',lwd=2 )
}
  
leg_tex<-vector()
for(j in order(aa, decreasing=T)){
  leg_tex=c(leg_tex,aa[j]) 
}
  
legend('topright', paste('Group', c(1:length(aa)) ) , col=order(aa, decreasing=T), lty=1, lwd=2)
  

dev.off()

############### Table 4 ############### 

weekend.percent = matrix(NA,nrow=length(t1_index),ncol=6)
for(k in c(2,5)){
  numstep=c()
  aa=(table(result_step[k,]))
  
  for(j in order(aa, decreasing=T)){
   numstep = c(numstep, mean(apply(step[,which(result_step[k,]==j)],2,sum)  )  )
   weekend.percent[k,j]=(round( length(which(is.wkend [which(result_step[k,]==j) ] ==T) )/   length(which(result_step[k,]==j))*100,1))

  }
  print( (table(result_step[k,]))[order(aa, decreasing=T) ] )
 print(round( numstep),2)
 print(weekend.percent[k,order(aa, decreasing=T) ])
}


##########   Only one person  ########## ########## ########## 
library(clusterSim)

dim(step.list[[67]])


for(tau in c(20,100)){
  Up_step = prepare_up(step.list[[67]],FUN=sqBound_es,t1,r1) ## Ensemble pen result
  
  fun.step <- function(x,k, t1=tau ){
    print(k)
    a = cluster_zits(K=k,Up=x, FUN=l1_mean ,t1,r1=0.2,is.median=T)
    clusters = rep(NA, ncol(step.list[[67]]))
    for(c in 1:k){
      clusters[(a$ls[[ which.min(a$obj) ]]$e)[[c]]] = c
    }
    return(list(cluster = clusters))
  }
  
  
  gap.step= clusGap(Up_step,fun.step,10)
  plot(gap.step)
}


###############  Plot  ############### 
ID=67
order(table(people.save))

table(is.wkend[which(people.save==ID)])

tp = seq(as.POSIXct("00:00:00",format="%H:%M"),
         as.POSIXct("23:59:00",format="%H:%M"), by = "1 min")
par(mfrow=c(2,2))
plot( tp, apply(step.list[[ID]][,which(is.wkend[which(people.save==ID)]==F)  ],1, mean), type='l', ylab='', xlab='time', col=1, 
          main='Weekday average',  lwd=2, ylim=c(0,25)) 
plot( tp, apply(step.list[[ID]][,which(is.wkend[which(people.save==ID)]==T)  ],1, mean), type='l', ylab='', xlab='time', col=1, 
      main='Weekend average',  lwd=2, ylim=c(0,25)) 


order( apply(step.list[[ID]],2, sum) )
which.max( apply(step.list[[ID]],2, sum) )

#par(mfrow=c(1,2))
plot( tp, step.list[[ID]][,order( apply(step.list[[ID]],2, sum) )[2] ], type='l', ylab='', xlab='time', col=1, 
      main='Minimum activity day',  lwd=2 , ylim=c(0,120)) 
plot( tp, step.list[[ID]][,which.max( apply(step.list[[ID]],2, sum) ) ], type='l', ylab='', xlab='time', col=1, 
      main='Maximum activity day',  lwd=2) 




K=4
t1_index=c(10,20,30,50,100)
result_1_step<-matrix(nrow=length(t1_index), ncol=ncol(step.list[[67]]))

for(k in c(2,5)){
  t1=t1_index[k]
  r1=0.2
  Up_step = prepare_up(step.list[[67]],FUN=sqBound_es,t1,r1) ## Ensemble pen result
  
  a1= cluster_zits(K=K,Up=log(Up_step), FUN=l1_mean ,t1,r1,is.median=T)
  yhat1=vector()
  for(j in 1:K){
    yhat1[a1$e[[j]]]=j
  }
  result_1_step[k,]=yhat1
  
}

table( result_1_step[k,])
dim(result_1_step)



###############  Plot 67th  ############### 
tp = seq(as.POSIXct("00:00:00",format="%H:%M"),
         as.POSIXct("23:59:00",format="%H:%M"), by = "1 min")

par(mfrow=c(1,K))
for(k in c(5)){
  
  aa=(table(result_1_step[k,]))
  
  for(j in 1:K){
  plot( tp, apply(step.list[[67]][,which(result_1_step[k,]==j)],1, mean), type='l', ylab='', xlab='time', col=1, 
        main=paste('Group - ',j, '(',aa[j] ,')' ),  lwd=2, ylim=c(0,26)) 
  }
  
  #legend('topright', paste(leg_tex), col=order(aa, decreasing=T), lty=1, lwd=2)
  
}

par(mfrow=c(1,2))
for(k in c(2,5)){
  
  aa=(table(result_1_step[k,]))
  
  
  plot( tp, apply(step.list[[67]][,which(result_1_step[k,]==1)],1, mean), type='l', ylab='', xlab='time', col=1, 
        main=paste('Thickness - ',t1_index[k] ),  lwd=2, ylim=c(0,30)) 
  #axis(side=1, at=seq(1, length(tp), length.out=7), labels=tp[seq(1, length(tp), length.out=7)])
  for(j in order(aa, decreasing=T)){ 
    lines( tp, apply(step.list[[67]][,which(result_1_step[k,]==j)],1, mean), col=j,type='l',lwd=2 )
  }
  leg_tex<-vector()
  for(j in order(aa, decreasing=T)){
    leg_tex=c(leg_tex,aa[j]) }
  legend('topright', paste(leg_tex), col=order(aa, decreasing=T), lty=1, lwd=2)
  
}



########
#Exploratory plots
########

# heatmap
people = c()
for(i in 1:length(step.list)){
  n=0
  for(p in step.list[[i]]){ n = n + length(p)/1440 }
  people=c(people,rep(i,n))
}
people.save=people


group.by.ppl = matrix(0,ncol=length(step.list), nrow=K)

ttt=0
aa= (table(result_step[2,]))
for(i in order(aa, decreasing=T)){
  ttt=ttt+1
  e30=which(result_step[2,]==i)
  group.by.ppl[ttt,as.integer(names(table(people[e30])))] = table(people[e30])
}

ttt=0
aa= (table(result_step[5,]))
group.by.ppl2 = matrix(0,ncol=length(step.list), nrow=K)
for(i in order(aa, decreasing=T)){
  ttt=ttt+1
  e100=which(result_step[5,]==i)
  group.by.ppl2[ttt,as.integer(names(table(people[e100])))] = table(people[e100])
}

par(mfrow=c(1,2))
par(mar=c(3.5,3.5,1.5,1))
par(oma=c(3,2,1,1.5))
par(mgp=c(2,.5,0))
image(group.by.ppl,axes = F,ylab="individual",xlab="(a) thickness 20")
axis(1, at=seq(0,1,0.2), labels=sprintf("group%d",1:K))
image(group.by.ppl2,axes = F,xlab="(b) thickness 100")
axis(1, at=seq(0,1,0.2), labels=sprintf("group%d",1:K))




par(mfrow=c(1,2))
case=  order( apply(step.list[[67]][1:(60*3),],2,sum) , decreasing = T)[1]
#case=  order( apply(step.list[[67]],2,sum) )[32]
#  plot( tp, step.list[[67]][,case] , type='l', ylab='', xlab='time', col=1,        main='',  lwd=1) 

  plot( tp, step.list[[67]][,case] , type='l', ylab='', xlab='time', col=1, 
        main='',  lwd=1, ylim=c(-10,180)) 
    polygon( c(tp,rev(tp)) ,c(sqBound_es(step.list[[67]][,case], 20, max), rep(0, 1440)),col="skyblue" )
    lines(tp,step.list[[67]][,case],type='l')

    plot( tp, step.list[[67]][,case] , type='l', ylab='', xlab='time', col=1, 
          main='',  lwd=1, ylim=c(-10,180)) 
    polygon( c(tp,rev(tp)) ,c(sqBound_es(step.list[[67]][,case], 100, max), rep(0, 1440)),col="skyblue" )
    lines(tp,step.list[[67]][,case],type='l')
    
  
  

###############  cluster validation  ############### 

clust_val<-matrix(nrow=4, ncol=4)
df2=dist(t(data))
for(j in 1:2){
 clust_stats <- cluster.stats(df2,  result_step[c(2,5)[j],]  )
 clust_val[j,]=round(c( clust_stats $dunn, clust_stats$avg.silwidth, clust_stats $entropy ,clust_stats $ch),3)
 }
 clust_stats <- cluster.stats((df2),  yhat_funfem)
 clust_val[3,]=round(c( clust_stats $dunn, clust_stats$avg.silwidth, clust_stats $entropy,clust_stats $ch ),3)
 
 
 clust_stats <- cluster.stats((df2),  yhat_hddc)
 clust_val[4,]=round(c( clust_stats $dunn, clust_stats$avg.silwidth, clust_stats $entropy,clust_stats $ch ),3)
 
 
 colnames(clust_val)=c('Dunn','silwidth','entropy','CH')
 rownames(clust_val)=c('tau20','tau100','FEM','HDDC')
 clust_val

