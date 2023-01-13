library(fda)
library(funFEM)
library(forecast)
library(funHDDC)
library(mclust)
library(cvTools)
library(clue)

source("/Users/yaejilim/Library/CloudStorage/GoogleDrive-yaeji.lim@stat.cau.ac.kr/My Drive/Papers/MIN_JI_clustering/Stepdata_code/yj_code/functions.R")

total.n.sim=20

TP_result1=TP_result2=matrix(nrow=5, ncol=total.n.sim)
TP_maj_result1<-TP_maj_result2<-vector()
TP_cv_result1=TP_cv_result2<-vector()
ccr.funhddc=ccr.funfem =ccr.dtw=vector()

ar.TP_result1=ar.TP_result2=matrix(nrow=5, ncol=total.n.sim)
ar.TP_maj_result1<-ar.TP_maj_result2<-vector()
ar.TP_cv_result1=ar.TP_cv_result2<-vector()
ar.funhddc=ar.funfem =ar.dtw=vector()


time_d<-matrix(nrow=total.n.sim, ncol=9)

model_type='1A'## 1A , 2A, 2B, block, 4A(sigma 0.1, 0.5)

for(sim in 1:total.n.sim){
  print(paste(sim,'th interation'))
  set.seed(sim)
  n=100
  TT=1500
  
  skip_sim <- FALSE 
  data_temp= data_generation(n=n, TT=TT, model=model_type)
 
 #  
 #  print(length(which(data_temp[[1]]==0))/length(data_temp[[1]]) )
 #  par(mfrow=c(1,2), mai=c(0.5,0.8,0.8,0.5))
 #  for(i in c(1: length(data_temp)) ){
 #  	sam=sample(c(1:n), 1)
 #    plot(data_temp[[i]][sam,], type='b',pch=20, cex=0.7, ylab='', xlab='',
 #         main=paste('Model 1 (g=', i,')' , sep=''), ylim=c(0,7) ) 
 #  abline(v=54,col=2,lty=2,lwd=2)
 # abline(v=128,col=2,lty=2,lwd=2)
 #  }

  
  
  data=data_temp[[1]]
  for(l in 2:length(data_temp)){
    data=rbind(data ,data_temp[[l]] )
  }
  data=t(data)  ##TT*n
  dim(data)
   
  K=length(data_temp)
  perm_num=nrow(perm(1:K))
  true_cl=matrix(ncol=perm_num, nrow=n*K)
  for(j in 1:perm_num){
    true_cl[,j]=rep(perm(1:K)[j,] , each=n)
  }
  
  tryCatch({
    

  ## various t value ####
  t1_index=c(10,20,30,50,100)
  maj_yhat1=maj_yhat2=matrix(nrow=ncol(data), ncol=length(t1_index))
  
  for(k in 1:length(t1_index)){
    start_time <- Sys.time()
    
    t1=t1_index[k]
    r1=0.1
    Up_sinu = prepare_sim(data,FUN=sqBound_es,t1,r1, smooth=0, smooth.ma.num = smooth.ma.num) ## Ensemble pen result
    
    a1= cluster_sim(K=K,data,Up=log(Up_sinu), FUN=log_mean ,t1,r1,is.median=1)
    yhat=vector()
    for(l in 1:K ){
      yhat[a1$e[[l]]]=l
    }
    
    maj_yhat1[,k]=yhat
    
    ccr=vector()
    for(j in 1:perm_num){
      ccr[j]=sum(diag(table(true_cl[,j],yhat)))/(ncol(data))
    }
    
    TP_result1[k,sim]=max(ccr)
    ar.TP_result1[k,sim]=adjustedRandIndex(true_cl[,1],yhat)
    

    end_time <- Sys.time()
    print(paste0("time : ", 
                 round(difftime(end_time, start_time, units = "secs"), 3),
                 " secs"))
    time_d[sim,k ] = round(difftime(end_time, start_time, units = "secs"), 3)
    
    
  }
  
  
  
  
 
  
  ## cv code ####
  start_time <- Sys.time()
  
  t1_index=c(10,20,30,50,100,150)
  cv_error<-matrix(nrow=2, ncol=length(t1_index))
  for(k in 1:length(t1_index)){
    
    t1=t1_index[k]
    r1=0.1
    
    K_cv=5
    folds <- cvFolds(ncol(data), K=K_cv)
    
    temp_error_ccr<-temp_error_adj<-vector()
    for(i in 1:K_cv){
      train <- data[,folds$subsets[folds$which != i] ] #Set the training set
      validation <- data[,folds$subsets[folds$which == i] ] #Set the validation set
      
      
      Up_sinu = prepare_sim(train,FUN=sqBound_es,t1,r1, smooth=0, smooth.ma.num = smooth.ma.num) ## Ensemble pen result
      
      a1= cluster_sim(K=K,train,Up=log(Up_sinu), FUN=log_mean ,t1,r1,is.median=1)
      for(l in 1:K){
        a1$pivot[,l]= apply( log(Up_sinu)[ ,a1$e[[l]]],1, median)}
      
      Up_test = prepare_sim(validation,FUN=sqBound_es,t1,r1, smooth=0, smooth.ma.num = smooth.ma.num) ## Ensemble pen result
      
      
      pmeasure=matrix(nrow=ncol(Up_test) , ncol=K)
      for(clstr in 1:K){
        for(j in 1:ncol(Up_test)){
          pmeasure[j,clstr] = log_mean(a1$pivot[,clstr],log(Up_test[,j]))
        }
      }
      
      yhat=apply(pmeasure,1, which.min)
      # print(yhat)
      true_cl2=true_cl[ folds$subsets[folds$which == i] ,]
      
      ccr=vector()
      for(j in 1:perm_num){
        ccr[j]=sum(diag(table(true_cl2[,j],yhat)))/(length(yhat))
      }
      
      temp_error_ccr[i]=max(ccr)
      temp_error_adj[i]=adjustedRandIndex(true_cl2[,1],yhat)
      
    }
    
    cv_error[1,k]=mean(temp_error_adj)
    cv_error[2,k]=mean(temp_error_ccr)
    
  }
 
  #par(mfrow=c(1,2)) 
  #plot(t1_index, cv_error[1,],  ylab='Adj. rand index', xlab='thickness', type='b')
  #plot(t1_index, cv_error[2,], type='b', ylab='CCR', xlab='thickness' )
  
  opt_k=which.max( cv_error[2,])
  
  t1=t1_index[opt_k]
  r1=0.1
  Up_sinu = prepare_sim(data,FUN=sqBound_es,t1,r1, smooth=0, smooth.ma.num = smooth.ma.num) ## Ensemble pen result
  
  
  a1= cluster_sim(K=K,data,Up=log(Up_sinu), FUN=log_mean ,t1,r1,is.median=1)
  yhat=vector()
  for(l in 1:K ){
    yhat[a1$e[[l]]]=l
  }
  ccr=vector()
  for(j in 1:perm_num){
    ccr[j]=sum(diag(table(true_cl[,j],yhat)))/(ncol(data))
  }
  
  TP_cv_result1[sim]=max(ccr)
  ar.TP_cv_result1[sim]= adjustedRandIndex(true_cl[,1],yhat)
  
  end_time <- Sys.time()
  print(paste0("time : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  time_d[sim,6 ] = round(difftime(end_time, start_time, units = "secs"), 3)
  
  
  # a2= cluster_sim(K=K,data,Up=(Up_sinu), FUN=sumsquare ,t1,r1,is.median=2)
  # yhat=vector()
  # for(l in 1:K ){
  #   yhat[a2$e[[l]]]=l
  # }
  # 
  # ccr=vector()
  # for(j in 1:perm_num){
  #   ccr[j]=sum(diag(table(true_cl[,j],yhat)))/(ncol(data))
  # }
  # TP_cv_result2[sim]=max(ccr)
  # ar.TP_cv_result2[sim]= adjustedRandIndex(true_cl[,1],yhat)
  # 
  # 
  # 
  
  
  
  }, error = function(e) { 
    print("TP error")
    print(e)
    skip_sim <<- TRUE
  })
  
  if (skip_sim == TRUE) {
    next
  }
  ################## Functional Clustering  ######################
  ################################################################
  start_time <- Sys.time()
  
  basis <- create.fourier.basis(c(0, nrow(data)), nbasis=11)
  fdobj <- smooth.basis(1:nrow(data), data,basis)$fd
  cl.sinu = funFEM(fdobj, K=K, init="kmeans", model="all")
  elements = list()
  for(c in 1:K){
    # get the point of cluster c
    elements[[c]] <- which(cl.sinu$cls == c)
  }
  yhat=vector()
  for(l in 1:K ){
    yhat[elements[[l]]]=l
  }
  
  ccr=vector()
  for(j in 1:perm_num){
    ccr[j]=sum(diag(table(true_cl[,j],yhat)))/(ncol(data))
  }
  
  ccr.funfem[sim] = max(ccr)
  ar.funfem[sim] = adjustedRandIndex(true_cl[,1],yhat)
  
  end_time <- Sys.time()
  print(paste0("time : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  time_d[sim,7 ] = round(difftime(end_time, start_time, units = "secs"), 3)
  
  
  start_time <- Sys.time()
  res.uni <- funHDDC(fdobj,K=K,model="AkBkQkDk",init="kmeans",threshold=0.2)
  
  yhat= res.uni$class
  ccr=vector()
  for(j in 1:perm_num){
    ccr[j]=sum(diag(table(true_cl[,j],yhat)))/(ncol(data))
  }
  ccr.funhddc[sim] = max(ccr)
  ar.funhddc[sim] = adjustedRandIndex(true_cl[,1],yhat)
  
  end_time <- Sys.time()
  print(paste0("time : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  time_d[sim,8 ] = round(difftime(end_time, start_time, units = "secs"), 3)
  
  
  
  ### DTW cluster
  start_time <- Sys.time()
  pc.l2 <-tsclust(t(data), k = K,  distance = "sdtw", centroid = "pam",seed = 3247, trace = FALSE,  control = partitional_control(nrep = 5L))
  #int.ts=which.max(sapply(pc.l2, cvi)[1,])
  #dtw_r=pc.l2 [[int.ts]]
  yhat=cl_class_ids(cl_medoid(pc.l2))
  ccr=vector()
  for(j in 1:perm_num){
    ccr[j]=sum(diag(table(true_cl[,j],yhat)))/(ncol(data))
  }
  ccr.dtw[sim] = max(ccr)
  ar.dtw[sim] = adjustedRandIndex(true_cl[,1],yhat)
  
  end_time <- Sys.time()
  print(paste0("time : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  time_d[sim,9 ] = round(difftime(end_time, start_time, units = "secs"), 3)
  
  
}

colnames(time_d)=c("TPT_10","TPT_20","TPT_30","TPT_50","TPT_100","TPT_CV","FUNFEM","FUNHDDC","DTW")
apply(time_d,2, mean)



#save(TP_result1,  TP_cv_result1, ccr.funfem ,ccr.funhddc,ar.TP_result1, ar.TP_cv_result1,  ccr.dtw,
   #  ar.funfem,ar.funhddc , ar.dtw,    file='/Users/yaejilim/Desktop/4A_3.RData')

#t1_index=c( 10,  20  ,30,  50 ,100, 150)

re1=round( c(  apply(TP_result1, 1, mean, na.rm=T)[c(2,3,4)],  mean(TP_cv_result1, na.rm=T), 
                mean(ccr.funfem, na.rm=T),mean(ccr.funhddc, na.rm=T) , mean(ccr.dtw, na.rm=T)       ),3)

re2=round( c(  apply(TP_result1, 1, sd, na.rm=T)[c(2,3,4)],  sd(TP_cv_result1, na.rm=T), 
               sd(ccr.funfem, na.rm=T),sd(ccr.funhddc, na.rm=T) ,sd(ccr.dtw, na.rm=T)   ) ,3)


paste(re1[1],'(',re2[1],') &',re1[2],'(',re2[2],')& ',re1[3],'(',re2[3],')&',  re1[4],'(',re2[4],')&',  re1[5],'(',re2[5],')&' ,  re1[6],'(',re2[6],')&',
      re1[7],'(',re2[7],')&',sep='' )


re1=round( c(  apply(ar.TP_result1, 1, mean, na.rm=T)[c(2,3,4)],  mean(ar.TP_cv_result1, na.rm=T), 
               mean(ar.funfem, na.rm=T),mean(ar.funhddc, na.rm=T),mean(ar.dtw, na.rm=T)   ),3)
re2=round( c(  apply(ar.TP_result1, 1, sd, na.rm=T)[c(2,3,4)],  sd(ar.TP_cv_result1, na.rm=T), 
               sd(ar.funfem, na.rm=T),sd(ar.funhddc, na.rm=T)   ,sd(ar.dtw, na.rm=T)   ),3)


paste(re1[1],'(',re2[1],') &',re1[2],'(',re2[2],')& ',re1[3],'(',re2[3],')&',  re1[4],'(',re2[4],')&',  re1[5],'(',re2[5],')&' ,  re1[6],'(',re2[6],')&' 
, re1[7],'(',re2[7],')&',sep='' )





re1=round( c(  apply(TP_result1,1, mean, na.rm=T)[c(3,5)] ,  mean(TP_cv_result1, na.rm=T),  mean(TP_maj_result1, na.rm=T), 
               apply(TP_result2,1, mean, na.rm=T)[c(3,5)] ,  mean(TP_cv_result2, na.rm=T),  mean(TP_maj_result2, na.rm=T),
               mean(ccr.funfem, na.rm=T),mean(ccr.funhddc, na.rm=T) ,mean(ccr.dtw, na.rm=T)   ),3)
re2=round( c(  apply(TP_result1,1, sd, na.rm=T)[c(3,5)] ,  sd(TP_cv_result1, na.rm=T),  sd(TP_maj_result1, na.rm=T), 
               apply(TP_result2,1, sd, na.rm=T)[c(3,5)] ,  sd(TP_cv_result2, na.rm=T),  sd(TP_maj_result2, na.rm=T),
               sd(ccr.funfem, na.rm=T),sd(ccr.funhddc, na.rm=T),   sd(ccr.dtw, na.rm=T)  ),3)

paste(re1[1],'(',re2[1],') &',re1[2],'(',re2[2],')& ',re1[3],'(',re2[3],')&',  re1[4],'(',re2[4],')&',  re1[5],'(',re2[5],')&' ,  re1[6],'(',re2[6],')&' 
      ,  re1[7],'(',re2[7],')&',  re1[8],'(',re2[8],')&' ,  re1[9],'(',re2[9],')&' ,  re1[10],'(',re2[10],')' 
)



re1=round( c(  apply(ar.TP_result1,1, mean, na.rm=T)[c(3,5)] ,  mean(ar.TP_cv_result1, na.rm=T),  mean(ar.TP_maj_result1, na.rm=T), 
               apply(ar.TP_result2,1, mean, na.rm=T)[c(3,5)] ,  mean(ar.TP_cv_result2, na.rm=T),  mean(ar.TP_maj_result2, na.rm=T),
               mean(ar.funfem, na.rm=T),mean(ar.funhddc, na.rm=T)),3)
re2=round( c(  apply(ar.TP_result1,1, sd, na.rm=T)[c(3,5)] ,  sd(ar.TP_cv_result1, na.rm=T),  sd(ar.TP_maj_result1, na.rm=T), 
               apply(ar.TP_result2,1, sd, na.rm=T)[c(3,5)] ,  sd(ar.TP_cv_result2, na.rm=T),  sd(ar.TP_maj_result2, na.rm=T),
               sd(ar.funfem, na.rm=T),sd(ar.funhddc, na.rm=T)),3)
paste(re1[1],'(',re2[1],') &',re1[2],'(',re2[2],')& ',re1[3],'(',re2[3],')&',  re1[4],'(',re2[4],')&',  re1[5],'(',re2[5],')&' ,  re1[6],'(',re2[6],')&' 
      ,  re1[7],'(',re2[7],')&',  re1[8],'(',re2[8],')&' ,  re1[9],'(',re2[9],')&' ,  re1[10],'(',re2[10],')' 
)


