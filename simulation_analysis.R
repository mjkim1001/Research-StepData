library(fda)
library(funFEM)
library(forecast)
library(funHDDC)
library(mclust)
library(cvTools)
library(clue)

source("./simulation_models.R")
source("./functions.R")

total.n.sim=100

TP_result1=TP_result2=matrix(nrow=5, ncol=total.n.sim)
TP_maj_result1<-TP_maj_result2<-vector()
TP_cv_result1=TP_cv_result2<-vector()
ccr.funhddc=ccr.funfem =ccr.dtw=vector()

ar.TP_result1=ar.TP_result2=matrix(nrow=5, ncol=total.n.sim)
ar.TP_maj_result1<-ar.TP_maj_result2<-vector()
ar.TP_cv_result1=ar.TP_cv_result2<-vector()
ar.funhddc=ar.funfem =ar.dtw=vector()

cv_opt=vector()

time_d<-matrix(nrow=total.n.sim, ncol=9)

model_type='1A'## 1A , 2A, 2B, block, 4A(sigma 0.1, 0.5)

###################################################
######### Simulation model figures ############ 
###################################################
n=100
TT=500

data_temp= data_generation(n=n, TT=TT, model=model_type)
print(length(which(data_temp[[1]]==0))/length(data_temp[[1]]) )
par(mfrow=c(1,2), mai=c(0.5,0.8,0.8,0.5))
for(i in c(1: length(data_temp)) ){
  sam=sample(c(1:n), 1)
  plot(data_temp[[i]][sam,], type='b',pch=20, cex=0.7, ylab='', xlab='',
       main=paste('Model 1 (g=', i,')' , sep=''), ylim=c(0,7) )
  abline(v=54,col=2,lty=2,lwd=2)
  abline(v=128,col=2,lty=2,lwd=2)
}

###################################################
######### Cross Validation for Thickness ############ 
###################################################
for(sim in 1:total.n.sim){
  print(paste(sim,'th interation'))
  set.seed(sim)
  n=100
  TT=500
  skip_sim <- FALSE
  #Generate Dataset
  data_temp= data_generation(n=n, TT=TT, model=model_type)
  data=data_temp[[1]]
  for(l in 2:length(data_temp)){
    data=rbind(data ,data_temp[[l]] )
  }
  data=t(data)  ##TT*(N=n*K) matrix
  dim(data)
  
  K=length(data_temp)
  perm_num=nrow(perm(1:K))
  true_cl=matrix(ncol=perm_num, nrow=n*K)
  for(j in 1:perm_num){
    true_cl[,j]=rep(perm(1:K)[j,] , each=n)
  }
  
  tryCatch({
    
    start_time <- Sys.time()
    ## various t value ####
    t1_index=c(10,20,30,50,100,150)
    cv_error<-matrix(nrow=2, ncol=length(t1_index))
    for(k in 1:length(t1_index)){
      
      t1=t1_index[k]
      r1=0.1
      ## Cross-validation for thickness code ####
      K_cv=5
      folds <- cvFolds(ncol(data), K=K_cv)
      
      temp_error_ccr<-temp_error_adj<-vector()
      for(i in 1:K_cv){
        train <- data[,folds$subsets[folds$which != i] ] #Set the training set
        validation <- data[,folds$subsets[folds$which == i] ] #Set the validation set

        Up_sim = prepare_up(train,FUN=sqBound_es,t1,r1) 
        
        a1= cluster_zits(K=K,Up=log(Up_sim), FUN=l1_mean ,t1,r1,is.median=T)
        
        Up_test = prepare_up(validation,FUN=sqBound_es,t1,r1)
        
        pmeasure=matrix(nrow=ncol(Up_test) , ncol=K)
        for(clstr in 1:K){
          for(j in 1:ncol(Up_test)){
            pmeasure[j,clstr] = l1_mean(a1$pivot[,clstr],log(Up_test[,j]))
          }
        }
        
        yhat=apply(pmeasure,1, which.min)
        #assign each validation set to the pivots derived using train set
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
    
    #Optimal Thickness based on 5 fold CV.
    opt_k=which.max( cv_error[2,])
    t1=t1_index[opt_k]
    print(t1)
    cv_opt[sim]=t1
    r1=0.1
    Up_sim = prepare_up(data,FUN=sqBound_es,t1,r1)
    
    a1= cluster_zits(K=K,Up=log(Up_sim), FUN=l1_mean ,t1,r1,is.median=1)
    yhat=vector()
    for(l in 1:K ){
      yhat[a1$clusters[[l]]]=l
    }
    ccr=vector()
    for(j in 1:perm_num){
      ccr[j]=sum(diag(table(true_cl[,j],yhat)))/(ncol(data))
    }
    
    TP_cv_result1[sim]=max(ccr)
    ar.TP_cv_result1[sim]= adjustedRandIndex(true_cl[,1],yhat)
    
    end_time <- Sys.time()
    time_d[sim,6 ] = round(difftime(end_time, start_time, units = "secs"), 3)
    
  }, error = function(e) { 
    print("TP error")
    print(e)
    skip_sim <<- TRUE
  })
  
  if (skip_sim == TRUE) {
    next
  }
}

###################################################
######### Simulation running 100 times ############ 
###################################################
#### Use parallel computing and compute time
library(parallel)
nCores = detectCores()

s <- system.time({
  run_sim <- mclapply(1:total.n.sim, function(sim){
    print(paste(sim,'th interation'))
    set.seed(sim)
    n=100
    TT=500
    #Generate Dataset
    data_temp= data_generation(n=n, TT=TT, model=model_type)
    data=data_temp[[1]]
    for(l in 2:length(data_temp)){
      data=rbind(data ,data_temp[[l]] )
    }
    data=t(data)  ##TT*(N=n*K) matrix
    dim(data)
     
    K=length(data_temp)
    perm_num=nrow(perm(1:K))
    true_cl=matrix(ncol=perm_num, nrow=n*K)
    for(j in 1:perm_num){
      true_cl[,j]=rep(perm(1:K)[j,] , each=n)
    }
    
    tryCatch({
    start_time <- Sys.time()
    
    t1_index=c(10,20,30,50,100,150)
    
    #par(mfrow=c(1,2)) 
    #plot(t1_index, cv_error[1,],  ylab='Adj. rand index', xlab='thickness', type='b')
    #plot(t1_index, cv_error[2,], type='b', ylab='CCR', xlab='thickness' )
    
    opt_k=4
    t1=t1_index[opt_k]
    print(t1)
    #cv_opt[sim]=t1
    r1=0.1
    Up_sim = prepare_up(data,FUN=sqBound_es,t1,r1)
    
    a1= cluster_zits(K=K,Up=log(Up_sim), FUN=l1_mean ,t1,r1,is.median=1)
    yhat=vector()
    for(l in 1:K ){
      yhat[a1$clusters[[l]]]=l
    }
    ccr=vector()
    for(j in 1:perm_num){
      ccr[j]=sum(diag(table(true_cl[,j],yhat)))/(ncol(data))
    }
    
    end_time <- Sys.time()
    return(list(obj = a1$obj ,ccr = max(ccr),aRand = adjustedRandIndex(true_cl[,1],yhat), time =round(difftime(end_time, start_time, units = "secs"), 3) ))  
    }, error = function(e) { 
      print("TP error")
      print(e)
      skip_sim <<- TRUE
    })
    
    if (skip_sim == TRUE) {
      next
    }
  }, mc.cores = nCores-1)
})
run_sim1 <- unlist(run_sim)
run_sim1 <- matrix(run_sim1,ncol=4,byrow = T)
colnames(run_sim1) <- c("obj","ccr","arand","time")
run_sim1 = as_tibble(run_sim1)

###################################################
######## Comparison methods running 100 times ####### 
###################################################
library(foreach)
library(doParallel)
nCores = detectCores()
cl <-makeCluster(nCores-1)
registerDoParallel(cl)

funfem_result<-foreach(sim = 1:total.n.sim, .combine = rbind) %dopar% {
  set.seed(sim)
  n=100
  TT=1500
  #Generate Dataset
  data_temp= data_generation(n=n, TT=TT, model=model_type)
  data=data_temp[[1]]
  for(l in 2:length(data_temp)){
    data=rbind(data ,data_temp[[l]] )
  }
  data=t(data)  ##TT*(N=n*K) matrix
  
  K=length(data_temp)
  perm_num=nrow(perm(1:K))
  true_cl=matrix(ncol=perm_num, nrow=n*K)
  for(j in 1:perm_num){
    true_cl[,j]=rep(perm(1:K)[j,] , each=n)
  }

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
  ####################################### 
  ###### funHDDC
  # res.uni <- funHDDC(fdobj,K=K,model="AkBkQkDk",init="kmeans",threshold=0.2)
  # yhat= res.uni$class
  ###### DTW
  #pc.l2 <-tsclust(t(data), k = K,  distance = "sdtw", centroid = "pam",seed = 3247, trace = FALSE,  control = partitional_control(nrep = 5L))
  #yhat=cl_class_ids(cl_medoid(pc.l2))
  #######################################
  #Compute CCR
  ccr=vector()
  for(j in 1:perm_num){
    ccr[j]=sum(diag(table(true_cl[,j],yhat)))/(ncol(data))
  }
  end_time <- Sys.time()
  
  (data.frame(ccr.funfem = max(ccr), ar.funfem= adjustedRandIndex(true_cl[,1],yhat), time = round(difftime(end_time, start_time, units = "secs"), 3)))
}

stopCluster(cl)
