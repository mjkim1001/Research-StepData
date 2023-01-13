######################################################
###### Piotr Fryzlewicz & Hernando Ombao 2009 ###########
######################################################
### 
library(dtwclust)
library(cluster)
library(ZIM)
perm <- function(v) {
  n <- length(v)
  if (n == 1) v
  else {
    X <- NULL
    for (i in 1:n) X <- rbind(X, cbind(v[i], perm(v[-i])))
    X
  }
}

data_generation=function(n=100, TT=256, model='1A'){
   
  if(model=='4A'){
    x1=x2=matrix(0, nrow=n, ncol=TT)
 
	 for(i in 1:n){
	 	sigma=0.1
	 	
	 	mu1=runif(1, 3,4)
	 	mu2=runif(1, 2,10)
	 	lambda1=rnorm(TT,mu1 , sigma)
	 	lambda2=rnorm(TT, mu2, sigma)
	 	si=runif(1, 0.4, 0.7)
	 	x1[i,]= rzip(n = TT,  lambda = lambda1, omega = si)
    	x2[i,]= rzip(n = TT, lambda = lambda2,  omega = si)
    	}
    	print(length( which(x1==0) ) / length(x1) )
    	#par(mfrow=c(1,2))
    	# plot(x1[1,], type='b', ylim=range(x1,x2))
     #    plot(x2[1,], type='b', ylim=range(x1,x2))	 
   x<-list()
    x[[1]]=x1
    x[[2]]=x2
  }
  
  if(model=='1A'){
    x1=x2=matrix(0, nrow=n, ncol=TT)
    
    
    phi1<-phi2<-vector()
    
    phi1[1:53]= 0.8
    phi1[54:128]= -0.9
    phi1[129:TT]=0.8
    phi2[1:TT]= -0.81
    
    for(i in 1:n){
      e  = rnorm(TT ,0,2)
      
      for(t in 3:TT){
        x1[i,t] = phi1[t] * x1[i,t-1] + phi2[t] * x1[i,t-2] +e[t]
      }
    }
    
    
    phi1<-phi2<-vector()
    
    phi1[1:53]= 0.8
    phi1[54:128]=1.6
    phi1[129:TT]=0.8
    phi2[1:TT]= -0.81
    
    for(i in 1:n){
      
      e  = rnorm(TT ,0,2)
      for(t in 3:TT){
        x2[i,t] = phi1[t] * x2[i,t-1] + phi2[t] * x2[i,t-2] +e[t]
      }
      
    }
    x1[x1<0]=0
    x2[x2<0]=0
    x<-list()
    x[[1]]=x1
    x[[2]]=x2
       
  }
  
  ### case 1B
  if(model=='1B'){
    
    x1=x2=matrix(0, nrow=n, ncol=TT)
    
    for(i in 1:n){
      phi1<-phi2<-vector()
      
      phi1[1:53]= 0.8
      phi1[54:700]= -0.9
      phi1[701:TT]= 0.8
      phi2[1:TT]= -0.81
      
      
      e  = rnorm(TT ,0,1)
      
      for(t in 3:TT){
        x1[i,t] = phi1[t] * x1[i,t-1] + phi2[t] * x1[i,t-2] +e[t]
      }
      
      
      phi1<-phi2<-vector()
      
      phi1[1:53]= 0.8
      phi1[54:700]= 1.6
      phi1[701:TT]= 0.8
      phi2[1:TT]= -0.81
      
      
      e  = rnorm(TT ,0,1)
      for(t in 3:TT){
        x2[i,t] = phi1[t] * x2[i,t-1] + phi2[t] * x2[i,t-2] +e[t]
      }
      
    }
    x1[x1<0]=0
    x2[x2<0]=0
    x<-list()
    x[[1]]=x1
    x[[2]]=x2
    
    }
  
  
  
  if(model=='1A_mean0'){
    x1=x2=matrix(0, nrow=n, ncol=TT)
    
    
    phi1<-phi2<-vector()
    
    phi1[1:53]= 0.8
    phi1[54:128]= -0.9
    phi1[129:TT]=0.8
    phi2[1:TT]= -0.81
    
    for(i in 1:n){
      e  = rnorm(TT ,0,1)
      
      for(t in 3:TT){
        x1[i,t] = phi1[t] * x1[i,t-1] + phi2[t] * x1[i,t-2] +e[t]
      }
    }
    
    
    phi1<-phi2<-vector()
    
    phi1[1:53]= 0.8
    phi1[54:128]=1.6
    phi1[129:TT]=0.8
    phi2[1:TT]= -0.81
    
    for(i in 1:n){
      
      e  = rnorm(TT ,0,1)
      for(t in 3:TT){
        x2[i,t] = phi1[t] * x2[i,t-1] + phi2[t] * x2[i,t-2] +e[t]
      }
      
    }
    x<-list()
    x[[1]]=x1
    x[[2]]=x2
    
  }
  

  
  
  
  ### case 2A
  
  if(model=="2A"){
    
    x1=x2=matrix(0, nrow=n, ncol=TT)
    
    for(i in 1:n){
      phi1<-phi2<-vector()
      
      for(t in 1:TT){
        phi1[t]= -0.8*(1-0.7*cos(pi*t/TT))
      }
      
      phi2[1:TT]= -0.81
      
      
      e  = rnorm(TT ,0,1)
      
      for(t in 3:TT){
        x1[i,t] = phi1[t] * x1[i,t-1] + phi2[t] * x1[i,t-2] +e[t]
      }
      
      
      phi1<-phi2<-vector()
      
      for(t in 1:TT){
        phi1[t]= -0.8*(1-0.001*cos(pi*t/TT))
      }
      phi2[1:TT]= -0.81
      
      
      e  = rnorm(TT ,0,1)
      for(t in 3:TT){
        x2[i,t] = phi1[t] * x2[i,t-1] + phi2[t] * x2[i,t-2] +e[t]
      }
      
    }
    x1[x1<0]=0
    x2[x2<0]=0
    x<-list()
    x[[1]]=x1
    x[[2]]=x2
    
  }
  ### case 2B
  if(model=="2B"){
    
    
    x1=x2=matrix(0, nrow=n, ncol=TT)
    
    for(i in 1:n){
      phi1<-phi2<-vector()
      
      for(t in 1:TT){
        phi1[t]= -0.8*(1-0.7*cos(pi*t/TT))
      }
      
      phi2[1:TT]= -0.81
      
      
      e  = rnorm(TT ,0,1)
      
      for(t in 3:TT){
        x1[i,t] = phi1[t] * x1[i,t-1] + phi2[t] * x1[i,t-2] +e[t]
      }
      
      
      phi1<-phi2<-vector()
      
      for(t in 1:TT){
        phi1[t]= -0.8*(1-0.1*cos(pi*t/TT))
      }
      phi2[1:TT]= -0.81
      
      
      e  = rnorm(TT ,0,1)
      for(t in 3:TT){
        x2[i,t] = phi1[t] * x2[i,t-1] + phi2[t] * x2[i,t-2] +e[t]
      }
      
    }
    x1[x1<0]=0
    x2[x2<0]=0
    x<-list()
    x[[1]]=x1
    x[[2]]=x2
  }
  
  ### case 3
  if(model=="3"){
    
    x<-list()
    sigma = c(0.01,0.1,0.7,1)
    sinu = rep(0,TT)
    for(sd in c(1:length(sigma))){
      x[[sd]]=matrix(0, nrow=n, ncol=TT)
      for(i in 1:n){
        x[[sd]][i,]= abs(sin(5*(0:(TT-1))/TT)+rnorm(TT,sd = sigma[sd]))
      }
    }
 
  }
  
  
  ### case 2A
  
  if(model=="2A_mean0"){
    
    x1=x2=matrix(0, nrow=n, ncol=TT)
    
    for(i in 1:n){
      phi1<-phi2<-vector()
      
      for(t in 1:TT){
        phi1[t]= -0.8*(1-0.7*cos(pi*t/TT))
      }
      
      phi2[1:TT]= -0.81
      
      
      e  = rnorm(TT ,0,1)
      
      for(t in 3:TT){
        x1[i,t] = phi1[t] * x1[i,t-1] + phi2[t] * x1[i,t-2] +e[t]
      }
      
      
      phi1<-phi2<-vector()
      
      for(t in 1:TT){
        phi1[t]= -0.8*(1-0.001*cos(pi*t/TT))
      }
      phi2[1:TT]= -0.81
      
      
      e  = rnorm(TT ,0,1)
      for(t in 3:TT){
        x2[i,t] = phi1[t] * x2[i,t-1] + phi2[t] * x2[i,t-2] +e[t]
      }
      
    }
    x<-list()
    x[[1]]=x1
    x[[2]]=x2
    
  }
  
  if(model=="block"){
  h = c(-1,1,-1,1)*runif(4,0,20) 
  h [5] = -sum(h[1:4]) 

  x<-list()
  for(k in 1:4){ 
    blk = rep(0,TT) 
    for(i in 1:n){ 
      xi = runif(5, (k-1)/5, (k+1)/5 ) 
      temp = rnorm(TT,sd = 3) 
      for( j in 1:5){ 
        temp = temp + h[j]*(1 + sign((0:(TT-1))/(TT) - xi[j]))/2 
      } 
    #  temp[which(temp<0.1)]=0
      temp[which(temp<0)]=0
       blk = cbind(blk,abs(temp)) 
    } 
    blk = blk[,-1]
    x[[k]]=t(blk)
    } 

 
  }
  return(x)
  
}


log_mean<-function(ux,uy){
  return(mean(abs(ux-uy)))
}
sumsquare<-function(ux,uy){
  return(sum((ux-uy)^2))
}

ccr<-function(asw,rslt){
  N = dim(asw)[2] * dim(asw)[1]
  c=NULL
  for(e in rslt){
    cr = rep(0,ncol(asw))
    for(i in 1:dim(asw)[2]){
      cr[i] = sum(e %in% asw[,i])
    }
    c = cbind(c,cr)
  }
  od = rbind(t(c(1:ncol(asw))) , allPerms(ncol(asw)))
  ccr_vec = rep(0, nrow(od))
  for(i in 1:nrow(od)){
    for(j in 1:ncol(c)){
      ccr_vec[i] = ccr_vec[i] + c[od[i,j],j]
    }
  }
  return(max(ccr_vec)/N)
}

library(forecast)


# function returns T*N matrices of upper thick pen bound
# data : data
# FUN : gets pen bound functions such as cirBound or sqBound_es
# t1: thickness, r1: scaling factor
# smooth : binary ( 1 or 0) factor for smoothing
prepare_sim<-function(data,FUN,t1,r1=0.5,smooth=1,smooth.ma.num=5){ 
  start_time = Sys.time()
  Te = dim(data)[1]
  N = dim(data)[2]
  if(smooth){
    for(i in 1:N){
      data[,i] = smooth.ma(data[,i],smooth.ma.num)
    }
  }
  Up = rep(0,Te) # matrix of (1)
  for(i in 1:N){
    Up=cbind(Up, matrix(FUN(data[,i],t1,max,r1),nrow=Te))
  }
  Up = Up[,-1]
  end_time = Sys.time()
 # print(end_time-start_time)
  return(Up)
}


prepare_sim_down<-function(data,FUN,t1,r1=0.5,smooth=1,smooth.ma.num=5){ 
  start_time = Sys.time()
  Te = dim(data)[1]
  N = dim(data)[2]
  if(smooth){
    for(i in 1:N){
      data[,i] = smooth.ma(data[,i],smooth.ma.num)
    }
  }
  Up = rep(0,Te) # matrix of (1)
  for(i in 1:N){
    Up=cbind(Up, matrix(FUN(data[,i],t1,min,r1),nrow=Te))
  }
  Up = Up[,-1]
  end_time = Sys.time()
  #print(end_time-start_time)
  return(Up)
}

#################################################################
################################################################
################## Pen shapes definition  ######################
################################################################
################################################################
# t : thickness, x: vector of time series, 
# FUN : max for an upper bound and min for an lower bound
# r : scaling factor

#square pen
sqBound<-function(x,t,FUN,r=1){
  t<-as.integer(t/2) ; s<-FUN(1,-1)*r
  bound<-c()
  n<-length(x)
  
  bound[1] = FUN(x[1:(1+t)])
  for(i in 2:n){
    prev = i-1
    if(i<=(t+1) || (i<=n-t && bound[prev] != x[i-t-1]) ){
      bound[i] <- FUN(bound[prev], x[i+t])
    }else{
      bound[i]<-FUN(x[(i-t):min(n,i+t)])
    }
  }
  return(bound+s*t)
}
#Circle pen
cirBound<-function(x,t,FUN,r=1){
  tr<-as.integer(t/2);  bound<-c()
  s<-FUN(1,-1)*r; n<-length(x)
  
  for(i in 1:n){
    temp<-c()
    if( (i-tr) < 1 ) { K <- (1:(i+tr))
    }else if( (i+tr) > n ){ K <- ((i-tr):n)
    }else{ K <- ((i-tr):(i+tr))}
    
    for (k in K){
      temp<- append(temp,x[k]+s*sqrt((t^2)/4-(k-i)^2))
    }
    bound[i]<-FUN(temp)
  }
  return(bound)
}
#Ensemble square pen
sqBound_es<-function(x,t,FUN,r=1){
  s<-FUN(1,-1)*r
  n<-length(x)
  Bound<-data.frame(1:n)
  
  for(esmb in 0:t){
    bound<-c()
    bound[1] = FUN(x[1:(1+t-esmb)])
    for(i in 2:n){
      prev = i-1
      
      if(i<=(esmb+1) || (i<=n-t+esmb && bound[prev] != x[i-esmb-1]) ){
        bound[i] <- FUN(bound[prev], x[i+t-esmb])
      }else{
        bound[i]<-FUN(x[(i-esmb):min(n,i+t-esmb)])
      }
    }
    Bound[,esmb+1]<-bound
  }
  bound<-apply(Bound,1,mean)
  return(bound+s*as.integer(t/2))
}

#moving average
smooth.ma <- function(x,h){  
  if(!h%%2){  
    y = as.vector(ma(x,h,F))  
    n = length(x)  
    for(i in 1:(h/2)){  
      y[n-i+1] = mean(x[(n-i+2 - h/2):n])  
      y[i] = mean(x[1:(i+h/2)])  
    }  
    return(y)  
  }  
  y = as.vector(ma(x,h))  
  n = length(x)  
  for(i in 1:((h-1)/2)){  
    y[i] = mean(x[1:(i+ (h-1)/2)])  
    y[n-i+1] = mean(x[(n-i+1-(h-1)/2 ):n])  
  }  
  return(y)  
}

# Up: upper thick pen bounds
# K : number of clusters
# FUN : distance functions
##log(Up)를 받아서 euclidean L1 optimization => FUN: log_mean
##Up를 받아서 euclidean L2 optimization => FUN: sumsquare, is.median=0
cluster_sim<-function(K,data,Up,FUN,t1,r1=0.5,maxIter=500, is.median=1){ 
  start_time = Sys.time() 
  Te = dim(data)[1]
  N = dim(data)[2]
  ####################################################################################3
  elements=list() # length K list, each contains vector of days assigned to the cluster
  
  #initialize cluster
  ###############3
  samp = sample(1:N,K)
  pmeasure = rep(0,N)  #각 pivot 날짜에 대해 모든 날들과의 measure를 저장
  means = data[,samp]
  pivot = Up[,samp]
  for(i in 1:K){
    mvec = rep(0,N)
    for(j in 1:N){
      mvec[j] = FUN(pivot[,i],Up[,j])
    }
    pmeasure = cbind(pmeasure,mvec)
    elements[[i]] = c(samp[i])
  }
  pmeasure = pmeasure[,-1]
  ###########clustering
  itr=1
  while(TRUE)
  {
    #assign
    clusters = apply(pmeasure,1,which.min)
    #update : elements와 means를 update
    ci = 1
    diff = TRUE
    means = NULL
    for(c in 1:K){
      # get the point of cluster c
      group <- which(clusters == c)
      if(length(group)==0){
        elements = elements[-ci]
        pmeasure = pmeasure[,-ci]
        pivot = pivot[,-ci]
        diff = FALSE
        next
      }else if(length(group)<=1){
        if ((length(group) != length(elements[[ci]]))|| sum(group != elements[[ci]])){
          diff = FALSE
        }
        elements[[ci]] = group
        means = cbind(means, Up[,group])
      }else{
        if ((length(group) != length(elements[[ci]]))|| sum(group != elements[[ci]])){
          diff = FALSE
        }
        elements[[ci]] = group
        # compute the mean point of all points in cluster c
        means = cbind(means, apply(Up[,group],1,mean))
      }
      ci = ci+1
    }
    K = length(elements)
    if(diff){
      #print(itr)
      end_time = Sys.time()
      #save(elements, pivot, file = sprintf("./fin/%d_run%d.rda",K,tau))
      #print(end_time-start_time)
      obj=0
      for(i in 1:K){
        obj = obj+sum(pmeasure[elements[[i]],i])
      }
      return(list(e = elements, obj = obj, pivot=pivot ))
    }
    itr=itr+1
    if(itr==maxIter){
      print("max iteration reached")
      end_time = Sys.time()
      #print(end_time-start_time)
      return(elements)
    }
    #plot
    #  print(summary(elements))
    #   plot(0,0,xlim =c(0,Te),ylim =c(0,ifelse(Te==1440,90,5)),type = "n")
    #    for(i in 1:K){
    #      lines(means[,i])
    #      lines(pivot[,i], lty=2)
    #  }
    #legend("topright", legend=sapply(elements,length), col=rainbow_hcl(K),lty=rep(1,K),title="#elements")
    
    #update prototype
    for(clstr in 1:K){
      group = elements[[clstr]]
      if(length(group)==1){
        pivot[,clstr] = Up[,group]
      }else{
        if(is.median){
          pivot[,clstr] = apply(Up[,group],1,median)
        }else{
          pivot[,clstr] = apply(Up[,group],1,mean)
        }
      }
      for(j in 1:N){
        pmeasure[j,clstr] = FUN(pivot[,clstr],Up[,j])
      }
    }
  }
}

prepare_sim<-function(data,FUN,t1,r1=0.5,smooth=1,smooth.ma.num=5){ 
  start_time = Sys.time()
  Te = dim(data)[1]
  N = dim(data)[2]
  if(smooth){
    for(i in 1:N){
      data[,i] = smooth.ma(data[,i],smooth.ma.num)
    }
  }
  Up = rep(0,Te) # matrix of (1)
  for(i in 1:N){
    Up=cbind(Up, matrix(FUN(data[,i],t1,max,r1),nrow=Te))
  }
  Up = Up[,-1]
  end_time = Sys.time()
  #print(end_time-start_time)
  return(Up)
}



################################################################
################## Real step 
##################################  ################## ##############################
prepare<-function(step_long, people, FUN,t1,r1=0.5, smooth=1){ 
  start_time = Sys.time()
  Up = rep(0,1440) # will be a matrix of upper bounds
  n = length(step_long)
  pb <- progress_bar$new(total = n)
  if(smooth){
    for(people in step_long){ # draw thick pen on concatenated(long) data 
      for(days in 1:length(people)){
        X <-smooth.ma(people[[days]],5)
        Up=cbind(Up, matrix(FUN(X,t1,max,r1),nrow=1440))
      }
      pb$tick()
    }
  }else{
    for(people in step_long){ # draw thick pen on concatenated(long) data 
      for(days in 1:length(people)){
        X <-people[[days]]
        Up=cbind(Up, matrix(FUN(X,t1,max,r1),nrow=1440))
      }
      pb$tick()
    }
  }
  Up = Up[,-1]
  end_time = Sys.time()
  print(end_time-start_time)
  return(Up)
}




################################################################
################## Clustering             ######################
##### function returns 1440*N matrices of upper thick pen bound
# FUN : gets pen bound functions such as cirBound or sqBound_es
# t1: thickness, r1: scaling factor
# smooth : binary ( 1 or 0) factor for smoothing
################################################################
############################################################

###function cluster
# Up: upper thick pen bounds
# K : number of clusters
# FUN : distance functions
cluster<-function(Up,K,FUN,maxIter=400){ #logUp 를 받아 L1 optimization
  start_time = Sys.time()
  Te = dim(Up)[1]
  N = dim(Up)[2]
  ####################################################################################3
  elements=list() # length K list, each contains vector of days assigned to the cluster
  
  #initialize cluster
  ###############3
  samp = sample(1:N,K)
  pmeasure = rep(0,N)  #각 pivot 날짜에 대해 모든 날들과의 measure를 저장
  means = Up[,samp]
  pivot = Up[,samp]
  for(i in 1:K){
    mvec = rep(0,N)
    for(j in 1:N){
      mvec[j] = FUN(pivot[,i],Up[,j])
    }
    pmeasure = cbind(pmeasure,mvec)
    elements[[i]] = c(samp[i])
  }
  pmeasure = pmeasure[,-1]
  ###########clustering start!
  itr=1
  while(TRUE)
  {
    #assign
    clusters = apply(pmeasure,1,which.min)
    #update : elements와 means를 update
    ci = 1
    diff = TRUE
    means = NULL
    for(c in 1:K){
      # get the point of cluster c
      group <- which(clusters == c)
      if(length(group)==0){ #remove one cluster
        elements = elements[-ci]
        pmeasure = pmeasure[,-ci]
        pivot = pivot[,-ci]
        diff = FALSE
        next
      }else if(length(group)<=1){
        if ((length(group) != length(elements[[ci]]))|| sum(group != elements[[ci]])){
          diff = FALSE
        }
        elements[[ci]] = group
        means = cbind(means, Up[,group])
      }else{
        if ((length(group) != length(elements[[ci]]))|| sum(group != elements[[ci]]) ){
          diff = FALSE
        }
        elements[[ci]] = group
        # compute the mean point of all points in cluster c
        means = cbind(means, apply(Up[,group],1,mean))
      }
      ci = ci+1
    }
    K = length(elements)
    if(diff){
      print(itr)
      end_time = Sys.time()
     # save(elements, pivot, file = sprintf("./fin/smoothma/%d_run%d.rda",K,tau))
      obj =0
      for(i in 1:K){
        obj = obj+sum(pmeasure[elements[[i]],i])
      }
      print(end_time-start_time)
      return(list(e = elements, obj = obj ))
    }
    itr=itr+1
    if(itr==maxIter){
      print("max iteration reached")
      end_time = Sys.time()
      print(end_time-start_time)
      return(elements)
    }
    #plot
    print(summary(elements))
    plot(0,0,xlim =c(0,Te),ylim =c(0,ifelse(Te==1440,90,5)),type = "n")
    for(i in 1:K){
      lines(exp(means[,i]),col=rainbow_hcl(K)[i])
      lines(exp(pivot[,i]), lty=2,col=rainbow_hcl(K)[i])
    }
    legend("topright", legend=sapply(elements,length), col=rainbow_hcl(K),lty=rep(1,K),title="#elements")
    
    #update prototype
    for(clstr in 1:K){
      group = elements[[clstr]]
      if(length(group)==1){
        pivot[,clstr] = Up[,group]
      }else{
        pivot[,clstr] = apply(Up[,group],1,median)
      }
      for(j in 1:N){
        pmeasure[j,clstr] = FUN(pivot[,clstr],Up[,j])
      }
    }
  }; print(itr)
}




thick_draw <-function(data,tau,FUN,r=1,time=FALSE,smooth=0, smooth.num=5,zeroBound=TRUE,ymin=NULL,ymax=NULL, xl,col='skyblue'){
  if(smooth){
    data = smooth.ma(data,smooth.num)
  }
  up = FUN(data, tau,max,r=r)
  if(zeroBound){low = rep(0,length(data))}else{low = FUN(data, tau,min,r=r)}
  if(time){
    x = seq(as.POSIXct("00:00:00",format="%H:%M"),
            as.POSIXct("23:59:00",format="%H:%M"), by = "1 min")
  }else{
    x = 1:length(data)
  }
  if(is.null(ymin)){
    plot(x,data,type='n',main=xl,ylab="", ylim=c(min(low),max(up)), xlab='')
  }else{
    plot(x,data,type='n',main=xl,ylab="",ylim=c(ymin,ymax), xlab='')
  }
  polygon( c(x,rev(x)) ,c(up, low),col=col )
  lines(x,data,type='l')
}

