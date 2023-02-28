######################################################
###### Piotr Fryzlewicz & Hernando Ombao 2009 ###########
######################################################

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
  
  if(model=='4B'){
    x1=x2=matrix(0, nrow=n, ncol=TT)
    
    for(i in 1:n){
      sigma=0.5
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
  
  return(x)
  
}
