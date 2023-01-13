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