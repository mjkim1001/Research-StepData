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

#################################################################
################################################################
### distance function used as an argument for clustering_zits ###
################################################################
################################################################

l1_mean<-function(ux,uy){
  return(mean(abs(ux-uy)))
}
l2_mean<-function(ux,uy){
  return(mean((ux-uy)^2))
}
#################################################################
################################################################
################## Clustering ZITS algorithm ######################
################################################################
################################################################

# function returns T times N matrix of upper thick pen bound
# data : data
# FUN : gets pen bound functions such as cirBound or sqBound_es
# t1: thickness, r1: scaling factor
prepare_up<-function(data,FUN,t1,r1=0.5){ 
  start_time = Sys.time()
  Te = dim(data)[1]
  N = dim(data)[2]
  
  Up = rep(0,Te) # matrix of (1)
  for(i in 1:N){
    Up=cbind(Up, matrix(FUN(data[,i],t1,max,r1),nrow=Te))
  }
  Up = Up[,-1]
  end_time = Sys.time()
  return(Up)
}

##### Input
# Up: upper thick pen bounds
# K : number of clusters
# FUN : distance functions
# To use our clustering algorithm, pass log(Up) as the input for the "Up" argument 
# and apply L1 optimization by using FUN<-l1_mean
##### Output 
# e : a list of elements assigned to each cluster
# obj: objective function value for the optimization problem
# pivot : T * K matrix, each column is a cluster medoid for each cluster
######
cluster_zits <- function(K, Up, FUN, t1, r1 = 0.5, maxIter = 500, is.median = 1) { 
  start_time <- Sys.time() 
  Te <- nrow(Up)
  N <- ncol(Up)
  clusters <- list() 
  # A length-K list, each element contains a vector of days assigned to the corresponding cluster
  
  # Initialize clusters
  ####################################################################################
  initial_centroids <- sample(1:N, K)
  pmeasure <- rep(0, N)  
  pivot <- Up[, initial_centroids]
  
  for(i in 1:K){
    distances_to_centroid <- rep(0, N)
    for(j in 1:N){
      distances_to_centroid[j] <- FUN(pivot[, i], Up[, j])
    }
    pmeasure <- cbind(pmeasure, distances_to_centroid)
    clusters[[i]] <- c(initial_centroids[i])
  }
  pmeasure <- pmeasure[, -1]
  
  ########### Clustering
  itr <- 1
  while(TRUE) {
    # Assign
    cluster_assignments <- apply(pmeasure, 1, which.min)
    
    # Update clusters
    cluster_index <- 1
    has_converged <- TRUE
    for(c in 1:K){
      # Get the points in the current cluster
      points_in_cluster <- which(cluster_assignments == c)
      
      if(length(points_in_cluster) == 0) {
        # remove empty clusters
        clusters <- clusters[-cluster_index]
        pmeasure <- pmeasure[, -cluster_index]
        pivot <- pivot[, -cluster_index]
        has_converged <- FALSE
        next
      } else if(length(points_in_cluster) <= 1) {
        if ((length(points_in_cluster) != length(clusters[[cluster_index]])) ||
            sum(points_in_cluster != clusters[[cluster_index]])) {
          has_converged <- FALSE
        }
        clusters[[cluster_index]] <- points_in_cluster
      } else {
        if ((length(points_in_cluster) != length(clusters[[cluster_index]])) ||
            sum(points_in_cluster != clusters[[cluster_index]])) {
          has_converged <- FALSE
        }
        clusters[[cluster_index]] <- points_in_cluster
      }
      cluster_index <- cluster_index + 1
    }
    K <- length(clusters)
    if(has_converged){# check for convergence
      end_time <- Sys.time()
      obj=0
      for(i in 1:K){
        obj = obj+sum(pmeasure[clusters[[i]],i])
      }
      return(list(clusters = clusters, obj = obj, pivot = pivot))
    }
    itr <- itr + 1
    if(itr == maxIter) {
      print("Max iteration reached")
      end_time <- Sys.time()
      return(clusters)
    }
    # Update prototype
    for(c in 1:K){
      points_in_cluster = clusters[[c]]
      if(length(points_in_cluster)==1){
        pivot[,c] = Up[,points_in_cluster]
      }else{
        # compute the mean/median point of all points in cluster c
        pivot[,c] <- apply(Up[,points_in_cluster], 1, if (is.median) median else mean)
      }
      for(j in 1:N){
        pmeasure[j,c] = FUN(pivot[,c],Up[,j])
      }
    }
  }
} 
