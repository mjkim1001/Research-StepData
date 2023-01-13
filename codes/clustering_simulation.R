################################################################
################## Data generator        ######################
################################################################
############## Simulation 1
Te=1024
sigma = c(0.01,0.1,0.7,1)
sinu = rep(0,Te)
for(sd in sigma){
  for(n in 1:50){
    sinu = cbind(sinu,abs(sin(5*(0:(Te-1))/Te)+rnorm(Te,sd = sd)))
  }
}
sinu = sinu[,-1]
############## Figure 5.1
par(mfrow=c(2,2)) 
for(i in 1:4){ 
  plot(sinu[,50*(i-1)+4], type='l',xlab="",ylab="", ylim=c(0,4)) 
} 
################################################################
############## Simulation 2
h = c(-1,1,-1,1)*runif(4,0,20) 
h [5] = -sum(h[1:4]) 
Te=512 
blk = rep(0,Te) 
for(k in 1:4){ 
  for(n in 1:50){ 
    xi = runif(5, (k-1)/5, (k+1)/5 ) 
    temp = rnorm(Te,sd = 3) 
    for( j in 1:5){ 
      temp = temp + h[j]*(1 + sign((0:(Te-1))/(Te) - xi[j]))/2 
    } 
    blk = cbind(blk,abs(temp)) 
  } 
} 
blk = blk[,-1]
############## Figure 5.2
par(mfrow=c(2,2)) 
for(i in 1:4){ 
  plot(blk[,50*(i-1)+1], type='l',xlab="",ylab="", ylim=c(0,40)) 
} 
################################################################
############## Simulation 3
Te=512 
pat = NULL 
#grp1 
h =c( c(1,-1,0)*runif(3,0,30), c(1,-1,-1,1)*runif(4,0,20)) 
h [3] = -sum(h[1:2]) 
h[8] = -sum(h[4:7]) 
for(n in 1:50){ 
  xi =c( runif(3,0/5,1/5), runif(5,2/5,3/5)) 
  temp = rnorm(Te,sd = 5) 
  for( j in 1:8){ 
    temp = temp + h[j]*(1 + sign((0:(Te-1))/(Te) - xi[j]))/2 
  } 
  pat = cbind(pat,abs(temp)) 
} 
#grp2
h =c( c(1,-1,-1,1,0)*runif(5,0,10), c(1,-1)*runif(2,0,5)) 
h [5] = -sum(h[1:4]) 
h[8] = -sum(h[6:7]) 
for(n in 1:50){ 
  xi =c( runif(5,0/5,2/5), runif(3,2/5,4/5)) 
  temp = rnorm(Te,sd = 3)
  for( j in 1:8){
    temp = temp + h[j]*(1 + sign((0:(Te-1))/(Te) - xi[j]))/2 
  } 
  pat = cbind(pat,abs(temp)) 
} 
#grp3
h =c( c(1,-1,1,0)*runif(4,0,20), c(1,-1,0)*runif(3,0,15),c(1,-1)*runif(2,0,20)) 
h [4] = -sum(h[1:3]) 
h[7] = -sum(h[5:6]) 
h[10]= -sum(h[8:9]) 
for(n in 1:50){ 
  xi =c( runif(4,1/5,2/5), runif(3,2/5,4/5),runif(3,4/5,1)) 
  temp = rnorm(Te,sd = 5) 
  for( j in 1:10){ 
    temp = temp + h[j]*(1 + sign((0:(Te-1))/(Te) - xi[j]))/2 
  } 
  pat = cbind(pat,abs(temp)) 
} 
############## Figure 5.3
par(mfrow=c(1,3)) 
for(i in 1:3){ 
  plot(pat[,50*(i-1)+4], type='l',xlab="",ylab="", ylim=c(0,30)) 
} 
################################################################

################################################################
################## Iterative Clustering   ######################
################## functions              ######################
################################################################
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
  print(end_time-start_time)
  return(Up)
}
Up_sinu = prepare_sim(sinu,sqBound_es,30,r1=0.01,smooth.ma.num = 3)

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
      print(itr)
      end_time = Sys.time()
      #save(elements, pivot, file = sprintf("./fin/%d_run%d.rda",K,tau))
      print(end_time-start_time)
      obj=0
      for(i in 1:K){
        obj = obj+sum(pmeasure[elements[[i]],i])
      }
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

#각각의 data, distance function 등에 대해 N번 반복 실행하여 결과를 리스트로 반환
simul_result <- function(K,data,FUN,t1,r1=0.5, N=12,smooth.ma.num = 3,is.l1=1, is.log = 1,is.up=1){#N번 돌려서 결과를 본다. #FUN : sqBound_es/ma
  start_time = Sys.time()
  obj_vec= c()
  sim_list= list()
  if(is.up){
    Up_sinu = prepare_sim(data,FUN,t1,r1,smooth = !is.null(smooth.ma.num),smooth.ma.num = smooth.ma.num)
    print(dim(Up_sinu))
    if(is.log){
      Up_sinu = log(Up_sinu)
    }
  }else{
    Up_sinu=data
  }
  if(is.l1){
    for(i in 1:N){
      sim_list[[i]] = cluster_sim(K,data,Up_sinu,log_mean,t1,r1,is.median = is.l1) ###min/max measure
      print(summary(sim_list[[i]]$e))
      obj_vec[i] = sim_list[[i]]$obj
    }
  }else{
    for(i in 1:N){
      sim_list[[i]] = cluster_sim(K,data,Up_sinu,sumsquare,t1,r1,is.median = is.l1) ###min/max measure
      print(summary(sim_list[[i]]$e))
      obj_vec[i] = sim_list[[i]]$obj
    }
  }
  
  end_time = Sys.time()
  print(end_time-start_time)
  return(list(ls = sim_list, obj = obj_vec))
}

################################################################
################## Run!                   ######################
################################################################
sinu_cl <-1:50
for(k in 1:3){
  sinu_cl= cbind(sinu_cl,(50*k+1):(50*(k+1)))
}#clustering answer for simulation set1
blk_cl <- sinu_cl #clustering answer for simulation set2
pat_cl <- sinu_cl[,1:3] #clustering answer for simulation set3

ccr.mat = matrix(NA,nrow=100,ncol=5)
for(i in 1:100){
  a = simul_result(4,sinu,sqBound_es,30,r1=0.01,N=20,smooth.ma.num = 3,is.l1 = 0, is.log = 0) # up->kmeans(TPL2)
  b = simul_result(4,sinu,sqBound_es,30,r1=0.01,N=20,smooth.ma.num = 3,is.l1 = 1) # up ->log->l1  #GOOD (TPMA0)
  c = simul_result(4,sinu,sqBound_es,30,r1=0.01,N=20,smooth.ma.num = 3,is.l1 = 1, is.log = 0)#up -> l1 #GOOD (TPL1)
  d = simul_result(4,sinu,sqBound_es,30,r1=0.01,N=20,smooth.ma.num = 3,is.l1 = 0, is.log = 0,is.up = 0) #sinu->l2 (L2)
  e = simul_result(4,sinu,sqBound_es,30,r1=0.01,N=20,smooth.ma.num = 3,is.l1 = 1, is.log = 0, is.up=0)#sinu -> l1 (L1)
  ccr.mat[i,] = c(ccr(sinu_cl,a$ls[[which.min(a$obj)]]$e ),
                  ccr(sinu_cl,b$ls[[which.min(b$obj)]]$e ),
                  ccr(sinu_cl,c$ls[[which.min(c$obj)]]$e ),
                  ccr(sinu_cl,d$ls[[which.min(d$obj)]]$e ),
                  ccr(sinu_cl,e$ls[[which.min(e$obj)]]$e ))
}
ccr.blk = matrix(NA,nrow=100,ncol=5)
for(i in 1:100){  
  a = simul_result(4,blk,sqBound_es,30,r1=0.1,N=20,smooth.ma.num = 3,is.l1 = 0, is.log = 0) # up->kmeans  
  b = simul_result(4,blk,sqBound_es,30,r1=0.1,N=20,smooth.ma.num = 3,is.l1 = 1) # up ->log->l1  #GOOD 
  c = simul_result(4,blk,sqBound_es,30,r1=0.1,N=20,smooth.ma.num = 3,is.l1 = 1, is.log = 0)#up -> l1 #GOOD  
  d = simul_result(4,blk,sqBound_es,30,r1=0.1,N=20,smooth.ma.num = 3,is.l1 = 0, is.log = 0,is.up = 0) #sinu->l2  
  e = simul_result(4,blk,sqBound_es,30,r1=0.1,N=20,smooth.ma.num = 3,is.l1 = 1, is.log = 0, is.up=0)#sinu -> l1  
  ccr.blk[i,] = c(ccr(blk_cl,a$ls[[which.min(a$obj)]]$e ),  
                  ccr(blk_cl,b$ls[[which.min(b$obj)]]$e ),  
                  ccr(blk_cl,c$ls[[which.min(c$obj)]]$e ),  
                  ccr(blk_cl,d$ls[[which.min(d$obj)]]$e ),  
                  ccr(blk_cl,e$ls[[which.min(e$obj)]]$e )) 
} 
ccr.pat = matrix(NA,nrow=100,ncol=5)
for(i in 1:100){  
  a = simul_result(3,pat,sqBound_es,30,r1=0.1,N=20,smooth.ma.num = 3,is.l1 = 0, is.log = 0) # up->kmeans  
  b = simul_result(3,pat,sqBound_es,30,r1=0.1,N=20,smooth.ma.num = 3,is.l1 = 1) # up ->log->l1  #GOOD 
  c = simul_result(3,pat,sqBound_es,30,r1=0.1,N=20,smooth.ma.num = 3,is.l1 = 1, is.log = 0)#up -> l1 #GOOD  
  d = simul_result(3,pat,sqBound_es,30,r1=0.1,N=20,smooth.ma.num = 3,is.l1 = 0, is.log = 0,is.up = 0) #sinu->l2  
  e = simul_result(3,pat,sqBound_es,30,r1=0.1,N=20,smooth.ma.num = 3,is.l1 = 1, is.log = 0, is.up=0)#sinu -> l1  
  ccr.pat[i,] = c(ccr(pat_cl,a$ls[[which.min(a$obj)]]$e ),  
                  ccr(pat_cl,b$ls[[which.min(b$obj)]]$e ),  
                  ccr(pat_cl,c$ls[[which.min(c$obj)]]$e ),  
                  ccr(pat_cl,d$ls[[which.min(d$obj)]]$e ),  
                  ccr(pat_cl,e$ls[[which.min(e$obj)]]$e )) 
}

################################################################
################## Functional Clustering  ######################
################################################################
library(fda)
library(funFEM)
basis <- create.fourier.basis(c(0, 1024), nbasis=11)
fdobj <- smooth.basis(1:1024, sinu,basis)$fd

basis <- create.bspline.basis(c(0, 512), nbasis=11)
fdobj2 <- smooth.basis(1:512, blk,basis)$fd

basis <- create.bspline.basis(c(0, 512), nbasis=11)
fdobj3 <- smooth.basis(1:512, pat,basis)$fd

ccr.funfem = matrix(NA,nrow=100,ncol=3)
for(i in 1:100){
  cl.sinu = funFEM(fdobj, K=4, init="kmeans", model="all")
  elements = list()
  for(c in 1:K){
    # get the point of cluster c
    elements[[c]] <- which(cl.sinu$cls == c)
  }
  ccr.funfem[i,1] = ccr(sinu_cl,elements )
}
for(i in 1:100){
  cl.blk = funFEM(fdobj2, K=4, init="kmeans", model="all")
  elements = list()
  for(c in 1:K){
    # get the point of cluster c
    elements[[c]] <- which(cl.blk$cls == c)
  }
  ccr.funfem[i,2] = ccr(blk_cl,elements )
}
for(i in 1:100){
  cl.pat = funFEM(fdobj3, K=3, init="kmeans", model="all")
  elements = list()
  for(c in 1:K){
    # get the point of cluster c
    elements[[c]] <- which(cl.pat$cls == c)
  }
  ccr.funfem[i,3] = ccr(pat_cl,elements )
}

ccr.sinu = cbind(ccr.sinu,ccr.funfem[,1])
ccr.blk = cbind(ccr.blk, ccr.funfem[,2] )
ccr.pat = cbind(ccr.pat, ccr.funfem[,3] )

sinu.mean= round(apply(ccr.sinu,2,mean),2)
sinu.std= round(sqrt(apply(ccr.sinu,2,var)),2)
blk.mean= round(apply(ccr.blk,2,mean),2)
blk.std= round(sqrt(apply(ccr.blk,2,var)),2)
pat.mean = round(apply(ccr.pat,2,mean),2)
pat.std= round(sqrt(apply(ccr.pat,2,var)),2)
#Table 5.1 Result
rbind( sinu.mean,blk.mean ,pat.mean)
rbind( sinu.std,blk.std ,pat.std)