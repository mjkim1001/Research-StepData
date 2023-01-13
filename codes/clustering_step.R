#example
Up = prepare(sqBound_es,30,0.2,smooth=1)
#distance functions
log_mean<-function(ux,uy){ #l1 distance
  return(mean(abs(ux-uy)))
}
sumsquare<-function(ux,uy){ #l2 distance
  return(sum((ux-uy)^2))
}

################################################################
################## Clustering             ######################
################################################################

###function cluster
# Up: upper thick pen bounds
# K : number of clusters
# FUN : distance functions
cluster<-function(Up,K,FUN,maxIter=500){ #logUp 를 받아 L1 optimization
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
      save(elements, pivot, file = sprintf("./fin/smoothma/%d_run%d.rda",K,tau))
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
  }
}

step_result <- function(Up,K,N=12,is.log=1){#N번 돌려서 결과를 본다. #FUN : sqBound_es/ma
  start_time = Sys.time()
  obj_vec= c()
  result= list()
  if(is.log){#Our method : L1 optimization after log transform
    for(i in 1:N){
      result[[i]] = cluster(log(Up),K,FUN=log_mean) ###min/max measure #log(Up)
      print(summary(result[[i]]$e))
      obj_vec[i] = result[[i]]$obj
    }
  }else{#L2 optimization with TPT
    for(i in 1:N){
      result[[i]] = cluster(Up,K,sumsquare) ###min/max measure,  #Up가 들어가야됨
      print(summary(result[[i]]$e))
      obj_vec[i] = result[[i]]$obj
    }
  }
  end_time = Sys.time()
  print(end_time-start_time)
  return(list(ls = result, obj = obj_vec)) #obj_vec는 길이 N짜리 obj vector.
} 

#run with different parameters
library(clusterSim)
fun.step <- function(x,k){
  print(k)
  a = step_result(x, k, N=20, is.log = 1)
  clusters = rep(NA, ncol(step))
  for(c in 1:k){
    clusters[(a$ls[[ which.min(a$obj) ]]$e)[[c]]] = c
  }
  return(list(cluster = clusters))
}
for(tau in c(30,100)){
  Up = prepare(sqBound_es,tau,0.2,smooth=1)
  gap.step= clusGap(Up,fun.step,10)
  plot(gap.step)
}

################################################################
################## Figures and Results   ######################
################################################################

######step figure
#result and obj_vec is the clustering result list with tau=30
#result2 and obj_vec2 is the clustering result list with tau=100
e30 = result[[which.min(obj_vec)]]$e 
e100 = result2[[which.min(obj_vec2)]]$e
for(i in 1:6){
  e100[[i]] = (result2[[which.min(obj_vec2)]]$e)[[c(5,1,2,3,6,4)[i]]]
  #manual 하게 e30 그룹 결과와 순서 비슷하도록 ordering 조정
}
#number of days (Table 4.1)
summary(e30)
summary(e100)
#Figure 4.1 and 4.2
plot_means(e30)
plot_means(e30,mode=2,Up)
plot_means(e100)
plot_means(e100,2,Up2)

#group by people (Figure 4.3)
people = c()
for(i in 1:length(step_long)){
  n=0
  for(p in step_long[[i]]){ n = n + length(p)/1440 }
  people=c(people,rep(i,n))
}
people.save=people
#people
par(mfrow=c(1,2))
par(mar=c(3.5,3.5,1.5,1))
par(oma=c(3,2,1,1.5))
par(mgp=c(2,.5,0))#축이름_눈금_제목
group.by.ppl = matrix(0,ncol=length(step_long), nrow=6)
for(i in 1:6){
  group.by.ppl[i,as.integer(names(table(people[e30[[i]]])))] = table(people[e30[[i]]])
}
group.by.ppl2 = matrix(0,ncol=length(step_long), nrow=6)
for(i in 1:6){
  group.by.ppl2[i,as.integer(names(table(people[e100[[i]]])))] = table(people[e100[[i]]])
}
image(group.by.ppl,axes = F,ylab="individual",xlab="(a) thickness 30")
axis(1, at=seq(0,1,0.2), labels=sprintf("grp%d",1:6))
image(group.by.ppl2,axes = F,xlab="(b) thickness 100")
axis(1, at=seq(0,1,0.2), labels=sprintf("grp%d",1:6))

#stepcount ( Table 4.1 )
numstep=c()
for(e in e30){
  numstep = c(numstep, mean(apply(step[,e],2,sum))  )
}
numstep2=c()
for(e in e100){
  numstep2 = c(numstep2, mean(apply(step[,e],2,sum))  )
}
numstep2
#weekend% ( Table 4.1 )
weekend.percent = matrix(NA,nrow=2,ncol=6)
for(i in 1:6){
  weekend.percent[1,i] = mean(is.wkend2[e30[[i]]])
  weekend.percent[2,i] = mean(is.wkend2[e100[[i]]])
}
weekend.percent

#function for figure 4.1 and 4.2
plot_means<-function(elements,mode=1,Up=NULL,median=FALSE){
  time = seq(as.POSIXct("00:00:00",format="%H:%M"),
             as.POSIXct("23:59:00",format="%H:%M"), by="1 min")
  #if(!is.null(Up)){Up = exp(Up)}
  K=length(elements)
  if(mode==1){#평균 한번에
    plot(time,step[,1],ylim =c(0,ifelse(Te==1440,25,5)),type = "n", xlab="",ylab="")
    for(i in 1:K){
      if(median){
        lines(time,apply(step[,elements[[i]]],1,median),col=rainbow(K)[i])
      }else{
        lines(time,apply(step[,elements[[i]]],1,mean),col=rainbow(K)[i])
      }
    }
    legend("topright", legend=sapply(elements,length), col=rainbow(K),lty=rep(1,K),title="#elements")
  }else if(mode==2){#하나씩
    par(mfrow=c(3,2))
    par(mar=c(3.5,3.5,1,1))#space around figure!!
    par(oma=c(3,2,1,1.5)) #outer margin
    par(mgp=c(2.5,1,0)) #축과의 거리
    for(i in 1:K){
      if(median){
        plot( time, apply(exp(Up[,elements[[i]]]),1,median), xlab=sprintf("(%s)",letters[i]) ,ylab='',type='l',lty=2,ylim=c(0,40+tau/2))
        lines(time, apply(step[,elements[[i]]],1,median))
      }else{
        plot( time, exp(apply(Up[,elements[[i]]],1,median)), xlab=sprintf("(%s)",letters[i]) ,ylab='',type='l',lty=2,ylim=c(0,40+tau/2))
        lines(time, apply(step[,elements[[i]]],1,mean))}
    }
  }
}
