#################################################################
################################################################
################## Read Data and          ######################
################## Data Preprocessing     ######################
################################################################
################################################################

read_step<-function(){
  step<-list()
  idx=1
  for(file in list.files(path="/Users/yaejilim/Library/CloudStorage/GoogleDrive-yaeji.lim@stat.cau.ac.kr/My Drive/Papers/MIN_JI_clustering/Stepdata_code/data/minji/")){
    if(endsWith(file,".RData")){
      load(paste("/Users/yaejilim/Library/CloudStorage/GoogleDrive-yaeji.lim@stat.cau.ac.kr/My Drive/Papers/MIN_JI_clustering/Stepdata_code/data/minji/",file,sep = ""))
      n<-length(foo$Steps)
      
      toMatrix=TRUE #check if each data ends at 11:59:00 pm
      
      for(l in 1:as.integer(n/1440)){
        if(!endsWith(toString(foo[1440*l,"ActivityHour"]),"59:00")){
          print(foo[1440*l,"ActivityHour"])
          toMatrix=FALSE
        }
      }
      
      if(toMatrix){
        A<-matrix(foo$Steps[1:(as.integer(n/1440)*1440)],nrow=1440)
        step[[idx]]<-A
      }
      idx = idx+1
    }
  }
  #delete days with all zero steps
  step.mat = NULL
  for(i in 1:length(step)){
    if(length(which(apply(step[[i]],2,sum)==0))){
      step[[i]]<-step[[i]][,-which(apply(step[[i]],2,sum)==0)]
    }
    step.mat=cbind(step.mat, matrix(step[[i]],nrow=1440))
  }
  
  return(list(step.mat = step.mat, step.list = step))
}
step0 = read_step()
step = step0$step.mat
step.list = step0$step.list
