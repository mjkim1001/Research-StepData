#################################################################
################################################################
################## Read Data and          ######################
################## Data Preprocessing     ######################
################################################################
################################################################
# step.list returns a list of step data for each individual.
# step.mat returns concateneated step data matrix from all participants.
library(lubridate)
read_step<-function(){
  step<-list()
  is.wkend = vector()
  step.mat = NULL
  idx=1
  for(file in list.files(path="./data2/minji/")){
    if(endsWith(file,".RData")){
      load(paste("./data2/minji/",file,sep = ""))
      n<-length(foo$Steps)
      
      toMatrix=TRUE #check if each data ends at 11:59:00 pm
      
      for(l in 1:as.integer(n/1440)){
        if(!endsWith(toString(foo[1440*l,"ActivityHour"]),"59:00")){
          print(foo[1440*l,"ActivityHour"])
          toMatrix=FALSE
        }
      }
      
      if(toMatrix){
        lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
        A<-matrix(foo$Steps[1:(as.integer(n/1440)*1440)],nrow=1440)
        #return if each date is weekends (Sun or Sat)
        wkend = !(wday(as.Date(substring(matrix(foo$ActivityHour[1:(as.integer(n/1440)*1440)],nrow=1440)[1,],1,7),format = "%d%h%y" )) %in% 2:6 )
        #delete days with all zero steps
        if(length(which(apply(A,2,sum)==0))){
          is.wkend = c(is.wkend, wkend[-which(apply(A,2,sum)==0)])
          A<-A[,-which(apply(A,2,sum)==0)]
        }else{
          is.wkend = c(is.wkend, wkend)
        }
        step[[idx]]<-A
        step.mat=cbind(step.mat, A)
      }
      
      idx = idx+1
    }
  }

  return(list(step.mat = step.mat, step.list = step, is.wkend = is.wkend))
}
step0 = read_step()
step = step0$step.mat
step.list = step0$step.list
is.wkend = step0$is.wkend
