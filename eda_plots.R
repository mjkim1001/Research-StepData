library(forecast)
library(colorspace)
library(fields)
library(progress)

#################################################################
################################################################
##################        Figures        ######################
################################################################
################################################################

source("./functions.R")

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


selected = c(16408,18595,556) # days selected for figures in the thesis

time = seq(as.POSIXct("00:00:00",format="%H:%M"),
           as.POSIXct("23:59:00",format="%H:%M"), by = "1 min")

#Fig2
par(mfrow=c(1,3))
for(i in c(1,200,556)){
  plot(time,step[,i],xlab="",ylab="", ylim=c(-5,150),type='b',pch=20, cex=0.7)
}  
  
# par(mfrow=c(1,3))
# par(mar=c(3.5,3,1.5,1))
# par(oma=c(3,2,1,1.5))
# par(mgp=c(2.5,1,0))
# par(mar=c(2,2,1,1))
# par(oma=c(1.5,1.5,1,1.5))
# par(mgp=c(1,1,0))
# for(i in selected){
#   plot(time,step[,i],type='l',xlab="",ylab="", ylim=c(0,120))
# }

#EDA plots
par(mfrow=c(2,3))
for(i in selected){
  plot(time,step[,i],type='l',xlab="",ylab="", ylim=c(-20,140))
  polygon( c(time,rev(time)) ,c(cirBound(step[,i], 30, max), cirBound(step[,i], 30, min)),col="skyblue" )
  lines(time,step[,i],type='l')
}
for(i in selected){
  plot(time,step[,i],type='l',xlab="",ylab="", ylim=c(-50,170))
  polygon( c(time,rev(time)) ,c(cirBound(step[,i], 100, max), cirBound(step[,i], 30, min)),col="skyblue" )
  lines(time,step[,i],type='l')
}


par(mfrow = c(1,1))
plot(time,step[,selected[3]],type='l',xlab="(a)",ylab="", ylim=c(-50,150))

#Fig4
par(mfrow=c(2,2),mar = c(2, 2, 4, 2))
thick_draw(step[,selected[3]], 30, sqBound, time=TRUE, smooth = 0, smooth.num = 0,zeroBound = FALSE,xl = "(a)",ymin=-50,ymax=180)
thick_draw(step[,selected[3]], 80, sqBound, time=TRUE, smooth = 0, smooth.num = 5,zeroBound = FALSE,xl = "(b)",ymin=-50,ymax=180)

thick_draw(step[,selected[3]], 30, sqBound_es, time=TRUE, smooth = 0, smooth.num = 5,zeroBound = FALSE,xl = "(c)",ymin=-50,ymax=180)
thick_draw(step[,selected[3]], 80, sqBound_es, time=TRUE, smooth = 0, smooth.num = 5,zeroBound = FALSE,xl = "(d)",ymin=-50,ymax=180)


#Fig 5

selected = c(16408,159,526, 14527) # days selected for figures in the thesis

id=4
par(mfrow=c(3,2),mar = c(2, 2, 4, 2))
thick_draw(step[,selected[2]], 30, sqBound_es, time=TRUE, smooth = 0, smooth.num = 5,zeroBound = T,xl = expression("(a) e-TPT"),ymin=-20,ymax=160,col=rgb(0,0,1,0.3))
thick_draw(step[,selected[id]], 30, sqBound_es, time=TRUE, smooth = 0, smooth.num = 5,zeroBound = T,xl = expression("(b) e-TPT"),ymin=-20,ymax=160,col=rgb(1,0,0,0.2))

plot(time,step[,selected[1]],type='n',main=expression("(c) overlap"),ylab="", ylim=c(-20,160), xlab='')
polygon( c(time,rev(time)) ,c(sqBound_es(step[,selected[2]], 30, max), rep(0,1440)),col=rgb(0,0,1,0.3) )
polygon( c(time,rev(time)) ,c(sqBound_es(step[,selected[id]], 30, max), rep(0,1440)),col=rgb(1,0,0,0.2) )
#lines(time,step[,selected[2]],type='l',col='blue', lty=2)
#lines(time,step[,selected[3]],type='l',col='red', lty=3)

plot( time, pmin(sqBound_es((step[,selected[2]]), 30, max),
                 sqBound_es((step[,selected[id]]), 30, max))/pmax(sqBound_es((step[,selected[2]]), 30, max),sqBound_es((step[,selected[id]]), 30, max)),      type='l' ,main=expression("(d) TPMA"[0] *" (e-TPT)"), ylab="" , ylim=c(0,1), xlab='')



#	TPMA0 based on TPT and TPMA based on e-TPT. 
# par(mfrow=c(1,3),mar = c(2, 2, 4, 2))
	# TPMA0_etpt =	pmin(sqBound_es((step[,selected[2]]), 30, max),sqBound_es((step[,selected[id]]), 30, max))/pmax(sqBound_es((step[,selected[2]]), 30, max),sqBound_es((step[,selected[id]]), 30, max))
# plot( time, TPMA0_etpt,   type='l' ,main= expression("(a) TPMA"[0] *" (e-TPT)"), ylab="" , ylim=c(0,1), xlab='')
	TPMA0_tpt =	pmin(sqBound((step[,selected[2]]), 30, max),sqBound((step[,selected[id]]), 30, max))/pmax(sqBound((step[,selected[2]]), 30, max), sqBound((step[,selected[id]]), 30, max))
plot( time, TPMA0_tpt,      type='l' ,main=expression("(e) TPMA"[0] *" (TPT)"), ylab="" , ylim=c(0,1), xlab='')

TPMA_etpt =	(pmin(sqBound_es((step[,selected[2]]), 30, max), sqBound_es((step[,selected[id]]), 30, max))- pmax(sqBound_es((step[,selected[2]]), 30, min), sqBound_es((step[,selected[id]]), 30, min)) )/(pmax(sqBound_es((step[,selected[2]]), 30, max), sqBound_es((step[,selected[id]]), 30, max))- pmin(sqBound_es((step[,selected[2]]), 30, min), sqBound_es((step[,selected[id]]), 30, min)) )
plot( time, TPMA_etpt,      type='l' ,main=expression("(f) TPMA"*" (e-TPT)"), ylab="" , ylim=c(0,1), xlab='')




