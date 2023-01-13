library(dtw)
#################################################################
################################################################
################## DTW part               ######################
################################################################
################################################################
########data############
data = read.table(file='./data/synthetic/synthetic_control.data')
data = t(data)

#figure 3.2
par(mar=c(3.5,3.5,1,1))#space around figure!!
par(oma=c(3,2,1,1.5)) #outer margin
par(mgp=c(2.5,1,0)) #축과의 거리
par(mfrow = c(2,3))
data.sel = cbind(data[,11:15], data[,101:105],data[,216] ,data[,207:210],
                                    data[,311:315], data[,411:415], data[,c(501,521:523,525)])
for(i in 1:6){
  plot( data.sel[,5*(i-1)+2] , xlab="", ylab="", type='l' )
}

########Table 3.1 results ##########
par(mfrow = c(1,1))
#1. Euclidean hierarchical clustering
l1.simul = hclust(dist(t(data.sel),method = function(x,y){mean((x-y)^2)}),method = "average")
plot(l1.simul)
rect.hclust(l1.simul,k=6)
cl.simul = cutree(l1.simul,k=6)
table(cl.simul)

#2. DTW + euclidean hierarchical clustering
dtw.simul = hclust(dist(t(data.sel),method = "dtw"),method = "average")
plot(dtw.simul)
rect.hclust(dtw.simul,k=6)
cl.simul = cutree(dtw.simul,k=6)
table(cl.simul)

#3. DTW + TPMA hierarchical clustering
library(dtwclust)
tp.dist = dist(t(data.sel),method = function(x,y){
  dtw(-tp_mat( smooth.ma(x,5) ,smooth.ma(y,5),sqBound_es,4 ), window.type = "sakoechiba",window.size=10)$distance
})
tp.simul = hclust( tp.dist ,method = "average")
plot(tp.simul,xlab = "",ylab="") ####figure 3.3
rect.hclust(tp.simul,k=6)
cl.simul = cutree(tp.simul,k=6)
table(cl.simul)
#function for TMPA distance matrix
tp_mat <- function(x,y, FUN,tau,r=1){
  ux = FUN(x,tau,max,r)
  uy = FUN(y,tau,max,r)
  lx = FUN(x,tau,min,r)
  ly = FUN(y,tau,min,r)
  ret.mat = matrix(NA, nrow = length(ux), ncol = length(uy))
  for(i in 1:length(ux)){
    for(j in 1:length(uy)){
      ret.mat[i,j] = {min(ux[i] ,uy[j] )-max(lx[i] ,ly[j] )}/{max(ux[i],uy[j])-min(lx[i],ly[j])}
    }
  }
  return(ret.mat)
}

########figure 3.4 result ##########
#dtw figure
alignment <- dtw(-tp_mat( smooth.ma(data.sel[,24],5) ,smooth.ma(data.sel[,25],5),sqBound_es,4 ),
                 window.type = "sakoechiba",window.size=10, keep=TRUE);
plot(smooth.ma(data.sel[,24],5) ,type='l',xlab="",ylab="")
lines(smooth.ma(data.sel[,25],5),type='l',col='blue')
for (i in 1:length(alignment$index1s)) {
  segments(alignment$index1s[i], smooth.ma(data.sel[,24],5)[alignment$index1s[i]],  alignment$index2s[i], (smooth.ma(data.sel[,25],5)[alignment$index2s[i]]), lty=4, col = 'darkgrey')
}

