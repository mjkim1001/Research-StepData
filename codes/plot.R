
TT=500
data_temp= data_generation(n=n, TT=TT, model='1A')


par(mfrow=c(1,2))
plot(data_temp[[1]][1,], type='l', ylab='', xlab='', main='Model 1 (g=1)')
abline(v=c(54,128),lty=2, lwd=2 ,col=2)
plot(data_temp[[2]][1,], type='l', ylab='', xlab='', main='Model 1 (g=2)')
abline(v=c(54,128),lty=2, lwd=2 ,col=2)



TT=500
data_temp= data_generation(n=n, TT=TT, model='2A')


par(mfrow=c(2,2))
plot(data_temp[[1]][1,], type='l', ylab='', xlab='', main='Model 2(a) (g=1)')
plot(data_temp[[2]][1,], type='l', ylab='', xlab='', main='Model 2(a) (g=2)')

data_temp= data_generation(n=n, TT=TT, model='2B')
plot(data_temp[[1]][1,], type='l', ylab='', xlab='', main='Model 2(b) (g=1)')
plot(data_temp[[2]][1,], type='l', ylab='', xlab='', main='Model 2(b) (g=2)')

