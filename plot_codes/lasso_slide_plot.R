### code to make the plot of exponential, laplace and normal densities

x<-seq(-5,5,by=0.01)
nx<-dnorm(x,0,1)
ex<-dexp(x,0.4)
require(VGAM)
lx<-dlaplace(x,0,1)

#postscript(file="lasso.eps",horizontal=FALSE, paper="special", width=8, height=7)
pdf(file="lasso.pdf")
plot(x,lx,xlab="x",ylab="density", main="plots of Normal, Exponential, and Laplace",type="l",col="red")
lines(x,nx,col="blue")
lines(x,ex)
legend(-5,0.4,c("Laplace", "Normal","Exponential"), lty=c(1,1,1),col = c("red","blue","black"))
dev.off()