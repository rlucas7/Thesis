#### code to make the half-cauchy plots 

require(LaplacesDemon)

x<-seq(-5,5,by=0.01)
nx1<-dnorm(x,mean=0,sd=1)
nx2<-dhalfcauchy(x,scale=1)
nx3<-dinvgamma(x,shape=1,scale=1)

postscript(file="half_cauchy_plot.eps",horizontal=FALSE, paper="special", width=8, height=7)
#pdf(file="half_cauchy_plot.pdf")
plot(x,nx1,main="Half-Cauchy prior",xlab="x", ylab="normal density",type="l")
lines(x,nx2,lty=2)
lines(x,nx3,lty=3)
egend(x=-4.5,y=.3,c("Half-Cauchy(1)","Inv-Gamma(1,1)","Normal(0,1)"), lty=c(1,1,1),col=c("black","red","blue"))dev.off()

nx2<-nx2*2
nx2[x<0]<-0
#postscript(file="half_cauchy_plot_color.eps",horizontal=FALSE, paper="special", width=8, height=7)
pdf(file="half_cauchy_plot_color.pdf")
plot(x,nx2,,xlab="x", ylab="density",type="l")
lines(x,nx3,col="red")
lines(x,nx1,col="blue")
legend(x=-4.5,y=.8,c("Half-Cauchy(1)","Inv-Gamma(1,1)","Normal(0,1)"), lty=c(1,1,1),col=c("black","red","blue"))
dev.off()