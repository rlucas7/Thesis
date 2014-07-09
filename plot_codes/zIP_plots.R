#### code to make the half-cauchy plots 

require(VGAM)

x<-seq(0,20,by=1)
nx1<-dpois(x,lambda=3)
nx2<-dnbinom(x,size=5, prob=0.5)
nx3<-dzipois(x,lambda=3,pstr0=0.2)

#postscript(file="ZIP_plot.eps",horizontal=FALSE, paper="special", width=8, height=7)
#pdf(file="ZIP_plot.pdf")
plot(x,nx3,xlab="x", ylab="density",type="l", xlim=c(0,20))
axis(1,at=seq(0,20,by=1))
axis(1,at=seq(12,20,by=1))
axis(1,at=seq(11,19,by=1))
lines(x,nx2,lty=2)
lines(x,nx1,lty=3)
legend(x=13,y=.15,c("ZIP(3,0.2)","NegBinom(5,1/2)","Poisson(3)"), lty=c(1,2,3))
dev.off()


postscript(file="ZIP_plot_color.eps",horizontal=FALSE, paper="special", width=8, height=7)
#pdf(file="ZIP_plot_color.pdf")
plot(x,nx3,xlab="x", ylab="density",type="l")
axis(1,at=seq(0,20,by=1))
axis(1,at=seq(12,20,by=1))
axis(1,at=seq(11,19,by=1))
lines(x,nx2,col="red")
lines(x,nx1,col="blue")
legend(x=13,y=.15,c("ZIP(3,0.2)","NegBinom(5,1/2)","Poisson(3)"), lty=c(1,1,1),col=c("black","red","blue"))
dev.off()