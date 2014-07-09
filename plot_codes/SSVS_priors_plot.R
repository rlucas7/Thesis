#### code to make the SSVS mixture priors 


x<-seq(-5,5,by=0.01)
nx1<-dnorm(x,mean=0,sd=1)
nx2<-dnorm(x,mean=0,sd=3)
nx3<-dnorm(x,mean=0,sd=10)

#postscript(file="SSVS_prior_plot.eps",horizontal=FALSE, paper="special", width=8, height=7)
#pdf(file="SSVS_prior_plot.pdf")
# plot(x,nx1,main="SSVS prior",xlab="x", ylab="normal density",type="l")
# lines(x,nx2,lty=2)
# lines(x,nx3,lty=3)
# legend(x=-4.5,y=.3,c("c=100","c=300","c=1000"), lty=c(1,2,3))
# dev.off()


postscript(file="SSVS_prior_plot_color.eps",horizontal=FALSE, paper="special", width=8, height=7)
#pdf(file="SSVS_prior_plot_color.pdf")
plot(x,nx1,main="SSVS prior",xlab="x", ylab="normal density",type="l")
lines(x,nx2,col="red")
lines(x,nx3,col="blue")
legend(x=-4.5,y=.3,c("c=100","c=300","c=1000"), lty=c(1,1,1),col=c("black","red","blue"))
dev.off()