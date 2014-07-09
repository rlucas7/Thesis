### code to make the moment prior plot  for writeup 
###
###	
###		Date: Tuesday March 26th, 2013		
###--------------------------------------------------------------------

# SECOND MOMENT PRIOR OF VALEN JOHNSON'S 2010 JRSS PAPER (G(X))
# RATIONAL FUNCTION ALTERNATIVE (MY SUGGESTION (HERE DENOTED f(X) ) )


x<-seq(from=-50,to=50, by=0.01)
fx<-  0.1*x^2/(0.009*x^2+0.27)^2
gx<-((x)^2)*exp(-0.5*(.09)*(x)^2)
gx[5001]<-4

hx<-0.1*x^2/(0.009*x^2+0.27)^3
hx<-hx/2
getwd()
pdf(file="moment_prior_plot_BW.pdf")
plot(x,fx,type="l", ylab="density",main="Moment prior density functions")
lines(x,gx,type="l", lty=46)
legend(-50,8,c( "t-moment(d.f.=2)","normal-moment"), lty = c(1,46))
dev.off()

pdf(file="moment_prior_plot_color.pdf")
plot(x,hx,type="l", ylab="density",main="Moment prior density functions")
lines(x,gx,type="l", col="red")
lines(x,fx,type="l",col="blue")
legend(-50,8,c( "t-moment(d.f.=5)","normal-moment","t-moment(d.f.=2)"), lty = c(1,1),col=c("black","red","blue"))
dev.off()


postscript(file="moment_prior_plot_color.eps",horizontal=FALSE, paper="special", width=8, height=7)
plot(x,hx,type="l", ylab="density",main="Moment prior density functions")
lines(x,gx,type="l", col="red")
lines(x,fx,type="l",col="blue")
legend(-50,8,c( "t-moment(d.f.=5)","normal-moment","t-moment(d.f.=2)"), lty = c(1,1),col=c("black","red","blue"))
dev.off()