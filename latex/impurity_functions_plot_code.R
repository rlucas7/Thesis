### code for impurity function plots

x <- seq(0.001,0.999,by=0.001)

f1<- -x*log(x) -(1-x)*log(1-x)
f2<- 2*x*(1-x)
f3<- pmin(x,1-x)

pdf(file='impurity_plot.pdf')
plot(x,f1,type='l')
lines(x,f2, lty=2)
lines(x,f3, lty=3)
dev.off()