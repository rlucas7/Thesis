#### code to illustrate the decision tree splitting graphic in ".eps" 
#### format
set.seed(1)
x<-rnorm(10, 0, 1)
y<-rnorm(10, 0, 1)


#postscript(file="partition_plot_color.eps",horizontal=FALSE, paper="special", width=8, height=7)
pdf(file="partition_plot_color.pdf")
plot(x,y,ylab="response", xlab="covariate/partition rule")
abline(v=1.2)
abline(v=0.65,col="red")
abline(v=0.53,col="blue")
abline(v=0.41,col="purple")
dev.off()
