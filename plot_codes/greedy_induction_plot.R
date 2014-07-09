#### code to make the plot illustrating optimizing the greedy objective
### when inducting a decision tree
set.seed(1)
x<-rnorm(10,0,1)
y<-rnorm(10,0,1)


postscript(file="greedy_induction.eps",horizontal=FALSE, paper="special", width=8, height=7)
#pdf(file="greedy_induction.pdf")
plot(x,y,main="Illustrating greedy induction")
abline(v=1.1)
abline(v=.66, lty=2)
abline(v=.53,lty=3)
abline(v=.42,lty=4)
legend(x=-.5,y=-1,c("1st", "2nd","3rd","4th"), lty=c(1,2,3,4))
dev.off()