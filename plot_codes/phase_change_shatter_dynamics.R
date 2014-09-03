### code to produce the shattering dynamics plot. 


n<-seq(0,20, by=0.01)
fn<-2^n
for(i in 1001:length(n)){
	fn[i]<-2^10+50*n[i]
}
getwd()
setwd("/Users/rlucas7/Desktop")
pdf(file="phase_change_shatter_dynamics3.pdf")
plot(n,fn,type="l", ylab="Shatter Coefficient growth with n")
segments(0,1535,10,1525,lty=2)
segments(10.05,0,10.05,1525,lty=2)
#arrows(15,1000,10.2,10)
#text(15,1030,'VC dimension')
dev.off()