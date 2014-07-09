### code to produce the shattering dynamics plot. 


n<-seq(0,20, by=0.01)
fn<-2^n
for(i in 1001:length(x)){
	fn[i]<-2^10+50*n[i]
}
getwd()
pdf(file="phase_change_shatter_dynamics.pdf")
plot(n,fn,type="l", main="Shatter Coefficient growth with n")
dev.off()