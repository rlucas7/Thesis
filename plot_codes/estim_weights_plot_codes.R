### code to plot the weights for the ads dataset analysis

getwd()
setwd("/Users/rlucas7/Desktop")
list.files()
estim_weights<-read.csv("estim_2_table_of_weights.csv")
estim_weights<-estim_weights[,-1]
estim_weights<-as.matrix(estim_weights)



postscript(file="weights_estim2.eps",horizontal=FALSE, paper="special", width=8, height=7)
plot(estim_weights[1,],las=1,cex.axis=1.3,mgp=c(0.1,0.5,0),xlab="",ylab="",pch=16,labels=FALSE)
axis(2,at=c(0.00,.05,.10,.15,.20,.25,.30,.35,.40),las=1,cex.axis=1.3)
axis(1,at=c(1,250,500,750,1000,1250,1500),las=1,cex.axis=1.3)
mtext("Covariates",1,line=3,cex=1.4)
mtext("Estimated Probabilities",3,line=1.5,cex=1.4)
dev.off()