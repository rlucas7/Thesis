###
###
###


# first run the "worse_tree_example.r" script, then after this we should run this script
min(dat3[,1])
max(dat3[,1])

min(dat3[,2])
max(dat3[,2])


### first region
postscript(file="proj_plot.eps",horizontal=FALSE, paper="special", width=8, height=7)

plot(c(0,1),c(0,1),xlab="", ylab="",col="lightgrey",las=1,cex.axis=1.2,mgp=c(0.1,0.5,0))
mtext("The True Tree",3,cex=1.75)
mtext("Covariate 2",2, line=3,cex=1.4)
mtext("Covariate 1",1,line=3,cex=1.4)

ind<-1:100
ind<-ind[dat3[1:100,6]==3]
points(dat3[ind,1],dat3[ind,2],pch=51)

ind<-1:100
ind<-ind[dat3[1:100,6]==0]
points(dat3[ind,1],dat3[ind,2],pch=48)

### second region (lower left)

ind<-101:120
ind<-ind[dat3[101:120,6]==2]
points(dat3[ind,1],dat3[ind,2],pch=50)

ind<-101:120
ind<-ind[dat3[101:120,6]==0]
points(dat3[ind,1],dat3[ind,2],pch=48)

### third region (upper left)

ind<-121:190
ind<-ind[dat3[121:190,6]==1]
points(dat3[ind,1],dat3[ind,2],pch=49)

ind<-121:190
ind<-ind[dat3[121:190,6]==0]
points(dat3[ind,1],dat3[ind,2],pch=48)


### fourth region (upper right)


ind<-191:250
ind<-ind[dat3[191:250,6]==2]
points(dat3[ind,1],dat3[ind,2],pch=50)

ind<-191:250
ind<-ind[dat3[191:250,6]==0]
points(dat3[ind,1],dat3[ind,2],pch=48)

abline(h=.5,col="red")

lines(c(0.55,0.55),c(0.5,1.2), col="blue")

lines(c(0.7,0.7),c(0,0.5), col="blue")

dev.off()