# code creates data from CGM simulated example from the CGM 1998 paper. 


#### requires the tree helper functions are loaded into the wd... 
set.seed(17)
c<-5		# max depth for trees .... 
nobs<-600	# the number of observations for the simulated data
nc<-100# the total number of covariates used in the 
##### CGM 1998 data creation .... 

x1<-matrix(NA, nrow=nobs, ncol=1)
x2<-matrix(NA, nrow=nobs, ncol=1)

for(i in 1:nobs){
x2[i]<-rgen_discrete_unif_four()
x1[i]<- runif(1,min=0, max=10)
}



y<-matrix(NA, nrow=nobs, ncol=1)

for(i in 1:nobs){
if(x1[i] <=5.0 && x2[i] ==1 ){ y[i]<-8}
if(x1[i] <=5.0 && x2[i] ==2 ){ y[i]<-8}

if(x1[i] >5.0 && x2[i] ==1 ){ y[i]<-2}
if(x1[i] >5.0 && x2[i] ==2 ){ y[i]<-2}

if(x1[i] <=3.0 && x2[i] ==3 ){ y[i]<-1}
if(x1[i] <=3.0 && x2[i] ==4 ){ y[i]<-1}

if(x1[i] <=3.0 && x2[i] ==3 ){ y[i]<-8}
if(x1[i] <=3.0 && x2[i] ==3 ){ y[i]<-8}

if(x1[i] >3.0 && x1[i] <=7 && x2[i] == 3 ){ y[i]<-5}
if(x1[i] >3.0 && x1[i] <=7 && x2[i] == 4 ){ y[i]<-5}

if( x1[i] > 7 && x2[i] == 3 ){ y[i]<-8}
if( x1[i] > 7 && x2[i] == 4 ){ y[i]<-8}

y[i]<-rnorm(1, mean=0, sd=2)
}


extra_data<-matrix(NA, nrow=nobs, ncol=nc+1)
for(i in 1:(nc-2)){
two_means<-runif(2,0,30)	
	for(j in 1:nobs){
	#x<-1+rbinom(1,1,0.5)	
			extra_data[j,i]<-rnorm(1,mean=0,sd=.5)
	}
}


head(extra_data)

### three columns left at the end to add x1, x2, and y. 


extra_data[,nc+1]<-y
extra_data[,nc-1]<-x1
extra_data[,nc]<-as.factor(x2)
extra_data<-as.data.frame(extra_data)
#extra_data$V100<-as.factor(x2)


names(extra_data)<-c(paste("x",1:nc,sep=""),"y")
head(extra_data)

## end of extra dimensions code 

f1<-as.formula(paste(names(extra_data)[101],paste(names(extra_data)[-101],collapse="+"), sep="~"))
f1

library(tree)
t1<-tree(f1,data=extra_data)
pdf(file="greedy_100_covar_tree.pdf")
plot(t1, type="unif")
text(t1)
dev.off()
getwd()
names(extra_data)

#t1<-prune(t1,data=extra_data,ycol=101)
#t1