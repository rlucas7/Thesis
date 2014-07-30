##### code to look at the relationship between interaction coefficient 
##### (sign and magintude) and how that impacts a resultant dtree.

set.seed(1)
library(tree)
p<-100
n<-500
non_zero <- 10
Y<- matrix(NA_real_, nrow=n, ncol=1)
X<- matrix(NA_real_, nrow=n, ncol=p)
true_beta <- c(rep(1 , non_zero), rep( 0 , p - non_zero))
interact <- 4
true_beta[p]<- interact
## simulate data 

for(i in 1:n){
	
	
	X[i,1:(p-1)] <-rnorm(p-1, mean = 1, sd =1)
	X[i,p] <- X[i,1]*X[i,2]
	
	Y[i] <- rnorm(1,mean = sum( X[i,]*true_beta  ), sd=1 )
	
	
}# end of for i in 1:n looop

lm_mod <- lm(Y~X)
X<-X[,1:(p-1)]
tree_mod <- tree(Y~X)

plot(tree_mod,type='uniform')
text(tree_mod)

## now if the interaction is ommitted 

X<-X[,1:(p-1)]
tree_mod <- tree(Y~X)

plot(tree_mod,type='uniform')
text(tree_mod)


####===========================SOC fit the divas model=========================== 



####===========================EOC to fit the divas model===========================



####===========================SOC fit the alovas model=========================== 



####===========================EOC to fit the alovas model===========================


####===========================SOC to do the hold out error calc=========================== 

## first simulate hold out data
Y2<- matrix(NA_real_, nrow=n, ncol=1)
X2<- matrix(NA_real_, nrow=n, ncol=p)
for(i in 1:n){
	
	
	X[i,1:(p-1)] <-rnorm(p-1, mean = 1, sd =1)
	X[i,p] <- X[i,1]*X[i,2]
	
	Y[i] <- rnorm(1,mean = sum( X[i,]*true_beta  ), sd=1 )
	
	
}# end of for i in 1:n looop

lm_mod2 <- lm(Y~X)
X<-X[,1:(p-1)]
tree_mod2 <- tree(Y~X)

plot(tree_mod2,type='uniform')
text(tree_mod2)

### now write the code to calculate the misclass/hould-out error 

# divas 

# alovas 

# cgm

# vimp (ishwaran)

 



####===========================EOC to do the hold out error calc=========================== 