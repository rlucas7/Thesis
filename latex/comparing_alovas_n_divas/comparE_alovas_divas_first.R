##### code to look at the relationship between interaction coefficient 
##### (sign and magintude) and how that impacts a resultant dtree.

set.seed(1)
library(tree)
p<-100
n<-500
non_zero <- 3
Y<- matrix(NA_real_, nrow=n, ncol=1)
X<- matrix(NA_real_, nrow=n, ncol=p)
true_beta <- c(rep(1 , non_zero), rep( 0 , p - non_zero))
interact <- 0.5
true_beta[p]<- interact
## simulate data 

for(i in 1:n){
	
	
	X[i,1:(p-1)] <-rnorm(p-1, mean = 1, sd =1)
	X[i,p] <- X[i,1]*X[i,2]
	
	Y[i] <- rnorm(1,mean = sum( X[i,]*true_beta  ), sd=1 )
	
	
}# end of for i in 1:n looop

lm_mod <- lm(Y~X)
tree_mod <- tree(Y~X)
getwd()
pdf(file="compare_divas_alovas.pdf")
par(fin=c(7,7))
plot(tree_mod,type='uniform')
text(tree_mod)
dev.off()
## now if the interaction is ommitted 

X<-X[,1:(p-1)]
tree_mod <- tree(Y~X)

pdf(file="compare_divas_alovas2.pdf")
plot(tree_mod,type='uniform')
text(tree_mod)
dev.off()

### load the stuff into wd needed 
source('/Users/rlucas7/ALoVaS/R_code/run_first_ALoVaS.r') # fur ALoVaS
### end of load the stuff 



####===========================SOC fit the divas model=========================== 
casper_data <- data.frame(Y,X)
	n=1000
	m=10
	data_dim=100
	
# SETTING UP DATA STRUCTURES....
t_lh<-matrix(NA,nrow=n*m,ncol=1) #matrix tracks the likelihood of the tree at the current iteration of the MCMC chain
output<-matrix(NA, nrow=n*m, ncol=2)
rtree<-vector(mode="list",length=n*m)

alpha_prime<-matrix(1,nrow=data_dim)
alpha_prime_store<-matrix(NA, nrow=n*m,ncol=length(alpha_prime))
weight<-matrix(NA, nrow=n*m,ncol=length(alpha_prime)-1)
count<-matrix(NA,nrow=n*m, ncol=length(alpha_prime))
#INITIALIZATIONS
f1 <- as.formula(paste(names(casper_data)[1],paste(names(casper_data)[-1],sep='', collapse='+'), collapse='~', sep='~') )

t1 <- tree(f1, data=casper_data)

t_lh[1]<-tree_lhood(t1,casper_data,theta=0,mult_penalty=FALSE)
output[1,1]<-length(unique(t1$frame[,1]))-1
output[1,2]<-length(get_tnodes(t1))
rtree[[1]]<-t1

Z<-matrix(NA, nrow=n*m, ncol=1) #here the "Z" vector represents the constraint multipliers $\widetilde{\alpha}$
t2<-t1 # the "optimal" tree found by the MCMC search   ////////// ????Is this line still necessary for the code ?????
Z[1]<-1
### initializations 

alpha_prime_store[1,]<-alpha_prime
weight[1,]<-1/sum(alpha_prime)

split_counts<-0

old_count<-0


for(j in 1:m){
	i<-1
	tree1<-t1
	tree2<-tree1 
	Z[i+1000*(j-1)]<-1
	alpha_prime_store[i+1000*(j-1) ,]<-alpha_prime
	weight[i+1000*(j-1),]<-1/length(alpha_prime)
 	t_lh[i+1000*(j-1)]<-tree_lhood(tree2,casper_data,theta=0,mult_penalty=FALSE)
 	split_counts<-0
for(i in 2:n){
	x<-rgen_discrete_unif()
			if( x == 1 ){ 
			# sample the grow
			tree2<-grow(tree1,casper_data,1,weight[1,])		
			tree2<-up(tree2,casper_data,1)
			}else{
				if(x ==2){
				tree2<-prune(tree1,casper_data,1)
				tree2<-up(tree2,casper_data,1)
			}else{ 
				if(x==3 ){
					tree2<-change(tree1,casper_data,1,weight[1000*(j-1)+i-1,])		
					tree2<-up(tree2,casper_data,1)
				}else{
					if(x==4 ){
						#rotate step...
						tree2<-rotate_k(tree1)
						tree2<-up(tree2,casper_data,1)	
						}else{ # x==5 case
						tree2<-swap(tree1,casper_data,1)
						tree2<-up(tree2,casper_data,1)
						#include swap function here, .....
						}
					}	
				}
			}
	
	t2_lh<-tree_lhood(tree2,casper_data,theta=0,mult_penalty=FALSE)	# calculate the proposed tree's l-hood	
	#proposed_tlh[i+1000*(j-1)]<-t2_lh
	### log ratio calculation 
	
	### old code:	
	#R<- t2_lh - t_lh[i+1000*(j-1)-1] 
	
	###NEW code: 
		### Start of likelihood ratio calculation 
	
	R<-log_ratio_calc(tree2,tree1,t2_lh, t_lh[i+1000*(j-1)-1],data=casper_data,weight=weight[i+1000*(j-1)-1,])
	
	#R<- t2_lh - c_t_lh[i+1000*(j-1)-1] 
	###End of likelihood ration calculation
	
	### end dof log ratio  calculation
	
	M_E<- -rexp(1,1)
	if(M_E <=R){ 
		t_lh[i+ (j-1)*1000 ]<- t2_lh # make the move to the proposed tree
			# stay at tree2 so do not update the tree2 object
			tree1<-tree2
			output[i+(j-1)*1000,1]<-length(unique(tree1$frame[,1]))-1
			output[i+(j-1)*1000,2]<-length(get_tnodes(tree1))
		}else{
			# stay at the current tree state 
			t_lh[i+(j-1)*1000]<-t_lh[1000*(j-1)+i-1]
			output[i+(j-1)*1000,1]<-length(unique(tree1$frame[,1]))-1
			output[i+(j-1)*1000,2]<-length(get_tnodes(tree1))			#num_term_nodes[i+(j-1)*1000-1]
		}
		rtree[[i+(j-1)*1000]]<-tree1
	if(i %%300 ==0 ){
		print(tree1)
		}
	
	#### now sample from the weights for each dimension.... 
	split_counts<-count_cov_split(tree=tree1,data=casper_data,1) + split_counts 
	Z[(j-1)*1000+i]<- data_dim/sum(split_counts)
	#split_counts<-(split_counts)*Z[(j-1)*1000+i]/10 # here the dividing number must change as the dimensionality of data does 
	weight[(j-1)*1000+i,]<-sample_weights(alpha=split_counts) ### the 1's are added inside the function.... 

	}# end of for(i 1:n) loop
} #end of for j in 1:m  loop

objects()
match(max(t_lh), t_lh)
str(rtree[8330], max.level=1)

rownames(tree1$frame )
rownames(rtree[[8330]]$frame) 
rtree[[8330]]
tree1
par(mfcol=c(1,1))
getwd()
pdf(file='greedy_tree_alovas_divas_compare.pdf')
plot(tree1, type='uniform')
text(tree1)
dev.off()

pdf(file='divas_tree_alovas_divas_compare.pdf')
plot(rtree[[8330]], type='uniform')
text(rtree[[8330]])
dev.off()



####===========================EOC to fit the divas model===========================



####===========================SOC fit the alovas model=========================== 


## make some data in the object format for the ALoVaS models

	simulation_results<-vector(mode='list',length=4)
	simulation_results[[1]]<-Y
	simulation_results[[2]]<-X
	simulation_results[[3]]<-true_beta
	simulation_results[[4]]<-NULL
	names(simulation_results)[1]<-'Y, the response vector'
	names(simulation_results)[2]<-'X, the design matrix'
	names(simulation_results)[3]<-'beta, the true coefficients'
	str(simulation_results, max.level=1)

## END OF MAKE SOME DATA IN THE OBJECT BLAH BLAH BLAH for the ALoVaS models



delta_grid<-seq(from=1,to=25, by=1)
MSE_grid<-seq(from=1,to=25, by=1)
n<-length(delta_grid)
n

results<-matrix(NA_real_, nrow=7, ncol=1)

ptm <- proc.time()

for(i in 2:n){
	dg <- delta_grid[i]
	ltree1<-lasso_tree(sim=simulation_results,n=1000,m=1,r_p =1, delta = dg, vocal=TRUE)
	MSE_grid[i] <- MSE_calc(ltree1,simulation_results)
}

mini <- match(min(MSE_grid),MSE_grid)
optimum <- delta_grid[mini]
# run one of these three

ltree1 <- lasso_tree(sim=obj,n=1000,m=1,r_p =1, delta = optimum, vocal=TRUE)
#stree1<-ssvs_tree(sim=obj,n=100,m=1,sig_sq_nz = rep(1,mz), sig_sq_z = rep(dg,mz) , vocal=TRUE)
#htree1 <- horseshoe_tree(sim=obj,n=1000,m=1,tau_sq = optimum, vocal=TRUE)

# write results to one of these three
 write.csv(results, file='lasso_compare_results.csv')
 #write.csv(results, file='lasso_compare_results.csv')
 #write.csv(results, file='lasso_compare_results.csv')

ind <- match( min(ltree1[[1]]), ltree1[[1]] )   
pdf(file='lasso_compare_tree2.pdf')
plot(ltree1[[2]][[ind]], type='uniform')
text(ltree1[[2]][[ind]])
 dev.off()
proc.time() - ptm
system('echo "this goes into email" | mail "done" rlucas7@vt.edu' )

getwd()
setwd("/Users/rlucas7/thesis/latex/comparing_alovas_n_divas")

pdf(file='elbow_plot_lasso_compare.pdf')

plot(delta_grid,MSE_grid/100,type='l')
dev.off()



####===========================EOC to fit the alovas model===========================


####===========================SOC to fit the cgm model===========================
casper_data <- data.frame(Y,X)
	n=1000
	m=10
	data_dim=101
	
# SETTING UP DATA STRUCTURES....
t_lh<-matrix(NA,nrow=n*m,ncol=1) #matrix tracks the likelihood of the tree at the current iteration of the MCMC chain
output<-matrix(NA, nrow=n*m, ncol=2)
rtree<-vector(mode="list",length=n*m)

#INITIALIZATIONS

t_lh[1]<-tree_lhood(t1,casper_data,theta=0,mult_penalty=FALSE)
output[1,1]<-length(unique(t1$frame[,1]))-1
output[1,2]<-length(get_tnodes(t1))
rtree[[1]]<-t1

t2<-t1 # the "optimal" tree found by the MCMC search   ////////// ????Is this line still necessary for the code ?????

for(j in 1:m){
	i<-1
	tree_num<-sample(1:hmm,1,prob=x_double_prime)
	tree1<-tree_init_mat[[tree_num]]
	tree2<-tree1 
	Z[i+1000*(j-1)]<-1
	alpha_prime_store[i+1000*(j-1) ,]<-alpha_prime
	weight[i+1000*(j-1),]<-1/length(alpha_prime)
 	t_lh[i+1000*(j-1)]<-tree_lhood(tree2,casper_data,theta=0,mult_penalty=FALSE)
 	split_counts<-0
for(i in 2:n){
	x<-rgen_discrete_unif()
			if( x == 1 ){ 
			# sample the grow
			tree2<-grow(tree1,casper_data,1,weight[1,])		
			tree2<-up(tree2,casper_data,1)
			}else{
				if(x ==2){
				tree2<-prune(tree1,casper_data,1)
				tree2<-up(tree2,casper_data,1)
			}else{ 
				if(x==3 ){
					tree2<-change(tree1,casper_data,1,weight[1000*(j-1)+i-1,])		
					tree2<-up(tree2,casper_data,1)
				}else{
					if(x==4 ){
						#rotate step...
						tree2<-rotate_k(tree1)
						tree2<-up(tree2,casper_data,1)	
						}else{ # x==5 case
						tree2<-swap(tree1,casper_data,1)
						tree2<-up(tree2,casper_data,1)
						#include swap function here, .....
						}
					}	
				}
			}
	
	t2_lh<-tree_lhood(tree2,casper_data,theta=0,mult_penalty=FALSE)	# calculate the proposed tree's l-hood	
	#proposed_tlh[i+1000*(j-1)]<-t2_lh
	### log ratio calculation 
	
	### old code:	
	#R<- t2_lh - t_lh[i+1000*(j-1)-1] 
	
	###NEW code: 
		### Start of likelihood ratio calculation 
	
	R<-log_ratio_calc(tree2,tree1,t2_lh, t_lh[i+1000*(j-1)-1],data=casper_data,weight=weight[i+1000*(j-1)-1,])
	
	#R<- t2_lh - c_t_lh[i+1000*(j-1)-1] 
	###End of likelihood ration calculation
	
	### end dof log ratio  calculation
	
	M_E<- -rexp(1,1)
	if(M_E <=R){ 
		t_lh[i+ (j-1)*1000 ]<- t2_lh # make the move to the proposed tree
			# stay at tree2 so do not update the tree2 object
			tree1<-tree2
			output[i+(j-1)*1000,1]<-length(unique(tree1$frame[,1]))-1
			output[i+(j-1)*1000,2]<-length(get_tnodes(tree1))
		}else{
			# stay at the current tree state 
			t_lh[i+(j-1)*1000]<-t_lh[1000*(j-1)+i-1]
			output[i+(j-1)*1000,1]<-length(unique(tree1$frame[,1]))-1
			output[i+(j-1)*1000,2]<-length(get_tnodes(tree1))			#num_term_nodes[i+(j-1)*1000-1]
		}
		rtree[[i+(j-1)*1000]]<-tree1
	if(i %%300 ==0 ){
		print(tree1)
		}
	}# end of for(i 1:n) loop
} #end of for j in 1:m  loop


####===========================EOC to fit the cgm model===========================



####===========================SOC to do the hold out error calc=========================== 

## first simulate hold out data
Y2<- matrix(NA_real_, nrow=n, ncol=1)
X2<- matrix(NA_real_, nrow=n, ncol=p)
for(i in 1:n){
	
	
	X[i,1:(p-1)] <-rnorm(p-1, mean = 1, sd =1)
	X[i,p] <- X[i,1]*X[i,2]
	
	Y[i] <- rnorm(1,mean = sum( X[i,]*true_beta  ), sd=1 )
	
	
}# end of for i in 1:n looop



### now write the code to calculate the misclass/hould-out error 

# divas 

# alovas 

# cgm

# vimp (ishwaran)

 



####===========================EOC to do the hold out error calc=========================== 