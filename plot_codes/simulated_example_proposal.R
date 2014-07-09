###
###			Code that implements the constraint suggested by Scotland. 
###
###
###
###
###		Date: Tuesday April 25 2013  Date Modified: 
###-----------------------------------------------------------------------------------------------------------------------------------------------

require(MCMCpack)

### first run the code in file /tree_functions/worse_tree_example.r then run this script...

f1<-as.formula(paste(c(names(dat3)[6],paste(names(dat3)[-6],collapse="+")),collapse="~"))
f2<--as.formula(paste(c(names(dat3)[6],paste(names(dat3)[1:2],collapse="+")),collapse="~"))

require(tree)
t1<-tree(formula=f1,data=dat3)
t2<-tree(formula=f2,data=dat3)
par(mfcol=c(1,2))

plot(t1,type="uniform")
text(t1)

plot(t2,type="uniform")
text(t2)


### Creating AND Initializing the data structures to store the output from the MCMC chain... 

### function to do initialization of data for MCMC sampler 

	n=1000
	m=10
	data_dim=99
	
# SETTING UP DATA STRUCTURES....
t_lh<-matrix(NA,nrow=n*m,ncol=1) #matrix tracks the likelihood of the tree at the current iteration of the MCMC chain
output<-matrix(NA, nrow=n*m, ncol=2)
rtree<-vector(mode="list",length=n*m)

alpha_prime<-matrix(1,nrow=data_dim)
alpha_prime_store<-matrix(NA, nrow=n*m,ncol=length(alpha_prime))
weight<-matrix(NA, nrow=n*m,ncol=length(alpha_prime))
count<-matrix(NA,nrow=n*m, ncol=length(alpha_prime))
#INITIALIZATIONS

t_lh[1]<-tree_lhood(t1,dat3,theta=0,mult_penalty=FALSE)
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
	tree1<-t2
	tree2<-tree1 
	Z[i+1000*(j-1)]<-1
	alpha_prime_store[i+1000*(j-1) ,]<-alpha_prime
	weight[i+1000*(j-1),]<-1/length(alpha_prime)
 	t_lh[i+1000*(j-1)]<-tree_lhood(tree2,dat3,theta=0,mult_penalty=FALSE)
 	split_counts<-0
for(i in 2:n){
	x<-rgen_discrete_unif()
			if( x == 1 ){ 
			# sample the grow
			tree2<-grow(tree1,dat3,6,weight[1,])		
			tree2<-up(tree2,dat3,6)
			}else{
				if(x ==2){
				tree2<-prune(tree1,dat3,6)
				tree2<-up(tree2,dat3,6)
			}else{ 
				if(x==3 ){
					tree2<-change(tree1,dat3,1,weight[1000*(j-1)+i-1,])		
					tree2<-up(tree2,dat3,6)
				}else{
					if(x==4 ){
						#rotate step...
						tree2<-rotate_k(tree1)
						tree2<-up(tree2,dat3,6)	
						}else{ # x==5 case
						tree2<-swap(tree1,dat3,6)
						tree2<-up(tree2,dat3,6)
						#include swap function here, .....
						}
					}	
				}
			}
	
	t2_lh<-tree_lhood(tree2,dat3,theta=0,mult_penalty=FALSE)	# calculate the proposed tree's l-hood	
	#proposed_tlh[i+1000*(j-1)]<-t2_lh
	### log ratio calculation 
	
	### old code:	
	#R<- t2_lh - t_lh[i+1000*(j-1)-1] 
	
	###NEW code: 
		### Start of likelihood ratio calculation 
	
	R<-log_ratio_calc(tree2,tree1,t2_lh, t_lh[i+1000*(j-1)-1],data=dat3,weight=weight[i+1000*(j-1)-1,])
	
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
	split_counts<-count_cov_split(tree=tree1,data=dat3,6) + split_counts 
	Z[(j-1)*1000+i]<- data_dim/sum(split_counts)
	#split_counts<-(split_counts)*Z[(j-1)*1000+i]/10 # here the dividing number must change as the dimensionality of data does 
	weight[(j-1)*1000+i,]<-sample_weights(alpha=split_counts) ### the 1's are added inside the function.... 

	}# end of for(i 1:n) loop
} #end of for j in 1:m  loop

# end of outer loop


system("say yo! Ive got your  output")

plot(t_lh,type="l")

postscript(file="weight_simulated_example_USE1.eps",horizontal=FALSE, paper="special", width=8, height=7)
boxplot.matrix(weight[9600:10000,],use.cols=TRUE,las=1,cex.axis=1.2,mgp=c(0.1,0.5,0))
mtext("Distribution of Weights Across Dimensions",3,cex=1.75)
mtext("Estimated probabilities",2, line=3,cex=1.4)
mtext("Covariates",1,line=3,cex=1.4)
arrows(50,.2,10, .28)
text(60, .18,labels="the correct covariates",cex=2)
dev.off()

postscript(file="my_method1.eps",horizontal=FALSE, paper="special", width=8, height=7)
boxplot.matrix(weight[6600:7000,],use.cols=TRUE,las=1,cex.axis=1.2,mgp=c(0.1,0.5,0))
mtext("Distribution of Weights Across Dimensions",3,cex=1.75)
mtext("Estimated probabilities",2, line=3,cex=1.4)
mtext("Covariates",1,line=3,cex=1.4)
arrows(50,.2,10, .28)
text(60, .18,labels="the correct covariates",cex=2)
dev.off()

# so the sampler stays at the max when stared there..., and doesn't appear to pick up the other tree when we allow for swithces in dimension probabilities.... 

###THIS SECTION OF CODE IMPLEMENTS THE CGM SAMPLER w/ROTATE FOR COMPARISON 

### function to do initialization of data for MCMC sampler 

	n=1000
	m=10
	data_dim=101
	
# SETTING UP DATA STRUCTURES....
t_lh<-matrix(NA,nrow=n*m,ncol=1) #matrix tracks the likelihood of the tree at the current iteration of the MCMC chain
output<-matrix(NA, nrow=n*m, ncol=2)
rtree<-vector(mode="list",length=n*m)

alpha_prime<-matrix(1,nrow=data_dim)
alpha_prime_store<-matrix(NA, nrow=n*m,ncol=length(alpha_prime))
weight<-matrix(NA, nrow=n*m,ncol=length(alpha_prime)-1)
count<-matrix(NA,nrow=n*m, ncol=length(alpha_prime))
#INITIALIZATIONS

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

for(k in 1:10){
	
for(j in 1:m){
	i<-1
	tree_num<-sample(1:hmm,1,prob=x_double_prime)
	tree1<-tree_init_mat[[tree_num]]
	tree2<-tree1 
	#Z[i+1000*(j-1)]<-1
	#alpha_prime_store[i+1000*(j-1) ,]<-alpha_prime
	#weight[i+1000*(j-1),]<-1/length(alpha_prime)
 	t_lh[i+1000*(j-1)]<-tree_lhood(tree2,casper_data,theta=0,mult_penalty=FALSE)
 	#split_counts<-0
for(i in 2:n){
	x<-rgen_discrete_unif()
			if( x == 1 ){ 
			# sample the grow
			tree2<-grow(tree1,casper_data,1)		
			tree2<-up(tree2,casper_data,1)
			}else{
				if(x ==2){
				tree2<-prune(tree1,casper_data,1)
				tree2<-up(tree2,casper_data,1)
			}else{ 
				if(x==3 ){
					tree2<-change(tree1,casper_data,1)		
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
	
	R<-log_ratio_calc(tree2,tree1,t2_lh, t_lh[i+1000*(j-1)-1],data=casper_data)
	
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
	#split_counts<-count_cov_split(tree=tree1,data=casper_data,1) + split_counts 
	#Z[(j-1)*1000+i]<- data_dim/sum(split_counts)
	#split_counts<-(split_counts)*Z[(j-1)*1000+i]/40 # here the dividing number must change as the dimensionality of data does 
	#weight[(j-1)*1000+i,]<-sample_weights(alpha=split_counts) ### the 1's are added inside the function.... 

	}# end of for(i 1:n) loop
} #end of for j in 1:m  loop

# end of outer loopcasper_data

### output in chain_count 
chain_count<-model_jump_count(nchains=10,chain_tree_list=rtree,true_tree_list=true_tree_list)

sum_out<-summary_jumps(chain_count)

old_count<-old_count+sum_out 
}
three<-exp(third_term(t_lh))
min(.5*( sum(abs(exp(-true_tree_lhood-max(best_Z(t_lh)))-.5))+three),.5*(sum(abs(exp(true_tree_lhood+max(best_Z(t_lh)))-.5))+three))
old_count
system("say yo! Ive got your  output")

plot(t_lh,type="l")

boxplot.matrix(weight,use.cols=TRUE, main="Distribution of Weights Across Dimensions")

colMeans(weight)
 
###THIS SECTION OF CODE IMPLEMENTS THE VANILLA CGM SAMPLER FOR COMPARISON 

### function to do initialization of data for MCMC sampler 

	n=1000
	m=10
	data_dim=101
	
# SETTING UP DATA STRUCTURES....
t_lh<-matrix(NA,nrow=n*m,ncol=1) #matrix tracks the likelihood of the tree at the current iteration of the MCMC chain
output<-matrix(NA, nrow=n*m, ncol=2)
rtree<-vector(mode="list",length=n*m)

alpha_prime<-matrix(1,nrow=data_dim)
alpha_prime_store<-matrix(NA, nrow=n*m,ncol=length(alpha_prime))
weight<-matrix(NA, nrow=n*m,ncol=length(alpha_prime)-1)
count<-matrix(NA,nrow=n*m, ncol=length(alpha_prime))
#INITIALIZATIONS

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

for(k in 1:10){
	
for(j in 1:m){
	i<-1
	tree_num<-sample(1:hmm,1,prob=x_double_prime)
	tree1<-tree_init_mat[[tree_num]]
	tree2<-tree1 
	#Z[i+1000*(j-1)]<-1
	#lpha_prime_store[i+1000*(j-1) ,]<-alpha_prime
	# weight[i+1000*(j-1),]<-1/length(alpha_prime)
 	t_lh[i+1000*(j-1)]<-tree_lhood(tree2,casper_data,theta=0,mult_penalty=FALSE)
 	#split_counts<-0
for(i in 2:n){
	x<-rgen_discrete_unif_four()
			if( x == 1 ){ 
			# sample the grow
			tree2<-grow(tree1,casper_data,1)		
			tree2<-up(tree2,casper_data,1)
			}else{
				if(x ==2){
				tree2<-prune(tree1,casper_data,1)
				tree2<-up(tree2,casper_data,1)
			}else{ 
				if(x==3 ){
					tree2<-change(tree1,casper_data,1)		
					tree2<-up(tree2,casper_data,1)
				}else{
					tree2<-swap(tree1,casper_data,1)
					tree2<-up(tree2,casper_data,1)
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
	
	R<-log_ratio_calc(tree2,tree1,t2_lh, t_lh[i+1000*(j-1)-1],data=casper_data)
	
	#R<- t2_lh - c_t_lh[i+1000*(j-1)-1] 
	###End of likelihood ration calculation
	
	### end dof log ratio  calculation
	
	M_E<- -rexp(1,1)
	if(M_E <=R){ 
		t_lh[i+ (j-1)*1000 ]<- t2_lh # make the move to the proposed tree
			# stay at tree2 so do not update the tree2 object
			tree1<-tree2
		#	output[i+(j-1)*1000,1]<-length(unique(tree1$frame[,1]))-1
		#	output[i+(j-1)*1000,2]<-length(get_tnodes(tree1))
		}else{
			# stay at the current tree state 
			t_lh[i+(j-1)*1000]<-t_lh[1000*(j-1)+i-1]
		#	output[i+(j-1)*1000,1]<-length(unique(tree1$frame[,1]))-1
	#		output[i+(j-1)*1000,2]<-length(get_tnodes(tree1))			#num_term_nodes[i+(j-1)*1000-1]
		}
		rtree[[i+(j-1)*1000]]<-tree1
	if(i %%300 ==0 ){
		print(tree1)
		}
	
	#### now sample from the weights for each dimension.... 
	# split_counts<-count_cov_split(tree=tree1,data=casper_data,1) + split_counts 
	# Z[(j-1)*1000+i]<- data_dim/sum(split_counts)
	# split_counts<-(split_counts)*Z[(j-1)*1000+i]/40 # here the dividing number must change as the dimensionality of data does 
	# weight[(j-1)*1000+i,]<-sample_weights(alpha=split_counts) ### the 1's are added inside the function.... 

	}# end of for(i 1:n) loop
} #end of for j in 1:m  loop

# end of outer loopcasper_data

### output in chain_count 
chain_count<-model_jump_count(nchains=10,chain_tree_list=rtree,true_tree_list=true_tree_list)

sum_out<-summary_jumps(chain_count)

old_count<-old_count+sum_out 
}
three<-third_term(t_lh)
min(.5*( sum(abs(exp(-true_tree_lhood-max(best_Z(t_lh)))-.5))+three),.5*(sum(abs(exp(true_tree_lhood+max(best_Z(t_lh)))-.5))+three))
old_count

system("say yo! Ive got your  output")
