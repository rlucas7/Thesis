#####
#####		This script runs a depth restricted tree sampler. Use this script with several values 
#####		of the maximum depth to perform simulation study of max depth and the sensitivity to 
#####		choice of this max depth value. 
#####
#####		Date FRiday dec 7th 2012			
#####-------------------------------------------------------------

####extra_data is the name of the dataframe that contains the CGM data with additional covariates 
getwd()
setwd("/Users/rlucas7/Desktop/depth_calibration_output_CGM_example")

c<-10 # max depth 
	n=1000
	m=3
	data_dim=100
M<-1000 # the max lag to calculate in the ACf output analysis.	
# SETTING UP DATA STRUCTURES....
t_lh<-matrix(NA,nrow=n*m,ncol=1) #matrix tracks the likelihood of the tree at the current iteration of the MCMC chain
output<-matrix(NA, nrow=n*m, ncol=3)
rtree<-vector(mode="list",length=n*m)

alpha_prime<-matrix(1,nrow=data_dim)
alpha_prime_store<-matrix(NA, nrow=n*m,ncol=length(alpha_prime))
weight<-matrix(1, nrow=m*n,ncol=length(alpha_prime))
#count<-matrix(NA,nrow=n*m, ncol=length(alpha_prime))
#INITIALIZATIONS

t_lh[1]<-tree_lhood(t1,extra_data,theta=0,mult_penalty=FALSE)
output[1,1]<-length(unique(t1$frame[,1]))-1
output[1,2]<-length(get_tnodes(t1))
rtree[[1]]<-t1

Z<-matrix(NA, nrow=n*m, ncol=1) #here the "Z" vector represents the constraint multipliers $\widetilde{\alpha}$
t2<-t1 # the "optimal" tree found by the MCMC search   ////////// ????Is this line still necessary for the code ?????
Z[1]<-1

alpha_prime[100]<-100
alpha_prime[99]<-100

alpha_prime_store[1,]<-alpha_prime
#weight[1,]<-1/sum(alpha_prime)

split_counts<-0

for(j in 1:m){
	i<-1
	#tree_num<-sample(1:hmm,1,prob=x_double_prime)
	#tree1<-tree_init_mat[[tree_num]]
	tree1<-snip.tree(t1,nodes=c(2,3))
	tree2<-tree1 
	split_counts<-alpha_prime
	Z[i+1000*(j-1)]<-1
	alpha_prime_store[i+1000*(j-1),]<-alpha_prime
	weight[i+1000*(j-1),]<-1/length(alpha_prime)
 	t_lh[i+1000*(j-1)]<-tree_lhood(tree2,extra_data,theta=0,mult_penalty=FALSE)
 
for(i in 2:n){
	x<-rgen_discrete_unif()
			if( x == 1 ){ 
			# sample the grow
			tree2<-grow(tree1,extra_data,101,weight[i-1+1000*(j-1),])		
			tree2<-up(tree2,extra_data,101)
			}else{
				if(x ==2){
				tree2<-prune(tree1,extra_data,101)
				tree2<-up(tree2,extra_data,101)
			}else{ 
				if(x==3 ){
					tree2<-change(tree1,extra_data,101,weight[i-1+1000*(j-1),])		
					tree2<-up(tree2,extra_data,101)
				}else{
					if(x==4 ){
						#rotate step...
						tree2<-rotate_k(tree1)
						tree2<-up(tree2,extra_data,101)	
						}else{ # x==5 case
						tree2<-swap(tree1,extra_data,101)
						tree2<-up(tree2,extra_data,101)
						#include swap function here, .....
						}
					}	
				}
			}
	
	t2_lh<-tree_lhood(tree2,extra_data,theta=0,mult_penalty=FALSE)	# calculate the proposed tree's l-hood	
	R<-log_ratio_calc(tree2,tree1,t2_lh, t_lh[i+1000*(j-1)-1],data=extra_data,weight=weight[i-1+1000*(j-1),])
	M_E<- -rexp(1,1)
	if(M_E <=R && depth_calc(tree2) <= c){ 
		t_lh[i+ (j-1)*1000 ]<- t2_lh # make the move to the proposed tree
			# stay at tree2 so do not update the tree2 object
			tree1<-tree2
			output[i+(j-1)*1000,1]<-length(unique(tree1$frame[,1]))-1
			output[i+(j-1)*1000,2]<-length(get_tnodes(tree1))
			output[i+(j-1)*1000,3]<-depth_calc(tree1)			#num_term_nodes[i+(j-1)*1000-1]
		}else{
			# stay at the current tree state 
			t_lh[i+(j-1)*1000]<-t_lh[1000*(j-1)+i-1]
			output[i+(j-1)*1000,1]<-length(unique(tree1$frame[,1]))-1
			output[i+(j-1)*1000,2]<-length(get_tnodes(tree1))			#num_term_nodes[i+(j-1)*1000-1]
			output[i+(j-1)*1000,3]<-depth_calc(tree1)			#num_term_nodes[i+(j-1)*1000-1]
		}
		rtree[[i+(j-1)*1000]]<-tree1
	if(i %%500 ==0 ){
		print(tree1)
		}
	#### now sample from the weights for each dimension.... 
	split_counts<-count_cov_split(tree=tree1,data=extra_data,101) + split_counts 
	#Z[(j-1)*1000+i]<- data_dim/sum(split_counts)
	#split_counts<-(split_counts)*Z[(j-1)*1000+i]/10 # here the dividing number must change as the dimensionality of data does 
	weight[(j-1)*1000+i,]<-sample_weights(alpha=split_counts) ### the 1's are added inside the function.... 
	#weight[i-1+1000*(j-1),]<-(split_counts-1)/(sum(split_counts)-1000)
	}# end of for(i 1:n) loop
} #end of for j in 1:m  loop

max_tree<-match(max(t_lh[-c(1:500,1000:1500,2000:2500)]),t_lh[-c(1:500,1000:1500,2000:2500)])
#par(mfcol=c(2,4))

pdf(file="depth_output_c=10.pdf")
par(mfcol=c(1,1))
plot(t_lh, type="l")
plot(rtree[[max_tree]],type="unif")
text(rtree[[max_tree]])
plot(output[,1], type="l",ylab="total number of nodes")
plot(output[,2], type="l", ylab="number of terminal nodes")
plot(output[,3],type="l", ylab="depth of tree")
boxplot.matrix(weight,use.cols=TRUE)
acf_out<-acf(t_lh,lag.max=3000)
dev.off()

z<-0
bf<-TRUE
for(i in 1:M){
if(sign(acf_out$acf[i]) !=sign(acf_out$acf[i+1]) && bf==TRUE){z<<-i; bf<-FALSE}
}
z
IACF<-0 # integrated auto-correlation function
if(z==0){IACF<-sum(acf_out$acf[1:1000])}else{IACF<-sum(acf_out$acf[1:z])}
IACF
write.csv(IACF,file="IACF_c=10.csv",append=FALSE)

colMeans(weight)

write.csv(colMeans(weight[c(500:1000,1500:2000,2500:3000),]),file="col_means_c=10.csv")

par(ask=FALSE)

for(i in seq(1000,3000,by=10)){
	plot(rtree[[i]],type="unif")
	text(rtree[[i]])
}