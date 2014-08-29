#### CGM 

niter<-nc

names(extra_data) # first entry is the response. 

tree_mat<-vector(mode="list", length=dim(extra_data)[2]-1)
for(i in 1:(dim(extra_data)[2]-1) ){
tree_mat[[i]]<-tree(paste(names(extra_data)[101],names(extra_data)[i],sep="~"),data=extra_data)
}

for(i in 1:niter){
	if(attr(tree_mat[[i]],"class")[1]=="tree"){
tree_mat[[i]]<-snip.tree(snip.tree(tree_mat[[i]],nodes=2),nodes=3)
	}
}

# bunch of warnings no big deal (nbd)
log_mat<-matrix(NA,nrow=niter, ncol=1)
for(i in 1:niter){
if(attr(tree_mat[[i]],"class")[1]=="tree"){log_mat[i]<-TRUE}else{log_mat[i]<-FALSE}
}

total_num_of_init_trees<-sum(1*log_mat)

tree_init_mat<-vector(mode="list", length=total_num_of_init_trees)

j<-1
for(i in 1:niter){
	if(log_mat[i]){
		tree_init_mat[[j]]<-tree_mat[[i]]
		 j<-j+1
		 }
}

hmm<-length(tree_init_mat)
#### now evolve each initial tree ahead 10 steps and use these as initializers.... 
n<-10
m<-1
j<-1
for(k in 1:hmm){
t_lh<-matrix(NA, 10,1)
rtree<-vector(mode="list", length=10)
tree1<-tree_init_mat[[k]]
tree2<-tree1
t_lh[1]<-tree_lhood(tree1,data=extra_data)
rtree[[1]]<-tree1
for(i in 2:n){
	x<-rgen_discrete_unif()
			if( x == 1 ){ 
			# sample the grow
			tree2<-grow(tree1,extra_data,101)		
			tree2<-up(tree2,extra_data,101)
			}else{
				if(x ==2){
				tree2<-prune(tree1,extra_data,101)
				tree2<-up(tree2,extra_data,101)
			}else{ 
				if(x==3 ){
					tree2<-change(tree1,extra_data,101)		
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
	R<-log_ratio_calc(tree2,tree1,t2_lh, t_lh[i+1000*(j-1)-1],data=extra_data)
	M_E<- -rexp(1,1)
	if(M_E <=R){ 
		t_lh[i+ (j-1)*1000 ]<- t2_lh # make the move to the proposed tree
			# stay at tree2 so do not update the tree2 object
			tree1<-tree2
		}else{
			# stay at the current tree state 
			t_lh[i+(j-1)*1000]<-t_lh[1000*(j-1)+i-1]
		}
		rtree[[i+(j-1)*1000]]<-tree1
	if(i %%2 ==0 ){
		print(tree1)
		}
	}# end of for(i 1:n) loop
tree_init_mat[[k]]<-rtree[[10]]
}



###calculate tree log likelihoods 
log_lh_mat<-matrix(NA, hmm,1)
for(i in 1:hmm){
log_lh_mat[i]<-tree_lhood(tree_init_mat[[i]], data=extra_data)
}

### calculation of initialization probabilities... 

x_plus<-sum(abs(log_lh_mat))
x_prime<-x_plus-abs(log_lh_mat)
x_double_prime<-x_prime/sum(x_prime)


### now to choose which tree to initialize with we do the sampling: 

sample(1:hmm,1,prob=x_double_prime)

f1<-as.formula("y~x100+x99")

tree(f1,data=extra_data)


