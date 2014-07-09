depth_calc<-function(tree_object){
	output<-floor(log(max(as.numeric(rownames(tree_object$frame))),base=2))
	return(output)
}