BAF2mBAF <- 
function(x, gt=NULL, het.gt=1, mBAF.thd=0.97, win.thd=0.8, w=1, k=2, median.adjust=FALSE, impute=FALSE)
{
	## require("RANN")
	n <- length(x)   # number of SNPs
	adjBAF <- NA
	## median adjust not recommended
	
	## without paired information
	if (is.null(gt)) {
		## reflect
    	mBAF <- abs(x-0.5) + 0.5
	    ## mBAF.thd
	    idx <- sort(which(mBAF < mBAF.thd))
	    idx.na <- sort(which(mBAF >= mBAF.thd))
	    
	    ## win.thd
	    mBAF1 <- mBAF[idx]
	    n1 <- length(mBAF1)
	    win.sum <- abs(mBAF1-0.5)
	    for (i in 1:w) {
	        win.sum[1:(n1-i)] <- win.sum[1:(n1-i)] + abs(mBAF1[1:(n1-i)]-mBAF1[(i+1):n1])
	        win.sum[(i+1):n1] <- win.sum[(i+1):n1] + abs(mBAF1[(i+1):n1]-mBAF1[1:(n1-i)])
	        }
	    idx1.na <- idx[win.sum >= win.thd]
	    idx.na <- sort(union(idx.na, idx1.na))
	    idx <- sort(setdiff(idx, idx1.na))
	}
	
	## with paired information
	if (!is.null(gt)) {
		## het gt criteria
		idx <- sort(which(gt %in% het.gt & x > 0.01 & x < 0.99))
	    idx.na <- sort(setdiff(1:n,idx))		
		
		## median adjust
        if (median.adjust) {
	        me <- median(x[idx]) - 0.5
	        x[idx] <- x[idx] - me
            adjBAF <- x
            }

        ## reflect
    	mBAF <- abs(x-0.5) + 0.5
	}
		
	## impute initial value for missing data
	if (impute) {
	    mBAF1 <- mBAF[idx]
	    sigma <- delta.sd(mBAF1)
	
        tmp <- nn2(data=cbind(idx),query=cbind(idx.na),k=k)   ## k-nearest neighbor
        idx.nn <- idx[tmp$nn.idx]
        mBAF[idx.na] <- rowMeans(matrix(mBAF[idx.nn],nrow=length(idx.na))) + rnorm(length(idx.na),0,sigma)
    }
    
    if (!impute) mBAF[idx.na] <- NA
  
	return(list(mBAF=mBAF, adjBAF=adjBAF, idx=idx, idx.na=idx.na))
}

