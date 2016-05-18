cnv.call <-
function(data, sample.id, segs.stat, maxL=NULL, N=1000,
         pvalue.cutoff=0.05, seed=NULL,
         do.manual.baseline=FALSE,
         log2mBAF.left=NULL, log2mBAF.right=NULL, log2ratio.bottom=NULL, log2ratio.up=NULL) 
         ##qvalue.cutoff=0.05 maybe used in the future
{ 
	## set seed for reproducibility
	if (!is.null(seed)) set.seed(seed)
	
	## functions
	seg.stat.df1.bootstrap <- function(x, n, stat1, sigma1)
    {
    	l1 <- x; r1 <- x+n-1
    	stat1.seg <- median(stat1[l1:r1], na.rm=TRUE)
    	tmp1 <- sigma1*sqrt(1/n)
    	return( (stat1.seg/tmp1)^2 )
    }
	
	seg.stat.df2.bootstrap <- function(x, n, stat1, stat2, sigma1, sigma2)
    {
    	l1 <- x; r1 <- x+n-1
    	stat1.seg <- median(stat1[l1:r1], na.rm=TRUE)
        stat2.seg <- median(stat2[l1:r1], na.rm=TRUE)
    	tmp1 <- sigma1*sqrt(1/n)
    	tmp2 <- sigma2*sqrt(1/n)
    	return( (stat1.seg/tmp1)^2 + (stat2.seg/tmp2)^2 )
    }

    compute.df1.pvalue <- function(seg1.stat, do.log2ratio=TRUE, data.null, maxL=maxL, N=N)
    {
        n <- ifelse(seg1.stat$numProbe>maxL, maxL, seg1.stat$numProbe)
        stat.tmp <- ifelse(do.log2ratio, seg1.stat$log2ratio.Median.adj, seg1.stat$log2mBAF.Median.adj )
        sigma.tmp <- ifelse(do.log2ratio, sigma.log2ratio*sqrt(1/n), sigma.log2mBAF*sqrt(1/n) )
        obs.stat <- (stat.tmp/sigma.tmp)^2
        x <- sample(x=1:(nrow(data.null)-n), size=N, replace=TRUE)
        if (do.log2ratio) {
            y <- sapply(X=x,FUN=seg.stat.df1.bootstrap, n=n, stat1=data.null$log2ratio.adj, sigma1=sigma.log2ratio)        	
        } else {
        	y <- sapply(X=x,FUN=seg.stat.df1.bootstrap, n=n, stat1=data.null$log2mBAF.adj, sigma1=sigma.log2mBAF)
        }
        return( sum(y>obs.stat)/N )
    }    

    compute.df2.pvalue <- function(seg1.stat, data.null, maxL=maxL, N=N)
    {
        n <- ifelse(seg1.stat$numProbe>maxL, maxL, seg1.stat$numProbe)
        log2ratio.1 <- seg1.stat$log2ratio.Median.adj
        log2mBAF.1  <- seg1.stat$log2mBAF.Median.adj
        tmp1 <- sigma.log2ratio*sqrt(1/n)
    	tmp2 <- sigma.log2mBAF*sqrt(1/n)
        obs.stat <- (log2ratio.1/tmp1)^2 + (log2mBAF.1/tmp2)^2
        x <- sample(x=1:(nrow(data.null)-n), size=N, replace=TRUE)
        y <- sapply(X=x,FUN=seg.stat.df2.bootstrap, n=n, stat1=data.null$log2ratio.adj, stat2=data.null$log2mBAF.adj, 
                                                         sigma1=sigma.log2ratio, sigma2=sigma.log2mBAF)
        return( sum(y>obs.stat)/N )
    }    
    
    ## compute log2ratio baseline
    sigma.log2ratio <- delta.sd(data$log2ratio)
    sigma.log2mBAF  <- delta.sd(data$log2mBAF)
    if (do.manual.baseline) {
    	idx.base <- compute.baseline.manual(segs.stat=segs.stat, data=data,
    	                                    log2mBAF.left=log2mBAF.left, 
    	                                    log2mBAF.right=log2mBAF.right, 
    	                                    log2ratio.bottom=log2ratio.bottom, 
    	                                    log2ratio.up=log2ratio.up)
    } else {
    	idx.base <- compute.baseline(segs.stat=segs.stat, data=data)
    }
    maxL <- ifelse(is.null(maxL), floor(nrow(data)/100), maxL)
    
    ## sample info
    segs.stat$Sample_ID  <- sample.id
    segs.stat$remark     <- ifelse(length(idx.base)>nrow(data)/10, 0, 1)
    segs.stat$log2ratio.base.Mean   <- mean(data$log2ratio[idx.base], na.rm=TRUE)
    segs.stat$log2ratio.base.Median <- median(data$log2ratio[idx.base], na.rm=TRUE)
    segs.stat$log2ratio.Sigma  <- sigma.log2ratio
    segs.stat$log2mBAF.base.Mean   <- mean(data$log2mBAF[idx.base], na.rm=TRUE)
    segs.stat$log2mBAF.base.Median <- median(data$log2mBAF[idx.base], na.rm=TRUE)
    segs.stat$log2mBAF.Sigma   <- sigma.log2mBAF
    segs.stat$log2ratio.Mean.adj   <- segs.stat$log2ratio.Mean   - segs.stat$log2ratio.base.Mean
    segs.stat$log2ratio.Median.adj <- segs.stat$log2ratio.Median - segs.stat$log2ratio.base.Median
    segs.stat$log2mBAF.Mean.adj    <- segs.stat$log2mBAF.Mean    - segs.stat$log2mBAF.base.Mean
    segs.stat$log2mBAF.Median.adj  <- segs.stat$log2mBAF.Median  - segs.stat$log2mBAF.base.Median

    data$log2ratio.adj <-  data$log2ratio - unique(segs.stat$log2ratio.base.Median)
    data$log2mBAF.adj  <-  data$log2mBAF  - unique(segs.stat$log2mBAF.base.Median)
    data.null <- data[idx.base,]
    if (nrow(data.null) <= max(segs.stat$numProbe)) maxL <- floor(nrow(data.null)/2)
    
    ## compute p-value
    segs.stat$log2ratio.p.value <- NA
    segs.stat$log2mBAF.p.value <- NA
    segs.stat$p.value <- NA
    for (i in 1:nrow(segs.stat)) {
    	segs.stat$log2ratio.p.value[i] <- compute.df1.pvalue(seg1.stat=segs.stat[i,], do.log2ratio=TRUE,  data.null=data.null, maxL=maxL, N=N)
    	segs.stat$log2mBAF.p.value[i]  <- compute.df1.pvalue(seg1.stat=segs.stat[i,], do.log2ratio=FALSE, data.null=data.null, maxL=maxL, N=N)    	
    	segs.stat$p.value[i] <- compute.df2.pvalue(seg1.stat=segs.stat[i,], data.null=data.null, maxL=maxL, N=N)
    }
    
    ## decide CNV status
    segs.stat$CNV <- NA
    for (i in 1:nrow(segs.stat)) {
    	if (segs.stat$log2mBAF.p.value[i] >= pvalue.cutoff) {
    		if (segs.stat$log2ratio.p.value[i] >= pvalue.cutoff) segs.stat$CNV[i] <- "normal"
            if (segs.stat$log2ratio.p.value[i] <  pvalue.cutoff) {
            	segs.stat$CNV[i] <- "undecided"
            	## reassign back to normal
            	if (abs(segs.stat$log2ratio.Median.adj[i]) < 0.5 * sigma.log2ratio  |  segs.stat$p.value[i] >= pvalue.cutoff) {
            		segs.stat$CNV[i] <- "normal"
            	}
            	## balanced gain
            	if (segs.stat$log2ratio.Median.adj[i] > 1.5 * sigma.log2ratio  &  segs.stat$p.value[i] < pvalue.cutoff) {
            		segs.stat$CNV[i] <- "gain"
            	}
            	## balanced loss in high ploidy genome
            	if (segs.stat$log2ratio.Median.adj[i] <  -1.5 * sigma.log2ratio  &  segs.stat$p.value[i] < pvalue.cutoff) {
            		segs.stat$CNV[i] <- "loss"
            	}
            }
        }
    	if (segs.stat$log2mBAF.p.value[i] < pvalue.cutoff & segs.stat$p.value[i] < pvalue.cutoff) {
    		if (segs.stat$log2mBAF.Median.adj[i] <= -0.5*sigma.log2mBAF) segs.stat$CNV[i] <- "undecided"
    		if (segs.stat$log2mBAF.Median.adj[i] > -0.5*sigma.log2mBAF & segs.stat$log2mBAF.Median.adj[i] <= 0) {
    			segs.stat$CNV[i] <- "undecided"
    			## reassign back to normal
            	if (abs(segs.stat$log2ratio.Median.adj[i]) < 0.5 * sigma.log2ratio  |  segs.stat$p.value[i] >= pvalue.cutoff) {
            		segs.stat$CNV[i] <- "normal"
            	}
            	## balanced gain
            	if (segs.stat$log2ratio.Median.adj[i] > 1.5 * sigma.log2ratio  &  segs.stat$p.value[i] < pvalue.cutoff) {
            		segs.stat$CNV[i] <- "gain"
            	}
            	## balanced loss in high ploidy genome
            	if (segs.stat$log2ratio.Median.adj[i] <  -1.5 * sigma.log2ratio  &  segs.stat$p.value[i] < pvalue.cutoff) {
            		segs.stat$CNV[i] <- "loss"
            	}
    	    }
    		if (segs.stat$log2mBAF.Median.adj[i] > 0 & segs.stat$log2ratio.p.value[i] >= pvalue.cutoff) {
    			segs.stat$CNV[i] <- "LOH"
    		    }
    		if (segs.stat$log2mBAF.Median.adj[i] > 0 & segs.stat$log2ratio.p.value[i] < pvalue.cutoff) {
    			if (segs.stat$log2ratio.Median.adj[i] > 0) segs.stat$CNV[i] <- "gain"
    			if (segs.stat$log2ratio.Median.adj[i] < 0) segs.stat$CNV[i] <- "loss"
    		    }
    	}
    } ## for i
    segs.stat$CNV[is.na(segs.stat$CNV)] <- "undecided"

    return(segs.stat)
}
