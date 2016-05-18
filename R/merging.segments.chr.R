merging.segments.chr <-
function(data.chr, data.null, segs.chr, chr, sigma.log2ratio, sigma.log2mBAF, N=1000, maxL=2000, p.value.cut.off=0.05)
{
	## segs.chr must be sorted
	seg.diff.bootstrap <- function(x, n1, n2, stat1, stat2, sigma1, sigma2)
    {
    	l1 <- x; r1 <- x+n1-1
    	l2 <- r1+1; r2 <- r1+n2
    	stat1.l <- median(stat1[l1:r1], na.rm=TRUE)
    	stat1.r <- median(stat1[l2:r2], na.rm=TRUE)
        stat2.l <- median(stat2[l1:r1], na.rm=TRUE)
    	stat2.r <- median(stat2[l2:r2], na.rm=TRUE)
    	tmp1 <- sigma1*sqrt(1/n1+1/n2)
    	tmp2 <- sigma2*sqrt(1/n1+1/n2)
    	return( ((stat1.l - stat1.r)/tmp1)^2 + ((stat2.l - stat2.r)/tmp2)^2 )
    }

    compute.delta.pvalue <- function(idx1.l, idx1.r, idx2.l, idx2.r, maxL)
    {
        n1 <- idx1.r - idx1.l + 1; n1 <- ifelse(n1>maxL, maxL, n1)
        n2 <- idx2.r - idx2.l + 1; n2 <- ifelse(n2>maxL, maxL, n2)
        log2ratio.1 <- median(data.chr$log2ratio[idx1.l:idx1.r], na.rm=TRUE)
        log2ratio.2 <- median(data.chr$log2ratio[idx2.l:idx2.r], na.rm=TRUE)
        log2mBAF.1  <- median(data.chr$log2mBAF[idx1.l:idx1.r], na.rm=TRUE)
        log2mBAF.2  <- median(data.chr$log2mBAF[idx2.l:idx2.r], na.rm=TRUE)
        tmp1 <- sigma.log2ratio*sqrt(1/n1+1/n2)
    	tmp2 <- sigma.log2mBAF*sqrt(1/n1+1/n2)
        obs.stat <- ((log2ratio.1 - log2ratio.2)/tmp1)^2 + 
                    ((log2mBAF.1 - log2mBAF.2)/tmp2)^2
        x <- sample(x=1:(nrow(data.null)-n1-n2), size=N, replace=TRUE)
        y <- sapply(X=x,FUN=seg.diff.bootstrap, n1=n1, n2=n2, stat1=data.null$log2ratio, stat2=data.null$log2mBAF, 
                                                sigma1=sigma.log2ratio, sigma2=sigma.log2mBAF)
        delta.pvalue <- sum(y>obs.stat)/N
        return(list(obs.stat=obs.stat, delta.pvalue=delta.pvalue))
    }

    ## initialize
    segs.chr$stat <- NA
    segs.chr$delta.pvalue <- NA
    
    if (nrow(segs.chr)>1) {
        for (i in 1:(nrow(segs.chr)-1) ) {
        	tmp <-  compute.delta.pvalue(idx1.l=segs.chr[i,1], idx1.r=segs.chr[i,2], 
        	                             idx2.l=segs.chr[i+1,1], idx2.r=segs.chr[i+1,2], maxL=maxL)
            segs.chr$stat[i] <- tmp$obs.stat
            segs.chr$delta.pvalue[i] <- tmp$delta.pvalue
        }
    
        ## merge recursively starting from the least significant two segments
        max.pvalue <- max(segs.chr$delta.pvalue, na.rm=TRUE)
        idx  <- which.max(segs.chr$delta.pvalue)
        while (nrow(segs.chr)>1 & max.pvalue>p.value.cut.off) {
        	segs.chr[idx, "chrIdxEnd"] <- segs.chr[idx+1, "chrIdxEnd"]
        	segs.chr <- segs.chr[-c(idx+1),]
        	## update stat and delta.pvalue before and after
        	if (idx>1) {
                tmp <-  compute.delta.pvalue(idx1.l=segs.chr[idx-1,1], idx1.r=segs.chr[idx-1,2], 
        	                                 idx2.l=segs.chr[idx,1], idx2.r=segs.chr[idx,2], maxL=maxL)
                segs.chr$stat[idx-1] <- tmp$obs.stat
                segs.chr$delta.pvalue[idx-1] <- tmp$delta.pvalue
        	}
        	if (idx<nrow(segs.chr)) {
                tmp <-  compute.delta.pvalue(idx1.l=segs.chr[idx,1], idx1.r=segs.chr[idx,2], 
        	                                 idx2.l=segs.chr[idx+1,1], idx2.r=segs.chr[idx+1,2], maxL=maxL)
                segs.chr$stat[idx] <- tmp$obs.stat
                segs.chr$delta.pvalue[idx] <- tmp$delta.pvalue
        	}
        	if (idx==nrow(segs.chr)) {
                segs.chr$stat[idx] <- NA
                segs.chr$delta.pvalue[idx] <- NA
        	}
            if (nrow(segs.chr)>1) {
                max.pvalue <- max(segs.chr$delta.pvalue, na.rm=TRUE)
                idx  <- which.max(segs.chr$delta.pvalue)
            }
        }
    } ## if nrow
 
    return(segs.chr)
}
