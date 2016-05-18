snp.refine.boundary <- function(data, segs.stat)
{	
    ## update change point function
    update.change.point <- function(xl, xm, xr, chptIdx)
    {
    	nxm <- length(xm)
    	if (nxm == 1) { return(chptIdx[1]) }
    	if (nxm >= 2) {
    		tstat <- rep(NA, nxm)
    		for (i in 1:nxm) {
    			left <- c(xl, xm[1:i])
    			right <- ifelse(i==nxm, xr, c(xm[(i+1):nxm], xr) )
    			tstat[i] <- abs(mean(left, na.rm=TRUE) - mean(right, na.rm=TRUE))/sqrt(1/length(left) + 1/length(right))
    		}
    		
    		return( chptIdx[which.max(tstat)[1]] )
    	}
    }
    
    ## assume data is sorted by chr and position
    chrs <- unique(segs.stat$chr)
    segs.stat.refine <- NULL
    for (chr in chrs) {
    	segs.stat.refine.chr <- segs.stat[segs.stat$chr==chr,]
    	data.chr <- data[data$chr==chr,]

        nseg.chr  <- nrow(segs.stat.refine.chr)
    	ndata.chr <- nrow(data.chr)
    	
    	## update index to whole data index
    	idx.seg <- which(data.chr$use.in.seg==1)
    	segs.stat.refine.chr$chrIdxStart <- idx.seg[segs.stat.refine.chr$chrIdxStart]
    	     segs.stat.refine.chr$chrIdxEnd   <- idx.seg[segs.stat.refine.chr$chrIdxEnd]
    	
    	if (nseg.chr==1) {
    		segs.stat.refine.chr$posStart[1]    <- data.chr$position[1]
    		segs.stat.refine.chr$posEnd[1]      <- data.chr$position[ndata.chr]
    		segs.stat.refine.chr$length[1]      <- segs.stat.refine.chr$posEnd[1] - segs.stat.refine.chr$posStart[1] + 1
    		segs.stat.refine.chr$chrIdxStart[1] <- 1
    		segs.stat.refine.chr$chrIdxEnd[1]   <- ndata.chr
    		segs.stat.refine.chr$numProbe[1]    <- ndata.chr
    	}
    	if (nseg.chr>=2) {
    		## start of segment 1
    		segs.stat.refine.chr$posStart[1]    <- data.chr$position[1]
    		segs.stat.refine.chr$chrIdxStart[1] <- 1
    		
    		## change points between consecutive segments    		
    		for (i in 1:(nseg.chr-1)) {
    			lstart <- segs.stat.refine.chr$chrIdxStart[i]
    			lend   <- segs.stat.refine.chr$chrIdxEnd[i] - 1
    			xl     <- data.chr$log2ratio[lstart:lend]
    			
    			mstart <- segs.stat.refine.chr$chrIdxEnd[i]
    			mend   <- segs.stat.refine.chr$chrIdxStart[i + 1] - 1
    			xm     <- data.chr$log2ratio[mstart:mend]
    			
    			rstart <- segs.stat.refine.chr$chrIdxStart[i + 1]
    			rend   <- segs.stat.refine.chr$chrIdxEnd[i + 1]
    			xr     <- data.chr$log2ratio[rstart:rend]
    			
    			chptIdx <- mstart:mend
    			updated.chpt <- update.change.point(xl=xl, xm=xm, xr=xr, chptIdx=chptIdx)
    			
    		    if (length(updated.chpt)==1) {
    		        ## update right end of segment i
    		        segs.stat.refine.chr$posEnd[i]      <- data.chr$position[updated.chpt]
    		        segs.stat.refine.chr$length[i]      <- segs.stat.refine.chr$posEnd[i] - segs.stat.refine.chr$posStart[i] + 1
    		        segs.stat.refine.chr$chrIdxEnd[i]   <- updated.chpt
    		        segs.stat.refine.chr$numProbe[i]    <- segs.stat.refine.chr$chrIdxEnd[i] - segs.stat.refine.chr$chrIdxStart[i] + 1
    		        
    		        ## update left end of segment i+1
    		        segs.stat.refine.chr$posStart[i+1]  <- data.chr$position[updated.chpt + 1]
    		        segs.stat.refine.chr$chrIdxStart[i+1] <- updated.chpt + 1
    		    } else {
    		    	## change point not change, just update length and numProbe for segment i
    		        segs.stat.refine.chr$length[i]      <- segs.stat.refine.chr$posEnd[i] - segs.stat.refine.chr$posStart[i] + 1
    		        segs.stat.refine.chr$numProbe[i]    <- segs.stat.refine.chr$chrIdxEnd[i] - segs.stat.refine.chr$chrIdxStart[i] + 1
    		    }
    		} ## for i

    		## end of segment nseg.chr
    		segs.stat.refine.chr$posEnd[nseg.chr]    <- data.chr$position[ndata.chr]
    		segs.stat.refine.chr$length[nseg.chr]    <- segs.stat.refine.chr$posEnd[nseg.chr] - segs.stat.refine.chr$posStart[nseg.chr] + 1
    		segs.stat.refine.chr$chrIdxEnd[nseg.chr] <- ndata.chr
    		segs.stat.refine.chr$numProbe[nseg.chr]  <- segs.stat.refine.chr$chrIdxEnd[nseg.chr] - segs.stat.refine.chr$chrIdxStart[nseg.chr] + 1
    	} ## if (nseg.chr>=2)

    segs.stat.refine <- rbind(segs.stat.refine, segs.stat.refine.chr)
    } ## for chr
    
   return( segs.stat.refine ) 
}
