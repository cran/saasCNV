merging.segments <- 
function(data, segs.stat,  use.null.data=TRUE, 
         N=1000, maxL=NULL, merge.pvalue.cutoff=0.05,
         do.manual.baseline=FALSE,
         log2mBAF.left=NULL, log2mBAF.right=NULL, log2ratio.bottom=NULL, log2ratio.up=NULL, 
         seed=NULL, verbose=TRUE)
{
    ## set seed for reproducibility
	if (!is.null(seed)) set.seed(seed)

  	cat("merge adjacent segments not far apart ...\n")
    sigma.log2ratio <- delta.sd(data$log2ratio)
    sigma.log2mBAF  <- delta.sd(data$log2mBAF)
    chrs <- sub("^chr","",unique(data$chr))
    maxL <- ifelse(is.null(maxL), floor(nrow(data)/100), maxL)
    
    if (use.null.data) {
    	if (do.manual.baseline) {
    	    idx.base <- compute.baseline.manual(segs.stat=segs.stat, data=data,
    	                                        log2mBAF.left=log2mBAF.left, 
    	                                        log2mBAF.right=log2mBAF.right, 
    	                                        log2ratio.bottom=log2ratio.bottom, 
    	                                        log2ratio.up=log2ratio.up)
    	    data.null <- data[idx.base,]
    	} else {
    	    idx.base <- compute.baseline(segs.stat=segs.stat, data=data)
    	    data.null <- data[idx.base,]    		
    	}
    } else {
    	data.null <- data
    }
    
    chrs.len <- as.integer(by(segs.stat$numProbe, segs.stat$chr, sum))
    if (nrow(data.null) <= max(chrs.len)) maxL <- floor(nrow(data.null)/3)
    
    ## boostrap L1 + L2 segments consecutively
    segs.stat.merge <- NULL
    for (chr in chrs) {
        data.chr <- data[data$chr==paste0("chr",chr),]
        segs.chr <- segs.stat[segs.stat$chr==paste0("chr",chr), c("chrIdxStart","chrIdxEnd")]
        segs.chr.merge <- merging.segments.chr(data.chr=data.chr, data.null=data.null,
                                               segs.chr=segs.chr, chr=chr, 
                                               sigma.log2ratio=sigma.log2ratio, 
                                               sigma.log2mBAF=sigma.log2mBAF, 
                                               N=N, maxL=maxL, p.value.cut.off=merge.pvalue.cutoff)
    
        ## summarize segmentation results
        segs.stat.chr <- seg.summary(data=data.chr, segs=segs.chr.merge[,1:2], chr=paste0("chr",chr))

        ## collect information
        segs.stat.merge <- rbind(segs.stat.merge, segs.stat.chr)
            
        if (verbose) cat("finish merging of chr =",chr,"\n")
    } ## for chr
    
    return(segs.stat.merge)
}

