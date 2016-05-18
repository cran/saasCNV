joint.segmentation <- 
function(data, min.snps=10, global.pval.cutoff=1e-4, max.chpts=30, verbose=TRUE)
{
    chrs <- sub("^chr","",unique(data$chr))

    segs.stat <- NULL
    for (chr in chrs) {
	    ## data in one chr
	    data.chr <- data[data$chr==paste0("chr",chr),]
    
        ## joint segmentation
        y <- cbind(data.chr$log2ratio, data.chr$log2mBAF)    
        y.mscbs <- fmscbs(y=y, win=nrow(y), 
                          MIN.SNPs=min.snps, 
                          GLOBAL.PVAL.CUTOFF=global.pval.cutoff, 
                          MAX.CHPTS=max.chpts)
        if (length(y.mscbs$chpts)>0) {
            segs <- y.mscbs$feature.regions
        } else {
    	    segs <- rbind(c(1,nrow(y)))
        }
    
        ## summarize first round segmentation results
        segs.stat.chr <- seg.summary(data=data.chr, segs=segs, chr=paste0("chr",chr))
    
        ## collect information
        segs.stat <- rbind(segs.stat, segs.stat.chr)
    
        if(verbose) cat("finish segmentation of chr =",chr,"\n")
    } ## for chr
    
    return(segs.stat)
}
