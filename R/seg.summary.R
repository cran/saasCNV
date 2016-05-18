seg.summary <-
function (data, segs, chr)
{ 
    cnv <- as.data.frame(matrix(NA,nrow(segs),15))
    names(cnv) <- c("chr","posStart","posEnd","length","chrIdxStart","chrIdxEnd","numProbe",
                    "log2ratio.Mean","log2ratio.SD","log2ratio.Median","log2ratio.MAD",
                    "log2mBAF.Mean","log2mBAF.SD","log2mBAF.Median","log2mBAF.MAD")

    for (i in 1:nrow(segs)) {
        idx.start <- segs[i,1]; idx.end <- segs[i,2]
        idx  <- idx.start:idx.end
        start <- data$pos[idx.start]; end <- data$pos[idx.end]
        
 	    ## seg info
 	    cnv$chr[i]          <- chr
 	    cnv$posStart[i]     <- start
 	    cnv$posEnd[i]       <- end
 	    cnv$length[i]       <- end - start + 1
 	    cnv$chrIdxStart[i]  <- idx.start
 	    cnv$chrIdxEnd[i]    <- idx.end
 	    cnv$numProbe[i]     <- length(idx)
    
 	    ## statistics
 	    cnv$log2ratio.Mean[i]    <- mean(data$log2ratio[idx],na.rm=T)
 	    cnv$log2ratio.SD[i]      <- sd(data$log2ratio[idx],na.rm=T) 
 	    cnv$log2ratio.Median[i]  <- median(data$log2ratio[idx],na.rm=T)
   	    cnv$log2ratio.MAD[i]     <- mad(data$log2ratio[idx],na.rm=T)
   	    cnv$log2mBAF.Mean[i]     <- mean(data$log2mBAF[idx],na.rm=T)
   	    cnv$log2mBAF.SD[i]       <- sd(data$log2mBAF[idx],na.rm=T)
   	    cnv$log2mBAF.Median[i]   <- median(data$log2mBAF[idx],na.rm=T)
   	    cnv$log2mBAF.MAD[i]      <- mad(data$log2mBAF[idx],na.rm=T)
    }

    return(cnv)
}
