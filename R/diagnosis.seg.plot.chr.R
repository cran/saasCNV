diagnosis.seg.plot.chr <-
function(data, segs, sample.id="Sample", chr=1, cex=0.3)
{
	data <- data[data$chr==paste0("chr",chr),]
	segs <- segs[segs$chr==paste0("chr",chr),]
	
	oldpar <- par(mfrow=c(4,1))
    ## log2 ratio
    par(mar=c(4,4,4,1))
    plot(data$position/1e6, data$log2ratio, pch=20, col="gray", cex=cex,
         xlab="Position (Mb)",ylab="log2ratio",
         main=paste0(sample.id," Chr ",chr)) 
    for (j in 1:nrow(segs)) {
    	idx.start <- segs[j,"chrIdxStart"]
    	idx.end   <- segs[j,"chrIdxEnd"]
    	x1 <- data$position[idx.start]/1e6
    	x2 <- data$position[idx.end]/1e6
    	y1 <- median(data$log2ratio[idx.start:idx.end], na.rm=TRUE)
        segments(x0=x1,y0=y1,x1=x2,y1=y1, col="red", lty=1, lwd=3)
        }
    abline(h=0, lty=2) 
    
    ## log2 tumor mBAF/normal mBAF
    par(mar=c(4,4,1,1))
    plot(data$position/1e6, data$log2mBAF, pch=20, col="gray", cex=cex,
         xlab="Position (Mb)",ylab="log2mBAF",
         main="") 
    for (j in 1:nrow(segs)) {
    	idx.start <- segs[j,"chrIdxStart"]
    	idx.end   <- segs[j,"chrIdxEnd"]
    	x1 <- data$position[idx.start]/1e6
    	x2 <- data$position[idx.end]/1e6
    	y1 <- median(data$log2mBAF[idx.start:idx.end], na.rm=TRUE)
        segments(x0=x1,y0=y1,x1=x2,y1=y1,col="red",lty=1,lwd=3)
        }
    abline(h=0, lty=2) 
    
    ## Tumor BAF
    par(mar=c(1,4,1,1))
    plot(data$position/1e6, data$tumor.BAF, pch=20, col="gray", cex=cex,
         xlab="Position (Mb)",ylab="Tumor BAF",main="") 
    abline(h=c(0,0.5,1), lty=2)
    
    ## Normal BAF
    par(mar=c(4,4,1,1))
    plot(data$position/1e6, data$normal.BAF, pch=20, col="gray", cex=cex,
         xlab="Position (Mb)",ylab="Normal BAF",main="") 
    abline(h=c(0,0.5,1), lty=2)
    par(oldpar)
}
